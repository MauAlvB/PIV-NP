"""Orquestación de un análisis completo (antes el programa principal de PIV-NP).

Flujo de cada paso (igual que el original):

1. Relocalizar partículas y acumular su masa nodal     (VELOCIDADES)
2. Cargar las medidas del instante en los nodos         (VELOCIDADES)
3. Corregir velocidades en el contorno                  (CONTOUR)
4. Cantidad de movimiento nodal                         (VELOCIDADES)
5. Interpolar velocidad/desplazamiento a partículas     (SOLMOV)
6. Deformaciones, energías y nueva posición             (SOLMOV)
7. Escribir resultados si toca                          (IMPRES_GiD)
"""

from __future__ import annotations

import logging
import math
import os
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .config import CaseConfig, ConfigError, find_file, load_case
from .contour import ContourContext, ContourCorrection, make_contour_correction
from .gid_writer import GidWriter
from .mesh import Grid, particle_grid
from .nodal import (
    accumulate_particle_mass,
    compute_nodal_momentum_v1,
    compute_nodal_momentum_v2,
    load_measurements,
)
from .particles import cell_centers, create_particles
from .pivlab_io import FrameSource, frame_interval_in_header, pivlab_to_node
from .restart import NodalState, RestartData, read_restart, write_restart
from .solver import advance_particles, count_nan_nodes, output_mask, update_strains
from .state import Nodes, Particles

log = logging.getLogger("pivnp")


@dataclass(frozen=True)
class RunOptions:
    prefetch: int = 4  # archivos PIVlab leídos por adelantado
    contour_min_neighbors: int = 3  # ICONTOUR=1: vecinos con dato necesarios
    contour_layers: int = 1  # ICONTOUR=1 y 3: capas de puntos a reconstruir
    contour_min_particles: int = 1  # ICONTOUR=2: partículas necesarias alrededor
    eol: str = os.linesep  # fin de línea de los archivos GiD (CRLF en Windows, como Intel)
    log_every: int = 10  # cada cuántos pasos se informa del archivo analizado
    #: Reproduce exactamente el comportamiento del Fortran original, para repetir con
    #: esta versión análisis hechos con él.
    legacy_compat: bool = False


DEFAULT_OPTIONS = RunOptions()


@dataclass
class RunSummary:
    case_name: str
    n_particles: int
    steps: int
    output_times: list[float] = field(default_factory=list)
    elapsed_s: float = 0.0


class Simulation:
    """Un análisis PIV-NP sobre los archivos de un directorio de caso."""

    def __init__(self, case_dir: Path, case_name: str, config: CaseConfig,
                 frames: FrameSource, contour: ContourCorrection,
                 options: RunOptions = DEFAULT_OPTIONS) -> None:
        self.case_dir = Path(case_dir)
        self.case_name = case_name
        self.config = config
        self.frames = frames
        self.contour = contour
        self.options = options

        self.grid: Grid = particle_grid(config)
        self.particles: Particles = create_particles(config, self.grid, options.legacy_compat)
        self.nodes = Nodes.zeros(config.n_nodes, self.grid.n_nodes)
        self.point_to_node = pivlab_to_node(config.n_cols, config.n_rows)
        # Centros de celda de la malla desplazada. El original solo los generaba cuando
        # no era un reinicio, así que en modo compatibilidad se dejan a cero.
        legacy_restart = config.restart and options.legacy_compat
        if config.mesh_version == 2 and not legacy_restart:
            self.centers = cell_centers(self.grid)
        else:
            self.centers = np.zeros((self.grid.n_cells, 2))

    @classmethod
    def from_directory(cls, case_dir: Path, case_name: str | None = None,
                       options: RunOptions = DEFAULT_OPTIONS) -> Simulation:
        name, config = load_case(case_dir, case_name)
        frames = FrameSource(case_dir, config.n_nodes, config.pivlab_format,
                             config.moisture, prefetch=options.prefetch)
        if config.moisture_from_images:
            frames.images = moisture_source(case_dir, name, frames)
            frames.moisture = False
        contour = make_contour_correction(config.contour, options.contour_min_neighbors,
                                          options.contour_layers,
                                          options.contour_min_particles)
        return cls(case_dir, name, config, frames, contour, options)

    @property
    def restart_path(self) -> Path:
        return self.case_dir / f"{self.case_name}.REC"

    def check_frame_interval(self) -> float | None:
        """Avisa si el DT del ``.PAR`` no coincide con el intervalo que usó PIVlab.

        Si no coinciden, los desplazamientos salen multiplicados por el cociente entre
        ambos: es un error que no da ningún síntoma salvo resultados a otra escala.
        """
        try:
            intervalo = frame_interval_in_header(self.frames.velocity_path(1))
        except (FileNotFoundError, OSError):
            return None
        if intervalo is None or math.isclose(intervalo, self.config.dt, rel_tol=1e-3):
            return intervalo
        log.warning(
            "El .PAR usa DT=%g s, pero los archivos PIVlab se exportaron con un intervalo "
            "entre imágenes de %g s: los desplazamientos saldrán multiplicados por %.4g. "
            "Revisa el DT del .PAR o el intervalo con el que exportaste desde PIVlab "
            "(si los archivos no vienen de PIVlab, ignora este aviso).",
            self.config.dt, intervalo, self.config.dt / intervalo)
        return intervalo

    def run(self) -> RunSummary:
        cfg = self.config
        started = time.perf_counter()
        summary = RunSummary(self.case_name, cfg.n_particles, cfg.total_steps)
        log.info("LEYENDO DATOS... caso %s: %d partículas, %d pasos", self.case_name,
                 cfg.n_particles, cfg.total_steps)
        self.check_frame_interval()

        step, t, resumed = 0, 0.0, False
        if cfg.restart:
            step, t, resumed = self._load_restart()
            if resumed:
                log.info("Continuando desde el paso %d (t = %g s)", step, t)

        writer = GidWriter(self.case_dir, self.case_name, self.options.eol,
                           self.options.legacy_compat)
        with writer:
            if cfg.restart and not resumed:
                # Con un .REC del Fortran original no se sabe por dónde iba el análisis:
                # se empieza una serie de resultados nueva, como hacía el original.
                self._write_output(writer, step, t, summary)

            steps = range(1, cfg.total_steps + 1)
            for frame in self.frames.frames(steps):
                step += 1
                t += cfg.dt
                self._update_nodes(frame, step)
                advance_particles(self.particles, self.nodes, self.grid, cfg, step,
                                  self.options.legacy_compat)
                update_strains(self.particles, self.nodes, self.grid, cfg, step,
                               self.options.legacy_compat)
                if step == 1 or step % cfg.print_every == 0:
                    self._write_output(writer, step, t, summary, append=resumed)

        self._save_restart(step, t)
        summary.elapsed_s = time.perf_counter() - started
        log.info("ANALYSIS FINISHED en %.1f s", summary.elapsed_s)
        return summary

    # --- pasos ---------------------------------------------------------------------------
    def _update_nodes(self, frame, step: int) -> None:
        cfg = self.config
        accumulate_particle_mass(self.particles, self.grid, self.nodes, step,
                                 accumulate=cfg.mesh_version == 2)
        load_measurements(frame, self.point_to_node, self.nodes, self.options.legacy_compat)
        self.contour.apply(ContourContext(self.nodes, self.grid, self.particles, cfg))
        if cfg.mesh_version == 1:
            compute_nodal_momentum_v1(self.nodes)
        else:
            compute_nodal_momentum_v2(self.nodes, self.grid, self.centers, self.particles, step,
                                      normalize=self.contour.normalizes_staggered)
        if step == 1 or step % self.options.log_every == 0:
            log.info("ANALYZING %s", frame.source.name)

    def _write_output(self, writer: GidWriter, step: int, t: float, summary: RunSummary,
                      append: bool = False) -> None:
        located = output_mask(self.particles, self.grid, step)
        count_nan_nodes(self.particles, self.nodes, self.grid, self.config.mesh_version, step)
        if not summary.output_times:  # primer instante impreso: malla y cabecera
            # La malla GiD lleva las posiciones iniciales, porque los desplazamientos se
            # acumulan desde ellas y GiD dibuja malla + desplazamiento. El original
            # escribía las posiciones ya movidas por el primer paso.
            positions = (self.particles.position if self.options.legacy_compat
                         else self.particles.initial_position)
            writer.write_mesh(positions, self.particles.nan_initial)
            # Al continuar un análisis, los resultados se añaden a los anteriores en vez
            # de sobrescribirlos.
            writer.start_results(append=append)
        writer.write_step(t, self.particles, self.nodes, located, self.config.moisture)
        summary.output_times.append(t)

    # --- reinicio ------------------------------------------------------------------------
    def _load_restart(self) -> tuple[int, float, bool]:
        """Carga el ``.REC``; devuelve (paso, instante, ``True`` si se continúa la serie).

        Un ``.REC`` con estado nodal (escrito por esta versión) permite continuar el
        análisis exactamente donde se quedó. Con uno del Fortran original, o en modo
        compatibilidad, se empieza en el paso 0 y el instante 0, como hacía el original.
        """
        cfg = self.config
        data = read_restart(find_file(self.case_dir, self.restart_path.name))
        if data.n_particles != cfg.n_particles:
            raise ConfigError(f"{self.restart_path.name}: tiene {data.n_particles} partículas "
                              f"y el .PAR define {cfg.n_particles}")
        if data.mesh_version != cfg.mesh_version:
            raise ConfigError(f"{self.restart_path.name}: IVERSION={data.mesh_version} "
                              f"distinto del .PAR ({cfg.mesh_version})")
        n = data.n_particles
        p = self.particles
        p.position[:n] = data.position
        # Con un .REC del original no se guardó la posición inicial: se reconstruye.
        p.initial_position[:n] = (data.position - data.displacement
                                  if data.initial_position is None else data.initial_position)
        p.displacement[:n] = data.displacement
        p.strain[:n] = data.strain
        p.eq_strain[:n] = data.eq_strain
        p.nan_initial[:n] = data.nan_initial
        if self.options.legacy_compat:
            # El original arrastraba el desplazamiento acumulado al "instantáneo".
            p.step_displacement[:n] = data.displacement
            return 0, 0.0, False
        if data.nodes is None:
            return 0, 0.0, False
        data.nodes.apply_to(self.nodes)
        return data.step, data.time, True

    def _save_restart(self, step: int, t: float) -> None:
        n = self.config.n_particles
        p = self.particles
        write_restart(self.restart_path, RestartData(
            mesh_version=self.config.mesh_version,
            position=p.position[:n], displacement=p.displacement[:n], strain=p.strain[:n],
            eq_strain=p.eq_strain[:n], nan_initial=p.nan_initial[:n],
            time=t, step=step, nodes=NodalState.from_nodes(self.nodes),
            initial_position=p.initial_position[:n],
        ), extended=not self.options.legacy_compat)


def moisture_source(case_dir: Path, case_name: str, frames: FrameSource):
    """Fuente de humedad calculada desde las imágenes del ensayo (MOISTER=2).

    La malla de PIV-NP se sitúa en la imagen con las coordenadas y el factor de conversión
    del primer archivo PIVlab, y con la máscara de nodos que PIVlab midió.
    """
    from .humedad.fuente import Malla, fuente_de_caso

    x, y, metros_por_pixel, con_dato = frames.mesh_in_metres(1)
    return fuente_de_caso(case_dir, case_name, Malla(x, y, metros_por_pixel), con_dato)


def write_moisture_files(case_dir: Path, case_name: str | None = None,
                         options: RunOptions = DEFAULT_OPTIONS) -> int:
    """Calcula la humedad desde las imágenes y la escribe en ``Moist_<n>.TXT``.

    No ejecuta el análisis: sirve para revisar los campos de humedad por separado, o para
    compararlos con los de análisis anteriores. Devuelve cuántos archivos se han escrito.
    """
    from .humedad.fuente import Malla, escribir_moist

    name, config = load_case(case_dir, case_name)
    case_dir = Path(case_dir)
    frames = FrameSource(case_dir, config.n_nodes, config.pivlab_format,
                         prefetch=options.prefetch)
    x, y, metros_por_pixel, _ = frames.mesh_in_metres(1)
    malla = Malla(x, y, metros_por_pixel)
    frames.images = moisture_source(case_dir, name, frames)

    escritos = 0
    for frame in frames.frames(range(1, config.total_steps + 1)):
        destino = case_dir / f"Moist_{frame.step}.TXT"
        escribir_moist(destino, malla, frame.moisture, frame.saturation)
        escritos += 1
        if frame.step == 1 or frame.step % options.log_every == 0:
            log.info("ESCRITO %s", destino.name)
    log.info("%d archivos de humedad escritos en %s", escritos, case_dir)
    return escritos


def run_case(case_dir: Path, case_name: str | None = None,
             options: RunOptions = DEFAULT_OPTIONS) -> RunSummary:
    """Ejecuta el caso descrito por ``PIV-NP.TXT`` (o ``case_name``) en ``case_dir``."""
    return Simulation.from_directory(case_dir, case_name, options).run()
