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
import os
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .config import CaseConfig, ConfigError, find_file, load_case
from .contour import ContourCorrection, make_contour_correction
from .gid_writer import GidWriter
from .mesh import Grid, particle_grid
from .nodal import (
    accumulate_particle_mass,
    compute_nodal_momentum_v1,
    compute_nodal_momentum_v2,
    load_measurements,
)
from .particles import cell_centers, create_particles
from .pivlab_io import FrameSource, pivlab_to_node
from .restart import RestartData, read_restart, write_restart
from .solver import advance_particles, output_mask, update_strains
from .state import Nodes, Particles

log = logging.getLogger("pivnp")


@dataclass(frozen=True)
class RunOptions:
    prefetch: int = 4  # archivos PIVlab leídos por adelantado
    contour_min_neighbors: int = 3
    contour_layers: int = 1
    eol: str = os.linesep  # fin de línea de los archivos GiD (CRLF en Windows, como Intel)
    log_every: int = 10  # cada cuántos pasos se informa del archivo analizado
    #: Reproduce los errores del Fortran original (H-01, H-04, H-08 y H-13) para poder
    #: repetir análisis antiguos. Ver docs/HALLAZGOS.md.
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
        self.particles: Particles = create_particles(config, self.grid)
        self.nodes = Nodes.zeros(config.n_nodes, self.grid.n_nodes)
        self.point_to_node = pivlab_to_node(config.n_cols, config.n_rows)
        # Centros de celda de la malla desplazada (H-13: el original solo los genera si no
        # es un reinicio, y entonces todas las velocidades acaban en la celda 1).
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
        contour = make_contour_correction(config.contour, options.contour_min_neighbors,
                                          options.contour_layers)
        return cls(case_dir, name, config, frames, contour, options)

    @property
    def restart_path(self) -> Path:
        return self.case_dir / f"{self.case_name}.REC"

    def run(self) -> RunSummary:
        cfg = self.config
        started = time.perf_counter()
        summary = RunSummary(self.case_name, cfg.n_particles, cfg.total_steps)
        log.info("LEYENDO DATOS... caso %s: %d partículas, %d pasos", self.case_name,
                 cfg.n_particles, cfg.total_steps)

        if cfg.restart:
            self._load_restart()

        step, t = 0, 0.0
        with GidWriter(self.case_dir, self.case_name, self.options.eol) as writer:
            if cfg.restart:
                self._write_output(writer, step, t, summary)

            steps = range(1, cfg.total_steps + 1)
            for frame in self.frames.frames(steps):
                step += 1
                t += cfg.dt
                self._update_nodes(frame, step)
                advance_particles(self.particles, self.nodes, self.grid, cfg, step)
                update_strains(self.particles, self.nodes, self.grid, cfg, step,
                               self.options.legacy_compat)
                if step == 1 or step % cfg.print_every == 0:
                    self._write_output(writer, step, t, summary)

        self._save_restart()
        summary.elapsed_s = time.perf_counter() - started
        log.info("ANALYSIS FINISHED en %.1f s", summary.elapsed_s)
        return summary

    # --- pasos ---------------------------------------------------------------------------
    def _update_nodes(self, frame, step: int) -> None:
        cfg = self.config
        accumulate_particle_mass(self.particles, self.grid, self.nodes, step,
                                 accumulate=cfg.mesh_version == 2)
        load_measurements(frame, self.point_to_node, self.nodes, self.options.legacy_compat)
        self.contour.apply(self.nodes, cfg.n_cols, cfg.n_rows)
        if cfg.mesh_version == 1:
            compute_nodal_momentum_v1(self.nodes)
        else:
            compute_nodal_momentum_v2(self.nodes, self.grid, self.centers, self.particles, step)
        if step == 1 or step % self.options.log_every == 0:
            log.info("ANALYZING %s", frame.source.name)

    def _write_output(self, writer: GidWriter, step: int, t: float,
                      summary: RunSummary) -> None:
        located = output_mask(self.particles, self.grid, step)
        if not summary.output_times:  # primer instante impreso: malla y cabecera
            # H-08: la malla GiD debe tener las posiciones iniciales, porque los
            # desplazamientos se acumulan desde ellas y GiD dibuja malla + desplazamiento.
            # El original escribía las posiciones ya movidas por el primer paso.
            positions = self.particles.position
            if not self.options.legacy_compat:
                positions = positions - self.particles.displacement
            writer.write_mesh(positions, self.particles.nan_initial)
            writer.start_results()
        writer.write_step(t, self.particles, self.nodes, located, self.config.moisture)
        summary.output_times.append(t)

    # --- reinicio ------------------------------------------------------------------------
    def _load_restart(self) -> None:
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
        p.displacement[:n] = data.displacement
        p.step_displacement[:n] = data.displacement
        p.strain[:n] = data.strain
        p.eq_strain[:n] = data.eq_strain
        p.nan_initial[:n] = data.nan_initial

    def _save_restart(self) -> None:
        n = self.config.n_particles
        p = self.particles
        write_restart(self.restart_path, RestartData(
            mesh_version=self.config.mesh_version,
            position=p.position[:n], displacement=p.displacement[:n], strain=p.strain[:n],
            eq_strain=p.eq_strain[:n], nan_initial=p.nan_initial[:n],
        ))


def run_case(case_dir: Path, case_name: str | None = None,
             options: RunOptions = DEFAULT_OPTIONS) -> RunSummary:
    """Ejecuta el caso descrito por ``PIV-NP.TXT`` (o ``case_name``) en ``case_dir``."""
    return Simulation.from_directory(case_dir, case_name, options).run()
