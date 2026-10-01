"""Orchestration of a full analysis (formerly the main program of PIV-NP).

Flow of every step (the same as the original):

1. Relocate particles and accumulate their nodal mass   (VELOCIDADES)
2. Load the measurements of the step into the nodes     (VELOCIDADES)
3. Correct the velocities at the boundary               (CONTOUR)
4. Nodal momentum                                       (VELOCIDADES)
5. Interpolate velocity/displacement to the particles   (SOLMOV)
6. Strains, energies and new position                   (SOLMOV)
7. Write results when due                               (IMPRES_GiD)
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
    AVERAGE,
    AVERAGE_2023,
    NO_AVERAGE,
    accumulate_particle_mass,
    compute_nodal_momentum_v1,
    compute_nodal_momentum_v2,
    load_measurements,
)
from .particles import cell_centers, create_particles
from .pivlab_io import pivlab_to_node
from .restart import NodalState, RestartData, read_restart, write_restart
from .solver import advance_particles, count_nan_nodes, output_mask, update_strains
from .sources import DisplacementSource, build_source
from .state import Nodes, Particles

log = logging.getLogger("pivnp")


@dataclass(frozen=True)
class RunOptions:
    #: Where the displacement data comes from; see :mod:`pivnp.sources`.
    source: str = "pivlab"
    prefetch: int = 4  # steps read ahead
    contour_min_neighbors: int = 3  # ICONTOUR=1: neighbours with data needed
    contour_layers: int = 1  # ICONTOUR=1 and 3: layers of points to rebuild
    contour_min_particles: int = 1  # ICONTOUR=2: particles needed around the point
    eol: str = os.linesep  # line ending of the GiD files (CRLF on Windows, like Intel)
    log_every: int = 10  # how often the file being analysed is reported
    #: Reproduce exactly the behaviour of the original Fortran, to repeat with this version
    #: analyses that were made with it.
    legacy_compat: bool = False
    #: With IVERSION=2, reproduce the distribution of the 2023 version: every staggered-grid
    #: node took the mean of the points contributing to it, instead of the sum of quarters.
    #: Needed to repeat the analyses made with that version.
    legacy_2023_average: bool = False


DEFAULT_OPTIONS = RunOptions()


@dataclass
class RunSummary:
    case_name: str
    n_particles: int
    steps: int
    output_times: list[float] = field(default_factory=list)
    elapsed_s: float = 0.0


class Simulation:
    """One PIV-NP analysis over the files of a case directory."""

    def __init__(self, case_dir: Path, case_name: str, config: CaseConfig,
                 frames: DisplacementSource, contour: ContourCorrection,
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
        # Cell centres of the staggered grid. The original only generated them when it was
        # not a restart, so in compatibility mode they are left at zero.
        legacy_restart = config.restart and options.legacy_compat
        if config.mesh_version == 2 and not legacy_restart:
            self.centers = cell_centers(self.grid)
        else:
            self.centers = np.zeros((self.grid.n_cells, 2))

    @classmethod
    def from_directory(cls, case_dir: Path, case_name: str | None = None,
                       options: RunOptions = DEFAULT_OPTIONS) -> Simulation:
        name, config = load_case(case_dir, case_name)
        frames = build_source(options.source, case_dir, config, options.prefetch)
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

    @property
    def staggered_average(self) -> int:
        """What to do with the boundary nodes of the staggered grid (IVERSION=2 only)."""
        if self.options.legacy_2023_average:
            return AVERAGE_2023
        return AVERAGE if self.contour.normalizes_staggered else NO_AVERAGE

    def check_frame_interval(self) -> float | None:
        """Warn when the DT of the ``.PAR`` does not match the interval of the source.

        If they do not match, the displacements come out multiplied by the ratio between the
        two: a mistake with no symptom other than results at a different scale. A source
        that does not know its interval returns ``None`` and nothing is checked.
        """
        interval = self.frames.frame_interval()
        if interval is None or math.isclose(interval, self.config.dt, rel_tol=1e-3):
            return interval
        log.warning(
            "The .PAR uses DT=%g s, but the PIVlab files were exported with an interval "
            "between images of %g s: displacements will come out multiplied by %.4g. Check "
            "the DT of the .PAR or the interval you exported the data with.",
            self.config.dt, interval, self.config.dt / interval)
        return interval

    def run(self) -> RunSummary:
        cfg = self.config
        started = time.perf_counter()
        summary = RunSummary(self.case_name, cfg.n_particles, cfg.total_steps)
        log.info("READING DATA... case %s: %d particles, %d steps", self.case_name,
                 cfg.n_particles, cfg.total_steps)
        self.check_frame_interval()

        step, t, resumed = 0, 0.0, False
        if cfg.restart:
            step, t, resumed = self._load_restart()
            if resumed:
                log.info("Continuing from step %d (t = %g s)", step, t)

        writer = GidWriter(self.case_dir, self.case_name, self.options.eol,
                           self.options.legacy_compat)
        with writer:
            if cfg.restart and not resumed:
                # With a .REC from the original Fortran there is no way to know how far the
                # analysis had got: a new series of results is started, as the original did.
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
        if self.frames.images is not None:
            log.info("Moisture: %s", self.frames.images.quality_summary())
        summary.elapsed_s = time.perf_counter() - started
        log.info("ANALYSIS FINISHED in %.1f s", summary.elapsed_s)
        return summary

    # --- steps ---------------------------------------------------------------------------
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
                                      normalize=self.staggered_average)
        if step == 1 or step % self.options.log_every == 0:
            log.info("ANALYZING %s", frame.label)

    def _write_output(self, writer: GidWriter, step: int, t: float, summary: RunSummary,
                      append: bool = False) -> None:
        located = output_mask(self.particles, self.grid, step)
        count_nan_nodes(self.particles, self.nodes, self.grid, self.config.mesh_version, step)
        if not summary.output_times:  # first printed step: mesh and header
            # The GiD mesh carries the initial positions, because the displacements are
            # accumulated from them and GiD draws mesh + displacement. The original wrote the
            # positions already moved by the first step.
            positions = (self.particles.position if self.options.legacy_compat
                         else self.particles.initial_position)
            writer.write_mesh(positions, self.particles.nan_initial)
            # When continuing an analysis, results are appended to the previous ones instead
            # of overwriting them.
            writer.start_results(append=append)
        writer.write_step(t, self.particles, self.nodes, located, self.config.moisture)
        summary.output_times.append(t)

    # --- restart -------------------------------------------------------------------------
    def _load_restart(self) -> tuple[int, float, bool]:
        """Load the ``.REC``; returns (step, time, ``True`` when the series continues).

        A ``.REC`` carrying nodal state (written by this version) makes it possible to
        continue the analysis exactly where it stopped. With one from the original Fortran,
        or in compatibility mode, it starts at step 0 and time 0, as the original did.
        """
        cfg = self.config
        data = read_restart(find_file(self.case_dir, self.restart_path.name))
        if data.n_particles != cfg.n_particles:
            raise ConfigError(f"{self.restart_path.name}: it holds {data.n_particles} "
                              f"particles and the .PAR defines {cfg.n_particles}")
        if data.mesh_version != cfg.mesh_version:
            raise ConfigError(f"{self.restart_path.name}: IVERSION={data.mesh_version} "
                              f"differs from the .PAR ({cfg.mesh_version})")
        n = data.n_particles
        p = self.particles
        p.position[:n] = data.position
        # A .REC from the original did not store the initial position: it is rebuilt.
        p.initial_position[:n] = (data.position - data.displacement
                                  if data.initial_position is None else data.initial_position)
        p.displacement[:n] = data.displacement
        p.strain[:n] = data.strain
        p.eq_strain[:n] = data.eq_strain
        p.nan_initial[:n] = data.nan_initial
        if self.options.legacy_compat:
            # The original carried the accumulated displacement into the "instantaneous" one.
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


def moisture_source(case_dir: Path, case_name: str, frames: DisplacementSource):
    """Moisture source computed from the test images (MOISTER=2).

    The PIV-NP grid is placed on the image with the coordinates and the conversion factor the
    source reports, and with the mask of points that have data.
    """
    from .moisture.source import Mesh, source_for_case

    x, y, metres_per_pixel, has_data = frames.mesh_in_metres(1)
    return source_for_case(case_dir, case_name, Mesh(x, y, metres_per_pixel), has_data)


def write_moisture_files(case_dir: Path, case_name: str | None = None,
                         options: RunOptions = DEFAULT_OPTIONS) -> int:
    """Compute the moisture from the images and write it to ``Moist_<n>.TXT``.

    It does not run the analysis: it is meant for looking at the moisture fields on their own,
    or for comparing them with those of earlier analyses. Returns how many files were written.
    """
    from .moisture.source import Mesh, write_moist

    name, config = load_case(case_dir, case_name)
    case_dir = Path(case_dir)
    frames = build_source(options.source, case_dir, config, options.prefetch)
    frames.moisture = False  # it is about to be computed from the images, not read
    x, y, metres_per_pixel, _ = frames.mesh_in_metres(1)
    mesh = Mesh(x, y, metres_per_pixel)
    frames.images = moisture_source(case_dir, name, frames)

    written = 0
    for frame in frames.frames(range(1, config.total_steps + 1)):
        target = case_dir / f"Moist_{frame.step}.TXT"
        write_moist(target, mesh, frame.moisture, frame.saturation)
        written += 1
        if frame.step == 1 or frame.step % options.log_every == 0:
            log.info("WRITTEN %s", target.name)
    log.info("%d moisture files written in %s", written, case_dir)
    log.info("Moisture: %s", frames.images.quality_summary())
    return written


def run_case(case_dir: Path, case_name: str | None = None,
             options: RunOptions = DEFAULT_OPTIONS) -> RunSummary:
    """Run the case described by ``PIV-NP.TXT`` (or by ``case_name``) in ``case_dir``."""
    return Simulation.from_directory(case_dir, case_name, options).run()
