import numpy as np
import pytest

from pivnp.config import parse_par
from pivnp.contour import (
    ContourContext,
    ExtrapolationCorrection,
    NeighborAverageCorrection,
    NoContourCorrection,
    ParticleAverageCorrection,
    make_contour_correction,
)
from pivnp.mesh import Grid, particle_grid
from pivnp.nodal import compute_nodal_momentum_v2
from pivnp.particles import cell_centers, create_particles
from pivnp.state import Nodes, Particles

# Grid of 4x3 PIVlab points (3x2 cells). "#" = no data.
#   row 2:  #  #  #  #
#   row 1:  1  2  #  #
#   row 0:  3  4  5  #
MASK = np.array([[0, 0, 0, 1], [0, 0, 1, 1], [1, 1, 1, 1]], dtype=np.int8).ravel()
VX = np.array([[3, 4, 5, 0], [1, 2, 0, 0], [0, 0, 0, 0]], dtype=float).ravel()

PAR = """case
b2
6 12 2 2 1.0 1.0
b3
0.5 4 1 0 {version} 1 {contour} 0 0
b4
2000 0.4
"""


def _context(version=1, contour=1, mask=MASK, vx=VX, particles=None):
    cfg = parse_par(PAR.format(version=version, contour=contour))
    grid = particle_grid(cfg)
    nodes = Nodes.zeros(cfg.n_nodes, grid.n_nodes)
    nodes.velocity[:, 0] = vx
    nodes.is_nan[:] = mask
    if particles is None:
        particles = Particles.zeros(cfg.n_particles, n_lost=cfg.n_nodes)
    return ContourContext(nodes, grid, particles, cfg)


def test_no_correction_changes_nothing():
    ctx = _context()
    NoContourCorrection().apply(ctx)
    np.testing.assert_array_equal(ctx.nodes.velocity[:, 0], VX)
    assert not ctx.nodes.filled.any()


# --- the property measured on real cases -------------------------------------------------
def test_rebuilding_a_hidden_node_beats_leaving_the_zero():
    """Hide a node that has data and check whether the rebuild beats leaving the zero.

    This is the test that was run on the real cases, here in miniature and on a smooth field:
    hide nodes that did have a measurement and compare against what PIVlab measured. Across
    the five cases tried, the neighbour average won every time; the particle average was the
    weakest.
    """
    rows, columns = 5, 6
    row, column = np.mgrid[0:rows, 0:columns]
    field = (2.0 + 0.5 * column + 0.25 * row).ravel()  # smooth field, no noise
    par = PAR.replace("6 12 2 2 1.0 1.0", f"{(rows - 1) * (columns - 1)} {rows * columns} "
                                          f"2 {rows - 1} 1.0 1.0")
    cfg = parse_par(par.format(version=1, contour=1))

    hidden = rows * columns // 2 + 1
    for method, tolerance in ((NeighborAverageCorrection(min_neighbors=1), 0.35),
                              (ExtrapolationCorrection(), 0.35)):
        grid = particle_grid(cfg)
        nodes = Nodes.zeros(cfg.n_nodes, grid.n_nodes)
        nodes.velocity[:, 0] = field
        nodes.velocity[hidden] = 0.0
        nodes.is_nan[hidden] = 1
        method.apply(ContourContext(nodes, grid, Particles.zeros(cfg.n_particles,
                                                                 n_lost=cfg.n_nodes), cfg))
        error = abs(nodes.velocity[hidden, 0] - field[hidden])
        assert nodes.filled[hidden], method
        assert error < tolerance * abs(field[hidden]), method  # better than leaving the zero


# --- ICONTOUR=1: neighbour average ------------------------------------------------------------
def test_neighbor_average_fills_boundary_points():
    ctx = _context()
    NeighborAverageCorrection(min_neighbors=3, layers=1).apply(ctx)
    vx = ctx.nodes.velocity[:, 0].reshape(3, 4)
    assert vx[1, 2] == pytest.approx(11 / 3)  # neighbours with data: 4, 5 and 2
    assert ctx.nodes.filled.reshape(3, 4)[1, 2]
    assert not ctx.nodes.filled.reshape(3, 4)[2, 3]  # not enough neighbours
    np.testing.assert_array_equal(ctx.nodes.is_nan, MASK)  # what was measured does not change


def test_more_layers_fill_further():
    one, two = _context(), _context()
    NeighborAverageCorrection(min_neighbors=2, layers=1).apply(one)
    NeighborAverageCorrection(min_neighbors=2, layers=3).apply(two)
    assert two.nodes.filled.sum() > one.nodes.filled.sum()
    measured = MASK == 0
    np.testing.assert_array_equal(two.nodes.velocity[measured, 0], VX[measured])


# --- ICONTOUR=2: average of the surrounding particles -----------------------------------------
def _particles_with_velocity(cfg, grid, vx):
    particles = create_particles(cfg, grid)
    particles.velocity[:, 0] = vx
    return particles


def test_particle_average_uses_surrounding_particles():
    cfg = parse_par(PAR.format(version=1, contour=2))
    grid = particle_grid(cfg)
    particles = _particles_with_velocity(cfg, grid, 7.0)
    ctx = _context(contour=2, particles=particles)
    ParticleAverageCorrection().apply(ctx)
    missing = ctx.nodes.is_nan != 0
    # every point without data that has particles around takes their velocity
    assert ctx.nodes.filled.sum() == missing.sum()
    np.testing.assert_allclose(ctx.nodes.velocity[missing, 0], 7.0)
    np.testing.assert_array_equal(ctx.nodes.velocity[~missing, 0], VX[~missing])


def test_particle_average_needs_enough_particles():
    cfg = parse_par(PAR.format(version=1, contour=2))
    grid = particle_grid(cfg)
    particles = _particles_with_velocity(cfg, grid, 7.0)
    particles.lost[:] = True  # no particle available
    ctx = _context(contour=2, particles=particles)
    ParticleAverageCorrection(min_particles=1).apply(ctx)
    assert not ctx.nodes.filled.any()
    np.testing.assert_array_equal(ctx.nodes.velocity[:, 0], VX)


def test_particle_average_in_staggered_mesh_uses_its_own_cell():
    cfg = parse_par(PAR.format(version=2, contour=2))
    grid = particle_grid(cfg)
    particles = create_particles(cfg, grid)
    cells = grid.locate(particles.position)
    particles.velocity[:, 0] = cells  # every particle carries the number of its cell
    ctx = _context(version=2, contour=2, particles=particles)
    ParticleAverageCorrection().apply(ctx)
    missing = np.flatnonzero(ctx.nodes.is_nan != 0)
    # point i is the centre of cell i, so it gets exactly its own number
    np.testing.assert_allclose(ctx.nodes.velocity[missing, 0], missing)


# --- ICONTOUR=3: extrapolation ----------------------------------------------------------------
def test_extrapolation_continues_a_linear_field():
    # v = 2·column on the two left columns; the rest without data
    mask = np.array([[0, 0, 1, 1]] * 3, dtype=np.int8).ravel()
    vx = np.array([[0.0, 2.0, 0.0, 0.0]] * 3).ravel()
    ctx = _context(mask=mask, vx=vx, contour=3)
    ExtrapolationCorrection(layers=2).apply(ctx)
    result = ctx.nodes.velocity[:, 0].reshape(3, 4)
    np.testing.assert_allclose(result[:, 2], 4.0)  # 2·2 − 0
    np.testing.assert_allclose(result[:, 3], 6.0)  # second layer
    assert ctx.nodes.filled.sum() == 6


def test_extrapolation_falls_back_to_the_average():
    ctx = _context(contour=3)
    ExtrapolationCorrection(layers=1).apply(ctx)
    assert ctx.nodes.filled.any()
    np.testing.assert_array_equal(ctx.nodes.is_nan, MASK)


# --- distribution on the staggered grid --------------------------------------------------
def test_staggered_normalization_only_changes_the_boundary():
    grid = Grid(3, 3, 1.0, 1.0, -0.5, -0.5)  # 9 PIVlab points -> 16 nodes
    particles = Particles.zeros(1, n_lost=9)
    plain, normalized = Nodes.zeros(9, grid.n_nodes), Nodes.zeros(9, grid.n_nodes)
    for nodes in (plain, normalized):
        nodes.velocity[:, 0] = 5.0
    centers = cell_centers(grid)
    compute_nodal_momentum_v2(plain, grid, centers, particles, step=2, normalize=False)
    compute_nodal_momentum_v2(normalized, grid, centers, particles, step=2, normalize=True)

    interior = normalized.active_count == 4
    np.testing.assert_array_equal(normalized.momentum[interior], plain.momentum[interior])
    np.testing.assert_allclose(normalized.momentum[interior, 0], 5.0)
    # at the boundary the original kept a fraction; now it is the mean
    boundary = (normalized.active_count > 0) & ~interior
    assert boundary.any()
    np.testing.assert_allclose(normalized.momentum[boundary, 0], 5.0)
    assert (plain.momentum[boundary, 0] < 5.0).all()


# --- selection -------------------------------------------------------------------------------
def test_factory_and_validation():
    assert isinstance(make_contour_correction(0), NoContourCorrection)
    assert make_contour_correction(1, 4, 2) == NeighborAverageCorrection(4, 2)
    assert make_contour_correction(2, min_particles=3) == ParticleAverageCorrection(3)
    assert make_contour_correction(3, layers=2) == ExtrapolationCorrection(2)
    assert not NoContourCorrection().normalizes_staggered
    assert NeighborAverageCorrection().normalizes_staggered
    with pytest.raises(ValueError):
        make_contour_correction(4)
    with pytest.raises(ValueError):
        NeighborAverageCorrection(min_neighbors=0)
    with pytest.raises(ValueError):
        ParticleAverageCorrection(min_particles=0)
    with pytest.raises(ValueError):
        ExtrapolationCorrection(layers=0)


def test_par_rejects_unknown_contour():
    with pytest.raises(ValueError, match="ICONTOUR"):
        parse_par(PAR.format(version=1, contour=9))
