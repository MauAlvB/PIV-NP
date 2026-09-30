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

# Malla de 4x3 puntos PIVlab (3x2 celdas). "#" = sin dato.
#   fila 2:  #  #  #  #
#   fila 1:  1  2  #  #
#   fila 0:  3  4  5  #
MASK = np.array([[0, 0, 0, 1], [0, 0, 1, 1], [1, 1, 1, 1]], dtype=np.int8).ravel()
VX = np.array([[3, 4, 5, 0], [1, 2, 0, 0], [0, 0, 0, 0]], dtype=float).ravel()

PAR = """caso
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


# --- la propiedad que se midió sobre casos reales ----------------------------------------
def test_rebuilding_a_hidden_node_beats_leaving_the_zero():
    """Se oculta un nodo con dato y se mira si la reconstrucción se acerca más que el cero.

    Es la prueba que se hizo sobre los casos reales, aquí en pequeño y con un campo suave:
    ocultar nodos que sí tenían medida y comparar con lo que medía PIVlab. Sobre los cinco
    casos probados, la media de vecinos ganó siempre; la media de partículas fue la peor.
    """
    filas, columnas = 5, 6
    fila, columna = np.mgrid[0:filas, 0:columnas]
    campo = (2.0 + 0.5 * columna + 0.25 * fila).ravel()  # campo suave, sin ruido
    par = PAR.replace("6 12 2 2 1.0 1.0", f"{(filas - 1) * (columnas - 1)} {filas * columnas} "
                                          f"2 {filas - 1} 1.0 1.0")
    cfg = parse_par(par.format(version=1, contour=1))

    oculto = filas * columnas // 2 + 1
    for metodo, tolerancia in ((NeighborAverageCorrection(min_neighbors=1), 0.35),
                               (ExtrapolationCorrection(), 0.35)):
        grid = particle_grid(cfg)
        nodes = Nodes.zeros(cfg.n_nodes, grid.n_nodes)
        nodes.velocity[:, 0] = campo
        nodes.velocity[oculto] = 0.0
        nodes.is_nan[oculto] = 1
        metodo.apply(ContourContext(nodes, grid, Particles.zeros(cfg.n_particles,
                                                                 n_lost=cfg.n_nodes), cfg))
        error = abs(nodes.velocity[oculto, 0] - campo[oculto])
        assert nodes.filled[oculto], metodo
        assert error < tolerancia * abs(campo[oculto]), metodo  # mejor que dejar el cero


# --- ICONTOUR=1: media de los vecinos ---------------------------------------------------------
def test_neighbor_average_fills_boundary_points():
    ctx = _context()
    NeighborAverageCorrection(min_neighbors=3, layers=1).apply(ctx)
    vx = ctx.nodes.velocity[:, 0].reshape(3, 4)
    assert vx[1, 2] == pytest.approx(11 / 3)  # vecinos con dato: 4, 5 y 2
    assert ctx.nodes.filled.reshape(3, 4)[1, 2]
    assert not ctx.nodes.filled.reshape(3, 4)[2, 3]  # sin vecinos suficientes
    np.testing.assert_array_equal(ctx.nodes.is_nan, MASK)  # no cambia qué se midió


def test_more_layers_fill_further():
    one, two = _context(), _context()
    NeighborAverageCorrection(min_neighbors=2, layers=1).apply(one)
    NeighborAverageCorrection(min_neighbors=2, layers=3).apply(two)
    assert two.nodes.filled.sum() > one.nodes.filled.sum()
    measured = MASK == 0
    np.testing.assert_array_equal(two.nodes.velocity[measured, 0], VX[measured])


# --- ICONTOUR=2: media de las partículas de alrededor -----------------------------------------
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
    # todos los puntos sin dato con partículas alrededor toman la velocidad de esas partículas
    assert ctx.nodes.filled.sum() == missing.sum()
    np.testing.assert_allclose(ctx.nodes.velocity[missing, 0], 7.0)
    np.testing.assert_array_equal(ctx.nodes.velocity[~missing, 0], VX[~missing])


def test_particle_average_needs_enough_particles():
    cfg = parse_par(PAR.format(version=1, contour=2))
    grid = particle_grid(cfg)
    particles = _particles_with_velocity(cfg, grid, 7.0)
    particles.lost[:] = True  # ninguna partícula disponible
    ctx = _context(contour=2, particles=particles)
    ParticleAverageCorrection(min_particles=1).apply(ctx)
    assert not ctx.nodes.filled.any()
    np.testing.assert_array_equal(ctx.nodes.velocity[:, 0], VX)


def test_particle_average_in_staggered_mesh_uses_its_own_cell():
    cfg = parse_par(PAR.format(version=2, contour=2))
    grid = particle_grid(cfg)
    particles = create_particles(cfg, grid)
    cells = grid.locate(particles.position)
    particles.velocity[:, 0] = cells  # cada partícula lleva el número de su celda
    ctx = _context(version=2, contour=2, particles=particles)
    ParticleAverageCorrection().apply(ctx)
    missing = np.flatnonzero(ctx.nodes.is_nan != 0)
    # el punto i es el centro de la celda i, así que recibe exactamente su número
    np.testing.assert_allclose(ctx.nodes.velocity[missing, 0], missing)


# --- ICONTOUR=3: extrapolación ----------------------------------------------------------------
def test_extrapolation_continues_a_linear_field():
    # v = 2·columna en las dos columnas de la izquierda; el resto sin dato
    mask = np.array([[0, 0, 1, 1]] * 3, dtype=np.int8).ravel()
    vx = np.array([[0.0, 2.0, 0.0, 0.0]] * 3).ravel()
    ctx = _context(mask=mask, vx=vx, contour=3)
    ExtrapolationCorrection(layers=2).apply(ctx)
    result = ctx.nodes.velocity[:, 0].reshape(3, 4)
    np.testing.assert_allclose(result[:, 2], 4.0)  # 2·2 − 0
    np.testing.assert_allclose(result[:, 3], 6.0)  # segunda capa
    assert ctx.nodes.filled.sum() == 6


def test_extrapolation_falls_back_to_the_average():
    ctx = _context(contour=3)
    ExtrapolationCorrection(layers=1).apply(ctx)
    assert ctx.nodes.filled.any()
    np.testing.assert_array_equal(ctx.nodes.is_nan, MASK)


# --- reparto en la malla desplazada -----------------------------------------------------
def test_staggered_normalization_only_changes_the_boundary():
    grid = Grid(3, 3, 1.0, 1.0, -0.5, -0.5)  # 9 puntos PIVlab -> 16 nodos
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
    # en el borde, el original se quedaba con una fracción; ahora es la media
    boundary = (normalized.active_count > 0) & ~interior
    assert boundary.any()
    np.testing.assert_allclose(normalized.momentum[boundary, 0], 5.0)
    assert (plain.momentum[boundary, 0] < 5.0).all()


# --- selección -------------------------------------------------------------------------------
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
