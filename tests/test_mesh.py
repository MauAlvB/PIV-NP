import math

import numpy as np
import pytest

from pivnp.config import parse_par
from pivnp.mesh import Grid, pivnp_grid, staggered_grid

from .legacy_reference import mesh_arrays, ucelda
from .test_config import PAR

NC, NFIL, AXC, AYC = 12, 3, 0.212115, 0.18


@pytest.fixture(params=[False, True], ids=["pivnp", "staggered"])
def grids(request):
    staggered = request.param
    nch = NC // NFIL
    if staggered:
        grid = Grid(nch + 1, NFIL + 1, AXC, AYC, -AXC / 2, -AYC / 2)
    else:
        grid = Grid(nch, NFIL, AXC, AYC, 0.0, 0.0)
    return grid, mesh_arrays(NC, NFIL, AXC, AYC, staggered)


def test_row_coordinates_match_legacy(grids):
    grid, arrays = grids
    assert grid.row_y.tolist() == arrays[6][1:]


def test_cell_nodes_match_legacy(grids):
    grid, arrays = grids
    for cell in range(grid.n_cells):
        x0 = grid.x0 + (cell % grid.n_cols + 0.5) * grid.dx
        y0 = grid.row_y[cell // grid.n_cols] + grid.dy / 2
        indc, nodes = ucelda(x0, y0, arrays, AXC)
        assert indc == cell + 1
        assert tuple(grid.cell_nodes(cell) + 1) == nodes


def _candidate_points(grid: Grid, rng) -> np.ndarray:
    xs = np.concatenate([
        rng.uniform(grid.x0 - grid.dx, grid.x0 + (grid.n_cols + 1) * grid.dx, 3000),
        grid.x0 + np.arange(-1, grid.n_cols + 2) * grid.dx,  # bordes exactos de columna
        np.nextafter(grid.x0 + np.arange(grid.n_cols + 1) * grid.dx, -np.inf),
    ])
    ys = np.concatenate([
        rng.uniform(grid.y0 - grid.dy, grid.row_y[-1] + grid.dy, xs.size - grid.n_rows - 1),
        grid.row_y,
    ])
    pts = np.stack([xs, rng.permutation(ys)], axis=1)
    extra = [[math.nan, 0.1], [0.1, math.nan], [grid.x0, grid.y0], [np.inf, 0.0]]
    return np.vstack([pts, extra])


def test_locate_matches_legacy_linear_search(grids):
    grid, arrays = grids
    points = _candidate_points(grid, np.random.default_rng(1))
    cells = grid.locate(points)
    for (x, y), cell in zip(points, cells, strict=True):
        found = ucelda(x, y, arrays, AXC)
        assert cell == (found[0] - 1 if found else -1), (x, y)


def test_grids_from_config():
    cfg = parse_par(PAR)
    g1, g2 = pivnp_grid(cfg), staggered_grid(cfg)
    assert (g1.n_cols, g1.n_rows, g1.n_nodes) == (59, 34, 2100)
    assert (g2.n_cols, g2.n_rows, g2.n_cells) == (60, 35, 2100)  # one cell per PIVlab point
    assert g2.x0 == -cfg.cell_width / 2 and g2.y0 == -cfg.cell_height / 2
