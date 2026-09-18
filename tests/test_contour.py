import numpy as np
import pytest

from pivnp.contour import (
    NeighborAverageCorrection,
    NoContourCorrection,
    make_contour_correction,
)
from pivnp.state import Nodes

# Malla de 4x3 nodos (3x2 celdas). "#" = sin datos.
#   fila 2:  #  #  #  #
#   fila 1:  1  2  #  #
#   fila 0:  3  4  5  #
MASK = np.array([[0, 0, 0, 1], [0, 0, 1, 1], [1, 1, 1, 1]], dtype=np.int8).ravel()
VX = np.array([[3, 4, 5, 0], [1, 2, 0, 0], [0, 0, 0, 0]], dtype=float).ravel()


def _nodes():
    nodes = Nodes.zeros(12, 12)
    nodes.velocity[:, 0] = VX
    nodes.is_nan[:] = MASK
    return nodes


def test_no_correction_changes_nothing():
    nodes = _nodes()
    NoContourCorrection().apply(nodes, 3, 2)
    np.testing.assert_array_equal(nodes.velocity[:, 0], VX)
    assert not nodes.filled.any()


def test_fills_boundary_nodes_with_neighbor_average():
    nodes = _nodes()
    NeighborAverageCorrection(min_neighbors=3, layers=1).apply(nodes, 3, 2)
    vx = nodes.velocity[:, 0].reshape(3, 4)
    # nodo (1,2): vecinos válidos 4, 5, 2 -> media 11/3
    assert vx[1, 2] == pytest.approx(11 / 3)
    # nodo (2,1): vecinos válidos 1 y 2 (y (1,0)...) -> solo 2 -> no llega al mínimo de 3
    assert nodes.filled.reshape(3, 4)[1, 2]
    assert not nodes.filled.reshape(3, 4)[2, 3]
    # las marcas de datos medidos no cambian: no se activan partículas en el aire
    np.testing.assert_array_equal(nodes.is_nan, MASK)


def test_more_layers_fill_further():
    one, two = _nodes(), _nodes()
    NeighborAverageCorrection(min_neighbors=2, layers=1).apply(one, 3, 2)
    NeighborAverageCorrection(min_neighbors=2, layers=3).apply(two, 3, 2)
    assert two.filled.sum() > one.filled.sum()
    measured = MASK == 0
    np.testing.assert_array_equal(two.velocity[measured, 0], VX[measured])


def test_factory_and_validation():
    assert isinstance(make_contour_correction(0), NoContourCorrection)
    assert make_contour_correction(1, 4, 2) == NeighborAverageCorrection(4, 2)
    with pytest.raises(ValueError):
        NeighborAverageCorrection(min_neighbors=0)
    with pytest.raises(ValueError):
        NeighborAverageCorrection(layers=0)
