import numpy as np
import pytest

from grain_growth_pf.mechanics.anisotropy import BoundaryLaw, LADDER, disorientation
from grain_growth_pf.mechanics.cahn_hoffman import variation, herring_force
from grain_growth_pf.entities.ordered_network import order_network
from grain_growth_pf.entities.conforming_network import triangle_interfaces, reconstruct_periodic


def test_symmetries_and_derivatives():
    law = BoundaryLaw()
    theta = np.linspace(-2, 2, 47)
    expected = law.evaluate(theta, .17, .39)
    for other in [law.evaluate(theta+1.31, .17+1.31, .39+1.31),
                  law.evaluate(theta+np.pi, .39, .17),
                  law.evaluate(theta, .17+np.pi/2, .39)]:
        for a, b in zip(expected, other):
            np.testing.assert_allclose(a, b, rtol=1e-12, atol=1e-12)
    eps = 1e-4
    f = lambda x: law.evaluate(x, .17, .39)[0]
    first = (f(theta-2*eps)-8*f(theta-eps)+8*f(theta+eps)-f(theta+2*eps))/(12*eps)
    second = (-f(theta+2*eps)+16*f(theta+eps)-30*f(theta)+16*f(theta-eps)-f(theta-2*eps))/(12*eps**2)
    np.testing.assert_allclose(first, expected[1], atol=1e-9)
    np.testing.assert_allclose(second, expected[2], atol=1e-7)
    assert np.all(expected[3] > 0)
    assert disorientation(0, np.pi/2) == 0


@pytest.mark.parametrize('closed', [False, True])
def test_exact_energy_gradient_and_force_closure(closed):
    theta = np.linspace(.1, 5.7, 16)
    points = np.column_stack((2*np.cos(theta), 1.4*np.sin(theta)))
    law = BoundaryLaw()
    result = variation(points, law, .1, .7, closed=closed)
    total = result.integrated_forces + result.endpoint_forces
    for eps in (1e-4, 1e-5, 1e-6):
        measured = np.zeros_like(points)
        for k in range(len(points)):
            for axis in (0, 1):
                plus, minus = points.copy(), points.copy()
                plus[k, axis] += eps
                minus[k, axis] -= eps
                measured[k, axis] = -(variation(plus, law, .1, .7, closed=closed).energy
                                      - variation(minus, law, .1, .7, closed=closed).energy)/(2*eps)
        np.testing.assert_allclose(total, measured, rtol=1e-4, atol=1e-7)
    np.testing.assert_allclose(total.sum(axis=0), 0, atol=1e-14)
    expected = np.zeros(2) if closed else result.line_vectors[-1]-result.line_vectors[0]
    np.testing.assert_allclose(result.integrated_forces.sum(axis=0), expected, atol=1e-14)


def test_isotropic_circle_and_straight_boundary():
    law = BoundaryLaw(LADDER['A0_ISOTROPIC'])
    theta = np.arange(64)*2*np.pi/64
    points = 3*np.column_stack((np.cos(theta), np.sin(theta)))
    result = variation(points, law, .13, .6, closed=True)
    np.testing.assert_allclose(result.normal_pressure, -1/3, atol=1e-13)
    np.testing.assert_allclose(result.integrated_forces.sum(axis=0), 0, atol=1e-14)
    straight = variation([[0, 0], [1, 2], [2, 4]], BoundaryLaw(), .13, .6)
    np.testing.assert_allclose(straight.normal_pressure, 0, atol=1e-14)
    g, dg, ddg, s, m = law.evaluate(theta, .13, .6)
    for actual, expected in [(g, 1), (dg, 0), (ddg, 0), (s, 1), (m, 4)]:
        np.testing.assert_allclose(actual, expected, atol=1e-14)


def test_junction_gradient_and_isotropic_herring():
    angles = np.arange(3)*2*np.pi/3
    arms = np.column_stack((np.cos(angles), np.sin(angles)))
    pairs = [(0, 1), (1, 2), (0, 2)]
    orientations = [.1, .3, .7]
    iso = BoundaryLaw(LADDER['A0_ISOTROPIC'])
    np.testing.assert_allclose(herring_force(arms, pairs, orientations, iso), 0, atol=1e-14)
    law = BoundaryLaw()
    force = herring_force(arms, pairs, orientations, law)
    def energy(j):
        return sum(variation([j, arm], law, orientations[i], orientations[k]).energy
                   for arm, (i, k) in zip(arms, pairs))
    eps = 1e-6
    for axis in range(2):
        direction = np.eye(2)[axis]*eps
        derivative = -(energy(direction)-energy(-direction))/(2*eps)
        assert force[axis] == pytest.approx(derivative, abs=1e-8)


def test_periodic_ordering_disconnected_pairs_and_split_energy():
    vertices = [[9, 0], [0, 0], [1, 0], [4, 1], [4, 2], [4, 3]]
    edges = [[1, 2], [4, 5], [0, 1], [3, 4]]
    paths = order_network(vertices, edges, [[0, 1]]*4, box=[10, 10])
    assert len(paths) == 2
    np.testing.assert_allclose(paths[0].points, [[9, 0], [10, 0], [11, 0]])
    whole = variation(paths[0].points, BoundaryLaw(), .1, .7)
    pieces = [variation(paths[0].points[k:k+2], BoundaryLaw(), .1, .7) for k in (0, 1)]
    assert sum(p.energy for p in pieces) == pytest.approx(whole.energy)
    np.testing.assert_allclose(pieces[0].endpoint_forces[-1]+pieces[1].endpoint_forces[0], whole.integrated_forces[1])
    with pytest.raises(ValueError, match='branch'):
        order_network([[0, 0], [1, 0], [0, 1], [-1, 0]], [[0, 1], [0, 2], [0, 3]], [[0, 1]]*3)


def test_periodic_winding_loop():
    paths = order_network([[0, 0], [2, 0], [4, 0], [6, 0], [8, 0]],
                         [[0, 1], [1, 2], [2, 3], [3, 4], [4, 0]], [[0, 1]]*5, box=[10, 10])
    path = paths[0]
    assert path.closed
    np.testing.assert_allclose(path.winding, [10, 0])
    v = variation(path.points, BoundaryLaw(), 0., .3, closed=True, closing_displacement=path.winding)
    np.testing.assert_allclose(v.normal_pressure, 0, atol=1e-14)


def test_invalid_geometry_fails_closed():
    with pytest.raises(ValueError, match='zero-length'):
        variation([[0, 0], [0, 0]], BoundaryLaw(), 0, 0)
    with pytest.raises(ValueError, match='ambiguous'):
        order_network([[0, 0], [5, 0]], [[0, 1]], [[0, 1]], box=[10, 10])
    with pytest.raises(ValueError, match='nonempty'):
        BoundaryLaw().normalize([], [], [], [])


def test_conforming_three_phase_junction():
    interfaces = triangle_interfaces(np.eye(3))
    assert len(interfaces) == 3
    for i, j, a, b in interfaces:
        assert min(np.linalg.norm(a-np.ones(3)/3), np.linalg.norm(b-np.ones(3)/3)) < 1e-14
        assert a[i] == pytest.approx(a[j])
        assert b[i] == pytest.approx(b[j])
    # A fourth phase can dominate in the interior without winning at any vertex.
    four = triangle_interfaces(np.vstack((np.eye(3), [.4, .4, .4])))
    assert len(four) == 6


def test_conforming_periodic_inclusion_and_label_exchange():
    y, x = np.mgrid[:12, :12]
    radius = np.hypot((x+6)%12-6, (y+6)%12-6)
    phase = np.clip(.5+(2.7-radius)/3, 0, 1)
    eta = np.stack((phase, 1-phase))
    vertices, edges, pairs, paths = reconstruct_periodic(eta)
    assert len(paths) == 1 and paths[0].closed
    assert len(edges) == len(paths[0].points)
    other = reconstruct_periodic(eta[::-1])
    def energy(ps):
        return sum(variation(p.points, BoundaryLaw(), .1, .7, closed=p.closed,
                             closing_displacement=p.winding if p.closed else None).energy for p in ps)
    assert energy(paths) == pytest.approx(energy(other[3]))


def test_loop_attached_to_junction_keeps_endpoint_forces():
    path = order_network([[0, 0], [1, 0], [1, 1], [0, 1]],
                         [[0, 1], [1, 2], [2, 3], [3, 0]], [[0, 1]]*4,
                         junctions=[0])[0]
    assert not path.closed
    assert path.endpoint_junctions == (0, 0)
    assert path.node_ids[0] == path.node_ids[-1]
    result = variation(path.points, BoundaryLaw(), .1, .7)
    np.testing.assert_allclose((result.integrated_forces+result.endpoint_forces).sum(axis=0), 0, atol=1e-14)


def test_wulff_weighted_curvature_static_mesh_convergence():
    law = BoundaryLaw()
    errors = []
    for n in (128, 256, 512):
        theta = np.arange(n)*2*np.pi/n
        points = law.vectors(theta, .1, .4)[0]
        result = variation(points, law, .1, .4, closed=True)
        errors.append(np.max(np.abs(result.normal_pressure+1)))
    assert errors[2] < .3*errors[1] < .1*errors[0]
    assert errors[-1] < .02
