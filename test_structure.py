"""Physics/geometry checks on reference structures; no fitted expected outputs."""
import unittest
import numpy as np
from ase import Atoms
from ase.build import bulk, graphene
from reproduce import random_carbon, geometry, layering, rdf_3d, ring_counts, AMU_A3_TO_G_CM3


class StructuralChecks(unittest.TestCase):
    def test_density_and_periodic_exclusion(self):
        atoms = random_carbon(80, 2.44, 42)
        self.assertAlmostEqual(sum(atoms.get_masses()) / atoms.get_volume() * AMU_A3_TO_G_CM3, 2.44, places=8)
        d = atoms.get_all_distances(mic=True)
        np.fill_diagonal(d, np.inf)
        self.assertGreaterEqual(d.min(), 1.2)

    def test_graphene_coordination_and_rings(self):
        atoms = graphene(a=2.46, size=(5, 5, 1), vacuum=5)
        atoms.pbc = True
        m = geometry(atoms)
        self.assertEqual(m['threefold_fraction'], 1)
        self.assertEqual(m['planar_threefold_fraction_proxy'], 1)
        self.assertAlmostEqual(m['normal_alignment_proxy'], 1)
        rings = ring_counts(atoms)
        self.assertEqual(rings['6'], len(atoms) // 2)
        self.assertEqual(sum(rings.values()), rings['6'])

    def test_diamond_is_not_planar_threefold(self):
        atoms = bulk('C', 'diamond', a=3.567, cubic=True).repeat((2,2,2))
        m = geometry(atoms)
        self.assertEqual(m['fourfold_fraction'], 1)
        self.assertEqual(m['planar_threefold_fraction_proxy'], 0)

    def test_periodic_translation_and_rotation(self):
        atoms = random_carbon(80, 2.44, 42)
        before = geometry(atoms), layering(atoms), ring_counts(atoms)
        atoms.translate([3.2, -5.3, 9.1]); atoms.wrap()
        atoms.rotate(33, 'x', rotate_cell=True)
        after = geometry(atoms), layering(atoms), ring_counts(atoms)
        for key in before[0]:
            self.assertAlmostEqual(before[0][key], after[0][key], places=10)
        self.assertAlmostEqual(before[1]['layering_amplitude_proxy'], after[1]['layering_amplitude_proxy'], places=10)
        self.assertEqual(before[2], after[2])

    def test_known_layer_spacing(self):
        rng = np.random.default_rng(1)
        positions = np.column_stack([rng.uniform(0, 12, 200), rng.uniform(0, 12, 200), np.repeat(np.arange(4)*3.1, 50)])
        atoms = Atoms('C'*200, positions=positions, cell=[12,12,12.4], pbc=True)
        m = layering(atoms)
        self.assertAlmostEqual(m['layering_amplitude_proxy'], 1)
        self.assertAlmostEqual(m['candidate_spacing_A'], 3.1)

    def test_rdf_ideal_gas_normalization(self):
        rng = np.random.default_rng(5)
        atoms = Atoms('C'*1500, positions=rng.uniform(0,20,(1500,3)), cell=[20]*3, pbc=True)
        r,g = rdf_3d(atoms, bins=30)
        self.assertLess(abs(g[r>2].mean()-1), 0.04)

    def test_boundary_crossing_hexagon(self):
        angles = np.arange(6)*np.pi/3
        p = np.column_stack([1.42*np.cos(angles),1.42*np.sin(angles),np.zeros(6)])
        atoms = Atoms('C6', positions=p, cell=[12]*3, pbc=True)
        atoms.wrap()
        self.assertEqual(ring_counts(atoms)['6'], 1)


if __name__ == '__main__':
    unittest.main()
