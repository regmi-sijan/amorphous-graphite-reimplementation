"""Small analytic tests; fixtures here are not simulated research results."""
import unittest
import numpy as np
from scipy.special import expit
from ase.build import bulk
from agworkflow.structure import graphite,identify_layers,in_plane_rdf,periodic_rings
from agworkflow.transport import fermi_level,transition_weights,conductivity_tensor,spc_row_density

class PeriodicStructureTests(unittest.TestCase):
    def test_covering_graph_rings_and_supercell(self):
        for rep in [(1,1,1),(2,1,1),(2,2,1)]:
            a=graphite(rep); rings=periodic_rings(a)
            self.assertEqual(rings['6'],len(a)//2)
            self.assertEqual(sum(rings.values()),len(a)//2)

    def test_in_plane_normalization_and_invariance(self):
        a=graphite(); layers=identify_layers(a)
        self.assertTrue(layers['detected'])
        r,g=in_plane_rdf(a,layers,cutoff=2,bins=100)
        self.assertAlmostEqual(r[np.argmax(g)],1.42,delta=.011)
        edges=np.linspace(0,2,101)
        sheet_density=len(a)*layers['spacing_A']/a.get_volume()
        self.assertAlmostEqual(float(np.sum(g*sheet_density*np.pi*np.diff(edges**2))),3)
        a.translate([.27,.16,1.23]);a.wrap();a.rotate(37,'x',rotate_cell=True)
        rotated=identify_layers(a)
        self.assertTrue(rotated['detected'])
        self.assertAlmostEqual(rotated['spacing_A'],layers['spacing_A'])
        np.testing.assert_allclose(in_plane_rdf(a,rotated,cutoff=2,bins=100)[1],g,atol=1e-12)

    def test_diamond_does_not_receive_sheet_rdf(self):
        a=bulk('C','diamond',a=3.567)
        layers=identify_layers(a)
        self.assertFalse(layers['detected'])
        with self.assertRaises(ValueError): in_plane_rdf(a,layers)

class TransportTests(unittest.TestCase):
    def test_electron_count(self):
        energies=np.array([[-2,-.1,.2,1],[-1.8,-.2,.3,1.1]])
        mu=fermi_level(energies,3,.05)
        self.assertAlmostEqual(2*expit((mu-energies)/.05).sum()/2,3,places=9)

    def test_degenerate_limit(self):
        e=np.array([0.,1e-11])
        w=transition_weights(e,0,.1,0,.05,10)
        self.assertTrue(np.all(np.isfinite(w)))
        self.assertGreater(w[0,0],0)
        self.assertAlmostEqual(w[0,1]/w[0,0],1,places=8)
        ee=np.array([0.,1e-7])
        self.assertAlmostEqual(transition_weights(ee,0,.1,0,.05,10)[0,1]/w[0,0],1,places=7)

    def test_tensor_gauge_rotation_and_positive_eigenvalues(self):
        rng=np.random.default_rng(73)
        z=rng.normal(size=(3,4,4))+1j*rng.normal(size=(3,4,4))
        p=(z+z.conj().transpose(0,2,1))/2
        w=transition_weights(np.array([-.3,-.1,.1,.4]),0,.05,.1,.2,20)
        tensor=conductivity_tensor(p,w)
        self.assertGreaterEqual(np.linalg.eigvalsh(tensor).min(),-1e-9)
        phase=np.exp(1j*rng.uniform(0,6,size=4))
        pp=p*phase.conj()[None,:,None]*phase[None,None,:]
        np.testing.assert_allclose(conductivity_tensor(pp,w),tensor,rtol=1e-12,atol=1e-8)
        q,_=np.linalg.qr(rng.normal(size=(3,3)))
        np.testing.assert_allclose(conductivity_tensor(np.einsum('ab,bij->aij',q,p),w),q@tensor@q.T,rtol=1e-12,atol=1e-8)

    def test_spc_against_dense_kernel_and_analytic_plane_waves(self):
        # Exact quadrature for a unitary mixture of periodic plane waves.
        n=32; length=2*np.pi; dv=length/n; x=np.arange(n)*dv
        modes=np.array([-1,0,1]); wave=np.exp(1j*x[:,None]*modes)/np.sqrt(length)
        rng=np.random.default_rng(4)
        u,_=np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))
        psi=wave@u; grad=np.zeros((3,n,3),complex)
        grad[0]=(wave*(1j*modes))@u
        p=np.zeros((3,3,3),complex);p[0]=u.conj().T@np.diag(modes)@u
        w=transition_weights(np.array([-.2,0,.25]),0,.08,.1,.12,length)
        row=spc_row_density(psi,grad,p,w)
        xi=psi.conj()[:,:,None]*(-1j*grad[0,:,None,:])
        kernel=np.einsum('ij,gij,hij->gh',w,xi,xi.conj())/3
        np.testing.assert_allclose(row,kernel.sum(axis=1)*dv,rtol=1e-12,atol=1e-8)
        sigma=np.trace(conductivity_tensor(p,w))/3
        self.assertAlmostEqual((row.sum()*dv).real/sigma,1,places=12)
        self.assertLess(abs((row.sum()*dv).imag)/sigma,1e-12)
        self.assertGreaterEqual(np.abs(row).sum()*dv,sigma*(1-1e-12))

if __name__=='__main__': unittest.main()
