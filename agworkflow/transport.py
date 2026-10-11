"""Broadened independent-particle Kubo-Greenwood and the SPC row-sum field.

Atomic units inside this module; exported conductivity is S/m. Restricted
orbitals have per-spin occupations f (0..1); the prefactor contains spin factor 2.
The Greenwood divided-difference form replaces 1/omega by 1/(ej-ei) on the exact
delta support. With finite broadening this is an explicit regularization. It has
a finite degenerate/DC limit and includes diagonal/intraband matrix elements.
DC values are finite-broadening estimates, not a converged transport prediction.
"""
import numpy as np
from scipy.special import expit
from scipy.optimize import brentq
from scipy.constants import elementary_charge, hbar, physical_constants

HA_EV=physical_constants['Hartree energy in eV'][0]
SIGMA_AU=elementary_charge**2/(hbar*physical_constants['Bohr radius'][0])


def fermi_level(energies, nelectrons, sigma_Ha):
    energies=np.asarray(energies)
    nk,nb=energies.shape
    if sigma_Ha<=0 or not 0<nelectrons<2*nb:
        raise ValueError('Need positive smearing and both occupied and empty states')
    def count(mu): return 2*np.sum(expit((mu-energies)/sigma_Ha))/nk-nelectrons
    return brentq(count,float(energies.min()-100*sigma_Ha-1),float(energies.max()+100*sigma_Ha+1))


def transition_weights(energy,mu,sigma_Ha,photon_Ha,eta_Ha,volume_bohr3,nk=1):
    """Positive weights for all ordered (i,j), including the degenerate limit."""
    if eta_Ha<=0 or sigma_Ha<=0 or photon_Ha<0 or volume_bohr3<=0:
        raise ValueError('Invalid transport units/parameters')
    f=expit((mu-np.asarray(energy))/sigma_Ha)
    delta=energy[None,:]-energy[:,None]
    diff=f[:,None]-f[None,:]
    degenerate=abs(delta)<1e-9
    divided=np.divide(diff,delta,out=np.zeros_like(delta),where=~degenerate)
    mid=expit((mu-(energy[:,None]+energy[None,:])/2)/sigma_Ha)
    divided[degenerate]=(mid*(1-mid)/sigma_Ha)[degenerate]
    gaussian=np.exp(-.5*((delta-photon_Ha)/eta_Ha)**2)/(np.sqrt(2*np.pi)*eta_Ha)
    return (2*np.pi/(volume_bohr3*nk))*divided*gaussian*SIGMA_AU


def conductivity_tensor(momentum,weights):
    """p[a,i,j]=<i|p_a|j>; symmetric dissipative tensor in S/m."""
    return np.einsum('ij,aij,bij->ab',weights,momentum,momentum.conj()).real


def spc_row_density(psi,gradient,momentum,weights):
    """Complex row sum per unit volume for isotropic sigma (S/m)/bohr^3.

    xi_ij,a = psi_i* (-i d_a psi_j). Summing xi*p_ij* avoids materializing the
    grid-by-grid conduction matrix. abs(row) is the authors' positive SPC field;
    only the integral of the complex/signed row sum recovers conductivity.
    """
    row=np.zeros(len(psi),complex)
    for a in range(3):
        row+=np.einsum('gi,ij,gj->g',psi.conj(),weights*momentum[a].conj(),
                       -1j*gradient[a],optimize=True)/3
    return row


def spectra(energies,momenta,mu,sigma_Ha,photon_eV,broadening_eV,volume):
    result=np.zeros((len(broadening_eV),len(photon_eV),3,3))
    for b,eta in enumerate(broadening_eV):
        for w,photon in enumerate(photon_eV):
            for ek,pk in zip(energies,momenta):
                weight=transition_weights(ek,mu,sigma_Ha,photon/HA_EV,eta/HA_EV,volume,len(energies))
                result[b,w]+=conductivity_tensor(pk,weight)
    return result
