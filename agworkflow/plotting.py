"""Plots from computed data only; no reference curves synthesized."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.ndimage import map_coordinates


def structural_plots(folder):
    folder=Path(folder)
    summary=json.loads((folder/'summary.json').read_text())
    fig,axes=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
    rdf=np.genfromtxt(folder/'rdf_3d.csv',delimiter=',',names=True)
    axes[0].plot(rdf['r_A'],rdf['g_3d'],label='3D RDF')
    if (folder/'rdf_in_plane.csv').exists():
        rdf=np.genfromtxt(folder/'rdf_in_plane.csv',delimiter=',',names=True)
        axes[0].plot(rdf['r_A'],rdf['g_2d'],label='In-plane 2D RDF')
    axes[0].set(xlabel='Distance (Å)',ylabel='g(r)');axes[0].legend()
    rings=summary['rings_chordless_periodic']
    axes[1].bar(list(rings),list(rings.values()))
    axes[1].set(xlabel='Ring size',ylabel='Periodic chordless ring count')
    fig.suptitle('Layers detected' if summary['layers']['detected'] else 'No validated graphitic layers')
    fig.savefig(folder/'structure.png',dpi=150);plt.close(fig)


def evolution_plot(folder):
    folder=Path(folder)
    data=json.loads((folder/'structural_time_series.json').read_text())
    times=[row['time_ps'] for row in data]
    fig,axes=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
    for key,label in [('threefold_fraction','Threefold coordination'),
                      ('planar_threefold_fraction_proxy','Planar threefold proxy')]:
        axes[0].plot(times,[row[key] for row in data],label=label)
    for key,label in [('layering_amplitude_proxy','Density modulation'),
                      ('normal_alignment_proxy','Local normal alignment')]:
        axes[1].plot(times,[row[key] for row in data],label=label)
    for ax in axes:
        ax.set(xlabel='Time (ps)',ylim=(-.05,1.05));ax.legend()
    fig.suptitle('Independent structural proxies; not the original paper order parameter')
    fig.savefig(folder/'structural_evolution.png',dpi=150);plt.close(fig)
    thermo=np.atleast_1d(np.genfromtxt(folder/'thermo.csv',delimiter=',',names=True))
    config=json.loads((folder/'config.json').read_text())
    fig,axes=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
    axes[0].plot(thermo['time_ps'],thermo['potential_eV_per_atom'])
    axes[0].set(xlabel='Time (ps)',ylabel='GAP potential energy (eV/atom)')
    axes[1].plot(thermo['time_ps'],thermo['temperature_K'])
    axes[1].axhline(config['temperature_K'],ls='--',color='gray',label='Thermostat target')
    axes[1].set(xlabel='Time (ps)',ylabel='Temperature (K)');axes[1].legend()
    fig.savefig(folder/'thermodynamics.png',dpi=150);plt.close(fig)


def electronic_plots(folder,layers):
    folder=Path(folder)
    fields=np.load(folder/'fields.npz')
    rho=fields['density_e_bohr3']; lattice=fields['lattice_bohr']
    # Density maps are planes through the periodic volumetric result.
    if layers['detected']:
        normal=np.array(layers['normal'])
        spacing=layers['spacing_A']/.529177210903
        origin=layers['phase_offset']*spacing*normal
        label='Detected neighboring plane family; periodic images included'
    else:
        normal=np.cross(lattice[0],lattice[1]);normal/=np.linalg.norm(normal)
        spacing=float(lattice[2]@normal);origin=np.zeros(3)
        label='Cell cross-sections; no graphitic layers detected'
    axis=np.eye(3)[np.argmin(abs(normal))]
    u=np.cross(normal,axis);u/=np.linalg.norm(u);v=np.cross(normal,u)
    span=max(np.linalg.norm(lattice,axis=1))
    x,y=np.meshgrid(np.linspace(0,span,120,endpoint=False),np.linspace(0,span,120,endpoint=False))
    fig,axes=plt.subplots(1,5,figsize=(14,3.2),constrained_layout=True)
    slices=[]
    for fraction in np.linspace(0,1,5):
        points=origin+x[...,None]*u+y[...,None]*v+fraction*spacing*normal
        frac=points@np.linalg.inv(lattice)
        grid=(frac%1)*np.array(rho.shape)
        slices.append(map_coordinates(rho,grid.reshape(-1,3).T,order=1,mode='grid-wrap').reshape(x.shape))
    vmax=max(float(np.quantile(s,.99)) for s in slices)
    for index,(ax,s) in enumerate(zip(axes,slices)):
        im=ax.imshow(s,origin='lower',vmin=0,vmax=vmax,cmap='viridis')
        ax.set_title(f'Plane {index+1}');ax.set_xticks([]);ax.set_yticks([])
    fig.colorbar(im,ax=axes,label='Electron density (e/bohr³)',shrink=.8)
    fig.suptitle(label)
    fig.savefig(folder/'density_slices.png',dpi=140);plt.close(fig)
    dos=np.genfromtxt(folder/'dos.csv',delimiter=',',names=True)
    fig,ax=plt.subplots(figsize=(6,4),constrained_layout=True)
    ax.plot(dos['energy_minus_mu_eV'],dos['states_per_eV_cell'])
    ax.axvline(0,color='gray',ls='--');ax.set(xlabel='E − μ (eV)',ylabel='DOS (states/eV/cell)')
    fig.savefig(folder/'dos.png',dpi=150);plt.close(fig)
    orbitals=fields['orbitals_e_bohr3'];bands=fields['orbital_band_indices']
    fig,axes=plt.subplots(1,4,figsize=(12,3.4),constrained_layout=True)
    for ax,data,band in zip(axes,orbitals,bands):
        im=ax.imshow(data.mean(axis=2).T,origin='lower',cmap='magma')
        ax.set_title(f'k₀ band {band}');ax.set_xticks([]);ax.set_yticks([])
        fig.colorbar(im,ax=ax,shrink=.7)
    fig.suptitle('Orbital |ψ|² averaged over third cell coordinate; no automatic π assignment')
    fig.savefig(folder/'orbitals.png',dpi=140);plt.close(fig)
    if (folder/'conductivity.npz').exists():
        c=np.load(folder/'conductivity.npz');iso=np.trace(c['tensor_S_m'],axis1=2,axis2=3)/3
        fig,ax=plt.subplots(figsize=(6,4),constrained_layout=True)
        for eta,row in zip(c['broadening_eV'],iso):ax.plot(c['photon_eV'],row,label=f'η={eta:g} eV')
        ax.set(xlabel='Photon energy (eV)',ylabel='Isotropic σ (S/m)');ax.legend()
        fig.savefig(folder/'conductivity.png',dpi=150);plt.close(fig)
        spc=np.load(folder/'spc.npz')['row_density']
        fig,axes=plt.subplots(1,2,figsize=(9,4),constrained_layout=True)
        for ax,data,title in zip(axes,[abs(spc).mean(axis=2),spc.real.mean(axis=2)],['Positive SPC','Signed conductivity contribution']):
            im=ax.imshow(data.T,origin='lower',cmap='magma' if title=='Positive SPC' else 'coolwarm')
            ax.set_title(title);fig.colorbar(im,ax=ax,shrink=.75)
        fig.suptitle('Projection along third cell coordinate; (S/m)/bohr³')
        fig.savefig(folder/'spc.png',dpi=150);plt.close(fig)
