"""Periodic layer geometry and correctly normalized in-plane RDF."""
import math
import numpy as np
from ase import Atoms
from ase.neighborlist import neighbor_list
from reproduce import geometry, layering, AMU_A3_TO_G_CM3


def graphite(repeat=(1, 1, 1), density=2.44):
    """AB graphite reference. Keep 1.42 A bonds; set c to requested density."""
    a = 1.42*np.sqrt(3)
    cell = [[a,0,0],[-a/2,a*np.sqrt(3)/2,0],[0,0,6.7]]
    atoms = Atoms('C4', scaled_positions=[[0,0,0],[1/3,2/3,0],
                 [0,0,.5],[2/3,1/3,.5]], cell=cell, pbc=True).repeat(repeat)
    target_volume = sum(atoms.get_masses())*AMU_A3_TO_G_CM3/density
    cell = atoms.cell.array.copy()
    cell[2] *= target_volume/atoms.get_volume()
    atoms.set_cell(cell, scale_atoms=True)
    return atoms


def identify_layers(atoms, cutoff=1.85, amplitude_min=.55, planar_fraction_min=.70,
                    rms_max_A=.30):
    """Test a commensurate plane family; thresholds are explicitly chosen here.

    A positive result requires density modulation AND locally planar bonding.
    Unlike clustering wrapped Cartesian z coordinates this also handles tilted
    layers. The number of plane phases is gcd(h,k,l), not an absolute sheet count.
    """
    m = layering(atoms)
    g = geometry(atoms, cutoff)
    result = {**m, **g, 'detected': False, 'method': 'reciprocal phase plus local geometry',
              'thresholds': {'amplitude_min': amplitude_min,
              'planar_fraction_min': planar_fraction_min, 'rms_max_A': rms_max_A}}
    if m['candidate_hkl'] is None:
        return result
    h = np.array(m['candidate_hkl'])
    q = h@atoms.cell.reciprocal()
    normal = q/np.linalg.norm(q)
    phase = atoms.get_scaled_positions()@h
    phase0 = np.angle(np.exp(2j*np.pi*phase).mean())/(2*np.pi)
    residual = (phase-phase0+.5)%1-.5
    rms = float(np.sqrt(np.mean(residual**2))*m['candidate_spacing_A'])
    detected = (m['layering_amplitude_proxy'] >= amplitude_min and
                g['planar_threefold_fraction_proxy'] >= planar_fraction_min and
                g['normal_alignment_proxy'] >= .70 and rms <= rms_max_A)
    result.update(detected=bool(detected), normal=normal.tolist(), plane_rms_A=rms,
                  phase_offset=float(phase0), periodic_phase_count=int(math.gcd(*h.tolist())))
    if detected:
        result['spacing_A'] = m['candidate_spacing_A']
    return result


def in_plane_rdf(atoms, layers, cutoff=4.0, bins=100, thickness_A=.7):
    """2D g(r) of same-plane pairs; refuses a structure without validated layers.

    Ordered pairs are normalized by N * (N*d/V) * annulus_area. The sheet number
    density N*d/V follows from the volume density and measured plane spacing.
    Neighbor vectors include periodic shifts; distances are projected in-plane.
    """
    if not layers['detected']:
        raise ValueError('No validated layers: an in-plane RDF is undefined')
    normal = np.array(layers['normal'])
    d = layers['spacing_A']
    if thickness_A >= d/2:
        raise ValueError('Layer thickness cutoff would mix adjacent planes')
    vectors = neighbor_list('D', atoms, np.sqrt(cutoff**2+thickness_A**2))
    z = vectors@normal
    radii = np.sqrt(np.maximum(0, np.sum(vectors**2,axis=1)-z*z))
    edges = np.linspace(0,cutoff,bins+1)
    count,_ = np.histogram(radii[(abs(z)<thickness_A)&(radii>1e-8)],edges)
    sigma = len(atoms)*d/atoms.get_volume()
    value = count/(len(atoms)*sigma*np.pi*np.diff(edges**2))
    return (edges[:-1]+edges[1:])/2,value


def periodic_rings(atoms, cutoff=1.85, max_size=8, visit_limit=5_000_000):
    """Chordless rings in the periodic covering graph, including small cells.

    Vertices contain an atom id AND lattice shift. Winding paths are not closed
    rings. Cycles are identified modulo translation and traversal direction.
    """
    i,j,s = neighbor_list('ijS',atoms,cutoff)
    adj=[[] for _ in atoms]
    for a,b,shift in zip(i,j,s): adj[a].append((int(b),tuple(int(v) for v in shift)))
    edge_sets=[set(x) for x in adj]
    def connected(a,b):
        return (b[0],tuple(b[k]-a[k] for k in range(1,4))) in edge_sets[a[0]]
    def canonical(path):
        variants=[]
        for reverse in (path,list(reversed(path))):
            for k in range(len(path)):
                rotated=reverse[k:]+reverse[:k]
                origin=rotated[0]
                variants.append(tuple((v[0],*(v[x]-origin[x] for x in range(1,4))) for v in rotated))
        return min(variants)
    rings=set(); visits=0
    def walk(path):
        nonlocal visits
        visits+=1
        if visits>visit_limit: raise RuntimeError('Ring search visit limit exceeded; increase visit_limit')
        last=path[-1]
        for atom,shift in adj[last[0]]:
            nxt=(atom,*(last[k+1]+shift[k] for k in range(3)))
            if nxt==path[0]:
                if len(path)>=3: rings.add(canonical(path))
                continue
            if nxt in path or len(path)>=max_size: continue
            if any(connected(nxt,v) for v in path[1:-1]): continue
            if len(path)>2 and connected(last,path[0]): continue
            walk(path+[nxt])
    for root in range(len(atoms)): walk([(root,0,0,0)])
    return {str(n):sum(len(r)==n for r in rings) for n in range(3,max_size+1)}
