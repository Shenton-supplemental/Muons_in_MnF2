from ase.io import read, write
from ase import Atoms
from ase.calculators.vasp import Vasp2


#import commands
import os
import sys


from shutil import copy2

import numpy as np
import re


from pymatgen.core import Element
from pymatgen.electronic_structure.core import Spin, OrbitalType
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.vasp import Vasprun

elementMn = Element('Mn')
elementF = Element('F')

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.gridspec import GridSpec
import matplotlib as mpl
from matplotlib.ticker import FormatStrFormatter
mpl.style.use('./plots_style.mplstyle')




def rgbline(ax, k, e, red, green, blue, alpha=1.):
    # creation of segments based on
    # http://nbviewer.ipython.org/urls/raw.github.com/dpsanders/matplotlib-examples/master/colorline.ipynb
    pts = np.array([k, e]).T.reshape(-1, 1, 2)
    seg = np.concatenate([pts[:-1], pts[1:]], axis=1)

    nseg = len(k) - 1
    r = [0.5 * (red[i] + red[i + 1]) for i in range(nseg)]
    g = [0.5 * (green[i] + green[i + 1]) for i in range(nseg)]
    b = [0.5 * (blue[i] + blue[i + 1]) for i in range(nseg)]
    a = np.ones(nseg, np.float) * alpha
    lc = LineCollection(seg, colors=list(zip(r, g, b, a)), linewidth=2)
    ax.add_collection(lc)


def get_mags(outcar, natoms = 81, ionic_step = -1):
    
    """Get the magnetic moments from a given ionic step"""
    
    re_mag  = re.compile("magnetization \(x\)")
    mags = []
    with open(outcar) as f:
        for line in f:
            if re_mag.search(line):
                next(f)
                next(f)
                next(f)
                for i in range(natoms):
                    mag = float(next(f).split()[-1])
                    mags.append(mag)
                next(f) # line of -------
                totmag = float(next(f).split()[-1])
    if ionic_step < 0:
        ionic_step -= 1 # because the final magmom is printed twice
    mags = np.array(mags).reshape(int(len(mags)/ natoms), natoms)
    mags = mags[ionic_step]
    return [mags, totmag]



def set_mags(atoms, magmoms):
    if len(atoms) == len(magmoms):
        for i, atom in enumerate(atoms):
            atom.magmom = magmoms[i]
    else:
        print('ERROR: the magmoms array should be the same length as the atoms object')




def get_hyperfine(outcar, natoms = 10, chosen_atom=-1, verbose=True):

    """Get the hyperfine tensor for each atom for the final ionic step.
    """

    re_totalspin  = re.compile("Total magnetic moment S")
    re_contact  = re.compile("Fermi contact \(isotropic\) hyperfine coupling parameter \(MHz\)")
    re_dipole   = re.compile("Dipolar hyperfine coupling parameters \(MHz\)")

    contacts = np.zeros(natoms)
    dipole_tensors = np.zeros((natoms, 3, 3))

    with open(outcar) as f:
        for line in f:
            if re_totalspin.search(line):
                totalS = float(line.split()[-1])

            if re_contact.search(line):
                next(f) # -------------------------------------------------------------
                next(f) #  ion      A_pw      A_1PS     A_1AE     A_1c      A_tot
                next(f) #  -------------------------------------------------------------

                for i in range(natoms):
                    # NB: core correction not included in A_tot
                    A_1c, A_tot  = np.array(next(f).split()).astype(float)[-2:] # core correction, A_tot
                    contacts[i] = A_tot + A_1c



            if re_dipole.search(line):
                next(f) # -------------------------------------------------------------
                next(f) #    ion      A_xx      A_yy      A_zz      A_xy      A_xz      A_yz
                next(f) #  -------------------------------------------------------------

                for i in range(natoms):
                    A_dip = np.array((next(f).split())).astype(float)

                    # diagonal:
                    dipole_tensors[i][0][0] = A_dip[1]
                    dipole_tensors[i][1][1] = A_dip[2]
                    dipole_tensors[i][2][2] = A_dip[3]

                    #off-diagonal
                    dipole_tensors[i][0][1] = A_dip[4]
                    dipole_tensors[i][1][0] = A_dip[4]

                    dipole_tensors[i][0][2] = A_dip[5]
                    dipole_tensors[i][2][0] = A_dip[5]

                    dipole_tensors[i][1][2] = A_dip[6]
                    dipole_tensors[i][2][1] = A_dip[6]



    total_hyperfine = dipole_tensors[chosen_atom] + np.eye(3)*contacts[chosen_atom]
    evals, evecs = np.linalg.eigh(total_hyperfine)
    
    if verbose:
        #total spin, fermi contact, dipolar tensors
        print('Fermi contact: {:6.3f} MHz'.format(contacts[chosen_atom]))
        print('Dipole tensor (MHz):')
        for row in dipole_tensors[chosen_atom]:
            print("{:10.3f} {:10.3f} {:10.3f}".format(*row))



        print('\nTotal hyperfine tensor: (MHz)')
        for row in total_hyperfine:
            print("{:10.3f} {:10.3f} {:10.3f}".format(*row))



        print('\nPrinciple values and directions of the hyperfine tensor:')
        for i in range(3):
            print("{:10.3f} Mhz    [{:6.3f} {:6.3f} {:6.3f}]".format(evals[i], *evecs[i]))


    return total_hyperfine
def extract_freqs(outcar_path, natoms, chosen_atom, gamma = 1):   
    # get the hyperfine tensor
    hf = get_hyperfine(outcar_path, natoms = natoms, chosen_atom=chosen_atom, verbose=False)

    # scale to muon precession frequency
    # (This is needed because I did not set NGYROMAG in the INCAR file.)
    hf *= gamma
    # high-field limit
    hf /= 2.0
    # split HF tensor into isotropic and anisotropic components
    iso = np.trace(hf) / 3.
    dipolar_tensor = hf - iso * np.eye(3)  
    # we want just the projections onto the z axis
    dip  =  np.linalg.norm(dipolar_tensor, axis=1)[2] * np.sign(dipolar_tensor[2][2])
    total = np.linalg.norm(hf, axis=1)[2] * np.sign(hf[2][2])

    return iso, dip, total

def get_max_force(atoms):
    # get the maximum force
    forces = atoms.get_forces()
    max_force = np.linalg.norm(forces, axis=1).max()
    return max_force


def plotbanddos(u, xc, dosruns_all, bandsruns_all, kpaths_all, sigma=0.02):
    dosruns   = dosruns_all[xc]
    bandsruns = bandsruns_all[xc]
    kpaths    = kpaths_all[xc]
    indU = int(round(u/2))
    print(indU)
    dosrun = dosruns[indU]
    
    chemsymbols = dosrun.atomic_symbols

    Mnindicies = [i for i, atom in enumerate(chemsymbols) if atom == "Mn"]
    Findicies  = [i for i, atom in enumerate(chemsymbols) if atom == "F"]
    
    pdos = dosrun.complete_dos.get_element_dos()
    print( 'Getting element projected DOS')

    # Get element projected DOS
    Mndos = pdos[elementMn].get_smeared_densities(sigma=sigma)
    odos  = pdos[elementF].get_smeared_densities(sigma=sigma)

    print( 'Importing bands run')
    kpath = kpaths[indU]
    bands   = bandsruns[indU].get_band_structure(kpath,
                                                 line_mode = False)
    
    
    # bands run path
    atoms = AseAtomsAdaptor.get_atoms(bandsruns[indU].final_structure)
    path = atoms.cell.bandpath(density=21)

    klinearcoord, ticks, labels = path.get_linear_kpoint_axis()

    pbands = bands.get_projection_on_elements()

    # set up matplotlib plot
    # ----------------------
    # set up 2 graph with aspec ration 2/1
    # plot 1: bands diagram
    # plot 2: Density of States
    gs = GridSpec(1, 2, width_ratios=[1.5, 1])
    fig = plt.figure()
    fig.suptitle(r"U = {0} eV".format(u))
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])  # , sharey=ax1)

    # set ylim for the plot
    emin = -8
    emax =  8
    ax1.set_ylim(emin, emax)
    ax2.set_ylim(emin, emax)

    # Band Diagram
    # ------------
    name = "Mn"
    
    # compute s, p, d normalized contributions
    contrib_up = np.zeros((bands.nb_bands, len(bands.kpoints), 3))
    contrib_dn = np.zeros((bands.nb_bands, len(bands.kpoints), 3))

    #----- spin up bands ---- #
    spin = Spin.up
    for b in range(bands.nb_bands):
        for k in range(len(bands.kpoints)):
            Mncontr = pbands[spin][b][k]["Mn"]**2
            Fcontr  = pbands[spin][b][k]["F"]**2
            tot = Mncontr + Fcontr
            if tot != 0.0:
                # Red
                contrib_up[b, k, 0] += Fcontr / tot
                #Green
                contrib_up[b, k, 1] += Mncontr / tot
                # Blue
                contrib_up[b, k, 2] = 0.0
    # plot bands using rgb mapping
    for b in range(bands.nb_bands):
        rgbline(ax1,
                klinearcoord,
                [e - bands.efermi for e in bands.bands[spin][b]],
                contrib_up[b, :, 0],
                contrib_up[b, :, 1],
                contrib_up[b, :, 2])

    #----- spin dn bands ---- #            
    spin = Spin.down            
    for b in range(bands.nb_bands):
        for k in range(len(bands.kpoints)):
            Mncontr = pbands[spin][b][k]["Mn"]**2
            Fcontr  = pbands[spin][b][k]["F"]**2
            tot = Mncontr + Fcontr
            if tot != 0.0:
                # Red
                contrib_dn[b, k, 0] += Fcontr / tot
                #Green
                contrib_dn[b, k, 1] = 0.0
                # Blue
                contrib_dn[b, k, 2] += Mncontr / tot

    # plot bands using rgb mapping
    for b in range(bands.nb_bands):
        rgbline(ax1,
                klinearcoord,
                [e - bands.efermi for e in bands.bands[spin][b]],
                contrib_dn[b, :, 0],
                contrib_dn[b, :, 1],
                contrib_dn[b, :, 2])

    # style
    ax1.set_xlabel("k-points")
    ax1.set_ylabel(r"$E - E_f$   /   eV")
    ax1.grid()

    # fermi level at 0
    ax1.hlines(y=0, xmin=0, xmax=len(bands.kpoints), color="k", lw=2)

    ax1.set_xticks(ticks)
    ax1.set_xticklabels(labels)

    ax1.set_xlim(0, klinearcoord[-1])

    # Density of states
    # -----------------

    ax2.set_yticklabels([])
    ax2.grid()
    ax2.set_xlim(1e-8, 15)
    ax2.set_xticklabels([])
    ax2.hlines(y=0, xmin=0, xmax=15, color='0.75', lw=2)
    ax2.set_xlabel("Density of States", labelpad=20)

    # spd contribution
    # --- Mn ---#
    ax2.plot(Mndos[Spin.up]+Mndos[Spin.down],
             dosrun.tdos.energies - dosrun.efermi,
             "b-", label=r"Mn", lw=2)

    # --- F ---#
    ax2.plot(odos[Spin.up]+odos[Spin.down],
             dosrun.tdos.energies - dosrun.efermi,
             "r-", label=r"F ", lw=2)

    # total dos
    ax2.fill_between(dosrun.tdos.get_smeared_densities(sigma=sigma)[Spin.up]
                    +dosrun.tdos.get_smeared_densities(sigma=sigma)[Spin.down],
                     0,
                     dosrun.tdos.energies - dosrun.efermi,
                     color=(0.7, 0.7, 0.7),
                     facecolor=(0.8, 0.8, 0.8))

    ax2.plot(dosrun.tdos.get_smeared_densities(sigma=sigma)[Spin.up]
            +dosrun.tdos.get_smeared_densities(sigma=sigma)[Spin.down],
             dosrun.tdos.energies - dosrun.efermi,
             color=(0.3, 0.3, 0.3),
             label="total DOS")

    # plot format style
    # -----------------
    ax2.legend(fancybox=True,
               shadow=True,
               ncol=2,
               loc=4) # 1-> upper right;4 => lower right
    plt.subplots_adjust(wspace=0)

    # plt.show()
    plt.savefig("images/{0}_banddos_Ueff{1}.png".format(xc, u),format='png', dpi=300)
    plt.close(fig)

def reorder_atoms(atoms_ref, atoms_to_reorder, check_species=True):
    """
    Takes in a reference ASE atoms object and another ASE atoms object to reorder
    
    Loops over reference atoms and, for each one, finds it's closest counterpart in atoms_to_reorder
    (taking into account the periodic boundaries via the minimum image convention).
    
    It finally reorders atoms_to_reorder such that they follow the order in the refernce one. 
    
    
    Returns the re-ordered atoms
    """
    atoms_to_reorder_new_indicies = []

    # loop over reference atoms
    for atom in atoms_ref:
        atoms_ref_index = atom.index

        # temporary copy of atoms_to_reorder
        temp = atoms_to_reorder.copy()
        
        # add atom to temp:
        temp.append(atom)

        # which atom is closest to the newly added atom in atoms_to_reorder
        distances = temp.get_distances(-1, range(len(atoms_to_reorder)), mic=True)

        
        if check_species:
            # loop over closest atoms to atom (in order of closeness)
            for closest_index in np.argsort(distances):
                if atom.symbol == atoms_to_reorder[closest_index].symbol:
                    # if the atom closest is of the same species, then
                    # we're done, otherwise keep going until you find 
                    # one that does match in species
                    atoms_to_reorder_index = closest_index
                    break
                else:
                    pass
            
        else:
            # get index of closest atom
            atoms_to_reorder_index = np.argsort(distances)[0]


        # append to array of new order
        atoms_to_reorder_new_indicies.append(atoms_to_reorder_index)

        
    # apply the new order
    return atoms_to_reorder[atoms_to_reorder_new_indicies]