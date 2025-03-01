# disorder.py: Encoding site-disorder in a ".vasp" format
# Since the format ignores any term after the coordinate line (UNLESS selective dynamics is switched on),
# We can include atomic species and occupancy terms after them. The result is still readable.
# Note: do not use with selective dynamics!!!!

# 8 May 2024: I starts to implement the Atomic Sudoku code

from poshcar.supercell import *
from poshcar.atomsub import *
from poshcar.distance import *
from poshcar.additive import *

from pathlib import Path
import ase.io.cif as asecif
from chgnet.model.model import CHGNet
from chgnet.model import StructOptimizer
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.transformations.standard_transformations import OrderDisorderedStructureTransformation

def cif2vasp_occ(ciffile, verbose = True):
    # Reads in a .cif file, output .vasp file with fractional occupancy information which mimics a Selective Dynamics fileexit
    
    filename_header = ciffile[:-4] #output to same place as input
    cifdata = asecif.read_cif(ciffile, index=-1, fractional_occupancies = True)
    
    #print(cifdata.__dict__['info']); cif data extractions
    occdict = cifdata.__dict__.get('info').get('occupancy') # dictionary containing species and occupancies
    sites = cifdata.__dict__.get('arrays').get('spacegroup_kinds') # array of sites, in sequence
    positions = cifdata.__dict__.get('arrays').get('positions') # array of site coordinates

    # write a vasp file here with ase, then read from it (yeah man that's dumb I know)
    cifdata.write(filename_header+'.vasp', format = 'vasp')
    vaspdata = readfile(filename_header+'.vasp')
    vaspdata.insert(7, "Site-Disordered Structure\n") # it starts with 'S' so VESTA can read this (so lame)
    vaspdata = vaspdata[:10] # keep header + 1 atom
    vaspdata[5] = " Ph \n"
    vaspdata[6] = "  1 \n" # Placeholder atom "Ph"

    cursor = 9
    for i in occdict:
        if int(i) in sites: j = int(i)
        else: j = max([x for x in sites if x < int(i)], default=None) 
        k = np.where(sites == j)[0][0] # extract coordinates
        toadd = next(iter(occdict[i].items())) # extract atoms
        vaspdata = addatom(vaspdata, toadd[0], positions[k], verbose = False)
        cursor += 1
        vaspdata[cursor] = vaspdata[cursor].rstrip() + ls + (toadd[0]+str(j)).ljust(5) + "  " + str(toadd[1]) + "\n"
    
    vaspdata = removeatom(vaspdata, [1], verbose = False) # Remove the placeholder atom
    # write data to vasp file again
    writefile(vaspdata, filename_header+'.vasp', verbose)
    return vaspdata


def composition(data, verbose = True):
    # Return vector containing fractional 
    df = elemindices(data, verbose = False)
    total = df['Occupancy'].sum()
    sum_by_spp = df.groupby('Species')['Occupancy'].sum()
    proportions = sum_by_spp/total
    if verbose: display(proportions)
    return proportions.values


def VirtualLibrary(path, target, pauling_weight = 1, bond_threshold = 0.1, chgnet = None, structure_opt = None, verbose = True):
    # Source: a folder (path) already populated by virtual cells from supercell software
    # No further filling required!!
    # Check for target composition (target is compulsory)
    # header = folder name or path
    
    idlist   = [] # header, for reference
    datalist = [] # init array
    complist = [] # composition
    cdlist   = [] # distance from target composition (Frobenius)
    bondlist = [] # bonding matrices
    bdlist1  = []
    bdlist2  = []
    bdlist   = [] # distance from target bonding profile
    E_list   = [] # chgnet total energies

    # Load CHGNET bc we defo need to use this
    if chgnet == None:
        chgnet = CHGNet.load()
    if structure_opt == True: # if we want to relax the structures, but don't pass in the relaxer
        relaxer = StructOptimizer()

    # Take care of the target file (could be source file input to supercell)
    # or maybe even a user-defined cell for correlated disorder
    bme_targ = matrix_bonding_average(target, 'element', bond_threshold, verbose = verbose)[1]
    original_comp = composition(target, verbose = verbose) # register real compositional data

    # Iterate over all .cif files in folder
    path = Path(path)  # Ensure the path is a Path object
    for ciffile in path.glob("*.cif"):
        filename_header = ciffile.stem  # Extract filename without extension
        idlist.append(filename_header)

        # Structural relaxation (if any)
        chgnet_structure = Structure.from_file(ciffile)

        # Relaxation step if specified
        if structure_opt:
            result = relaxer.relax(chgnet_structure, verbose=False)
            print(f"\nCHGNet took {len(result['trajectory'])} steps. Relaxed structure:")
            final_structure = result['final_structure']
        else:
            final_structure = chgnet_structure

        # Convert cif to vasp for matrices
        filepath = path / f"{filename_header}.vasp"
        chgnet_structure.to(filename=filepath, fmt="poscar")

        # Read the VASP file
        virtual = readfile(str(filepath))
        if verbose: printvaspdata(virtual)

        # Run matrix bonding average
        runiq, bme, unaverageness = matrix_bonding_average(virtual, 'element', bond_threshold, bme_correlated=bme_targ, verbose=False)

        # Append to data lists
        datalist.append(virtual)

        # Append to bond lists
        if verbose: display(bme)
        bondlist.append(bme)

        # Bonding matrix distance
        bd_metric = np.linalg.norm(bme - bme_targ, 'fro')
        print("Bonding matrix distance: ", bd_metric, " Unaverageness: ", unaverageness)
        bdlist1.append(bd_metric)
        bdlist2.append(unaverageness)
        bdlist.append(bd_metric + pauling_weight * unaverageness)

        # Append to composition lists
        virtual_comp = composition(virtual, verbose=False)
        complist.append(virtual_comp)

        try: 
            compdist = distance(original_comp, virtual_comp)
        except: 
            print("Elements mismatch! (possibly all atoms of certain species removed?)")
            compdist = "err"  # Unlikely due to compounded probabilities!
        cdlist.append(compdist)

        # CHGNET total energy prediction
        prediction = chgnet.predict_structure(final_structure)
        totalenergy = float(prediction['e'])
        if verbose: 
            print(f"CHGNet-predicted energy (eV/atom):\n{totalenergy}\n")
        E_list.append(totalenergy)

    summary = pd.DataFrame({'Name': idlist, 'BondingProfile': bondlist, 'MatrixDistance' : bdlist1, 'Unaverageness' : bdlist2, 'BondDiff' : bdlist, 'Composition': complist, 'CompDiff': cdlist, 'FormationEnergy': E_list})
    summary.to_csv(Path(path) / 'VirtualLibrary.csv')
    if verbose: display(summary)
    print("Summary of results written to VirtualLibrary.csv")
    return datalist, summary
