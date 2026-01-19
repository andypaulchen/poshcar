# distance: Find distance between sites in cell, in quoted Cartesian coordinates

# Andy Paul Chen, Monday, 9 August 2021, Little Italy, Cleveland, Ohio (National Day of Singapore)
# 2 Nov 2023: I can't believe I am back again
# 19 Jan 2026: Back here to speed up these functions like 60 timwa

from poshcar.cartesian import * # Cartesian coordinates package

# Cite: https://doi.org/10.1039/B801115J; https://doi.org/10.1002/chem.200800987 (Bk onwards; single bond)
covalent_radius_pm = [31+5, 28, 128+7, 96+3, 84+3, 76+1, 71+1, 66+2, 57+3, 58, 166+9, 141+7, 121+4, 111+2, 107+3, 105+3, 102+4, 106+10, 203+12, 176+10, 170+7, 160+8, 153+8, 139+5, 161+8, 152+6, 150+7, 124+4, 132+4, 122+4, 122+3, 120+4, 119+4, 120+4, 120+3, 116+4, 220+9, 195+10, 190+7, 175+7, 164+6, 154+5, 147+7, 146+7, 142+7, 139+6, 145+5, 144+9, 142+5, 139+4, 139+5, 138+4, 139+3, 140+9, 244+11, 215+11, 207+8, 204+9, 203+7, 201+6, 199, 198+8, 198+6, 196+6, 194+5, 192+7, 192+7, 189+6, 190+10, 187+8, 175+10, 187+8, 170+8, 162+7, 151+7, 144+4, 141+6, 136+5, 136+6, 132+5, 145+7, 146+5, 148+4, 140+4, 150, 150, 260, 221+2, 215, 206+6, 200, 196+7, 190+1, 187+1, 180+6, 169+3, 168, 168, 165, 167, 173, 176, 161, 157, 149, 143, 141, 134, 129, 128, 121, 122, 136, 143, 162, 175, 165, 157] # add one SD

nobond = [0,2,1,1,0,0,0,0,0,2,1,1,1,0,0,0,0,2,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2]
# 1: Metal (M-M bond forbidden)
# 2: Noble gas (any bond forbidden)

def is_nearest(data, c, verbose = True):
    # read c (array of three floats)
    # return index of closest site to coordinates
    # display distance
    # allcoord extraction
    ns, allcoord = images(data)
    
    # analysis
    for i in range(ns):
        for virtual in range(27):
            dist = distance(np.array(c), allcoord[virtual][i])
            if i==0 and virtual==0:
                mindex = 1
                mindist = dist
            else:
                if dist < mindist:
                    mindex = i+1
                    mindist = dist
                    
    if verbose: print("Minimum distance agreement: ", flpr.format(mindist), ", Index #", mindex)
    return mindex

# Andy Paul Chen, 2 November 2023
# My first coding after surviving war in Israel. They shot as us with rockets!
# Can people here ever catch a break??

def images(data):
    # Generate 27 images of an atomic site around an original unit cell (Cartesian coordinates)
    # Useful dor dealing with periodic boundary artefacts
    # Output to a unique 27-cell data structure - list of 27 lists of coordinates
    # Convert to Cartesian coordinates
    if not is_cart(data): data = switchcart(data, verbose = False)
    # Choose starting file line for coordinates according to site-disorder (selective dynamics) tag
    dcindex = 8 if is_seldyn(data) else 7
    B = basis(data) # Read lattice vectors
    
    # Read atomic coordinates
    monocoord = []
    for line in range(len(data)):
        if line > dcindex:
            dirs = np.array(re.findall(r"-?\d+\.\d+", data[line].strip()))[:3]
            dirs = dirs.astype(float) # convert from string to float
            monocoord.append(dirs)
    ns = len(monocoord) # Number of atoms in this unit cell

    # Let me write this down here real quicc. We are going to make a supercell
    # but not follow through all the processes of the supercell function. We will
    # extend the list of atoms to the 27 surrounding unit cells only. Closest
    # neighbours which are images on the list will be mapped to the original ion.
    
    # Make 26 shadow-coordinate cells
    allcoord = [deepcopy(monocoord)]
    for it in range(26): allcoord += [deepcopy(monocoord)]
    
    # Move 26 shadow-coordinate cells to their respective positions
    for i in [9,10,11,12,13,14,15,16,17]: 
        for j in range(ns): allcoord[i][j] = allcoord[i][j] + B[0] # +a
    for i in [18,19,20,21,22,23,24,25,26]: 
        for j in range(ns): allcoord[i][j] = allcoord[i][j] - B[0] # -a
    for i in [3,4,5,12,13,14,21,22,23]: 
        for j in range(ns): allcoord[i][j] = allcoord[i][j] + B[1] # +b
    for i in [6,7,8,15,16,17,24,25,26]: 
        for j in range(ns): allcoord[i][j] = allcoord[i][j] - B[1] # -b
    for i in [1,4,7,10,13,16,19,22,25]: 
        for j in range(ns): allcoord[i][j] = allcoord[i][j] + B[2] # +c
    for i in [2,5,8,11,14,17,20,23,26]: 
        for j in range(ns): allcoord[i][j] = allcoord[i][j] - B[2] # -c
            
    return ns, allcoord

def elemindices(data, verbose = True):
    if is_seldyn(data) and data[7].lower() == 'site-disordered structure\n': # site-disordered structures require site-specific symbols
        if verbose: print("Site-disordered structure detected!")
        dcindex = 8
        spp = []
        res = []
        occ = []
        for line in range(len(data)):
            if line > dcindex:
                spp.append(re.sub(r'[^a-zA-Z]', '', np.array(re.findall(r'\S+', data[line].strip()))[3]))
                res.append(np.array(re.findall(r'\S+', data[line].strip()))[3])
                occ.append(np.array(re.findall(r'\S+', data[line].strip()))[4])
    else:
        numatoms = np.array(re.findall('\d+', data[6].strip()))
        numatoms = numatoms.astype(int) # convert all elements to int
        elemlist = data[5].split()
        spp = list(itertools.chain.from_iterable(itertools.repeat(elemlist[i], numatoms[i]) for i in range(len(elemlist))))
        res = spp
        occ = np.ones(len(res))
    res_indexed = deepcopy(res)
    for i in range(len(res)): res_indexed[i] = res[i] + "-" + str(i+1)
    atomspp = pd.DataFrame({'Species': spp, 'Wyckoff Site': res, 'POSCAR Site': res_indexed, 'Occupancy': np.array(occ).astype(float)})
    return atomspp

def matrix_distances(data, verbose=True):
    # Distances between atom pairs (in terms of vectors and also distance)
    
    atomspp = elemindices(data, verbose=verbose)
    ns, allcoord = images(data)  # allcoord: (27, ns, 3)

    R0 = np.asarray(allcoord[0])        # (ns, 3)
    Rimg = np.asarray(allcoord)         # (27, ns, 3)

    # vectors: (27, ns, ns, 3)
    vectors_matrix = Rimg[:, None, :, :] - R0[None, :, None, :]
    # distances: (27, ns, ns)
    distance_matrix = np.linalg.norm(vectors_matrix, axis=-1)

    if verbose:
        df = pd.DataFrame(distance_matrix[0], columns=list(atomspp['POSCAR Site']), index=list(atomspp['POSCAR Site']))
        print("Distances between atoms (native cell only):")
        display(df)

    return vectors_matrix, distance_matrix

def matrix_bonding(data, tolerance = 0.1, verbose=True):
    # Use distances, covalent atomic radii and thresholding to identify "bonding" pairs of atoms
    
    distance_matrix = matrix_distances(data, verbose=False)[1]  # (27, ns, ns)
    atomspp = elemindices(data, verbose=verbose)
    species = np.asarray(atomspp['Species'])
    ns = species.size

    # Atomic numbers (vectorized)
    atomnos = np.array([periodic_table.index(s) for s in species], dtype=int)

    covrad = np.asarray(covalent_radius_pm) / 100.0
    radii = covrad[atomnos]                        # (ns,)
    nob = np.asarray(nobond)[atomnos]              # (ns,)

    # Pairwise radius sums
    rsum = radii[:, None] + radii[None, :]         # (ns, ns)
    cutoff = rsum * (1.0 + tolerance)

    # Distance condition
    bonded = distance_matrix <= cutoff[None, :, :]  # (27, ns, ns)

    # Remove self-bonds
    eye = np.eye(ns, dtype=bool)
    bonded &= ~eye[None, :, :]

    # nobond rules
    bonded &= ~((nob[:, None] == 1) & (nob[None, :] == 1))[None, :, :]
    bonded &= ~((nob[:, None] == 2) | (nob[None, :] == 2))[None, :, :]
    bonding_matrix = bonded.astype(np.int8)

    if verbose:
        df = pd.DataFrame(bonding_matrix[0], columns=list(atomspp['POSCAR Site']), index=list(atomspp['POSCAR Site']))
        print("Bonding between atoms (native cell only):")
        display(df) # Display example (native cell)

    return bonding_matrix

def matrix_bonding_average(data, mode, tolerance = 0.1, bme_correlated = 'amaiwana', verbose = True):
    # Average bonding matrix, can also be used for site-disordered stuff
    # mode: string input, first letter for [s]ite classification or [e]lement species classification
    # Edit of 10 Sep 2024: deviations from the mean are recorded in a separate matrix to output
    # bme_correlated specified if another averaged matrix needs to be used (use this in the case of correlated disorder)
    
    atomspp = elemindices(data, verbose=verbose)
    weights = np.asarray(atomspp['Occupancy'])

    # Sum bonding matrices once
    bms = matrix_bonding(data, tolerance, verbose=False).sum(axis=0)

    # Classification
    classification = mode[0].upper()
    if classification == 'E':
        labels = np.asarray(atomspp['Species'])
    elif classification == 'S':
        labels = np.asarray(atomspp['Wyckoff Site'])
    else:
        labels = np.asarray(atomspp['Species'])

    # Unique labels, preserve order
    uniq, first_idx = np.unique(labels, return_index=True)
    order = np.argsort(first_idx)
    uniq = uniq[order]

    us = len(uniq)

    # Build index groups
    groups = [np.where(labels == u)[0] for u in uniq]

    # Preallocate
    bme = np.zeros((us, us))

    # Core computation (vectorized blocks)
    for i, gi in enumerate(groups):
        wi = weights[gi]
        wi_sum = wi.sum()

        for j, gj in enumerate(groups):
            wj = weights[gj]
            Bij = bms[np.ix_(gi, gj)]

            bme[i, j] = np.einsum("i,ij,j->", wi, Bij, wj) / wi_sum

    # Unaverageness without explicit bma/bmc
    if isinstance(bme_correlated, str):
        target = bme
    else:
        target = np.asarray(bme_correlated)

    bme_expanded = np.zeros_like(bms)
    for i, gi in enumerate(groups):
        for j, gj in enumerate(groups):
            bme_expanded[np.ix_(gi, gj)] = target[i, j]

    unaverageness = np.linalg.norm(bme_expanded - bms, ord="fro")

    if verbose:
        df = pd.DataFrame(bme, index=uniq, columns=uniq)
        print("bme - average local coordination by environment:")
        display(df)
        print("Unaverageness:", unaverageness)

    return uniq, bme, unaverageness

def crashtest(data, tolerance = 0.1, verbose = True):
    # Determine if an atom pair is too close to one another
    # tolerance is a percentage value (0.0-1.0)
    # 26 Jul 2024: returns crash matrix instead of bool output
    # This function is not yet numpy-optimized (too useless for now)

    tolerance = min(abs(tolerance),1.0)
    if verbose: print ("Tolerance: ", tolerance)
    accept = True
    
    vector_matrix, distance_matrix = matrix_distances(data, verbose) # Extract distance matrix, atom indices
    atomspp = elemindices(data, verbose)
    ns = len(list(atomspp['POSCAR Site']))
    covalent_radius_a = [x/100 for x in covalent_radius_pm]
    # Retrieve atomic numbers
    atomnos = np.zeros(ns).astype(int)
    for k in range(ns): atomnos[k] = periodic_table.index(list(atomspp['Species'])[k]) # extract atomic numbers
    
    # Construct bonding matrix
    # interatomic distances smaller than (radius(1) + radius(2))*(1+tolerance factor) are bonded
    crashing_matrix = np.zeros((ns,ns)).astype(int)
    
    for i in range(ns): # i: from atom
        for j in range(ns): # j: to atom
            ei = atomnos[i] # atomic number of species i
            ej = atomnos[j] # atomic number of species j
            for virtual in range(27): # for each virtual image
                if distance_matrix[virtual][i][j] <= (covalent_radius_a[ei]+covalent_radius_a[ej])*(1-tolerance):
                    if (i != j):
                        crashing_matrix[i][j] = 1
                        accept = False
                    
    df = pd.DataFrame(crashing_matrix, columns=list(atomspp['POSCAR Site']), index=list(atomspp['POSCAR Site']))
    if verbose:
        if accept == False: 
            print("Crashing between atoms:")
            display(df) # Display example
        else: print("No crashing between atoms!")
    return crashing_matrix