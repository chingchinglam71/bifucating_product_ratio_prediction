from rdkit import Chem
from collections import Counter
import argparse
import math
import os
import sys
import glob

'''
A script for carrying out the bond length difference analysis as described in 

Relationships between Product Ratios in Ambimodal Pericyclic Reactions and Bond Lengths in Transition Structures
Zhongyue Yang, Xiaofei Dong, Yanmin Yu, Peiyuan Yu, Yingzi Li, Cooper Jamieson, and K. N. Houk
Journal of the American Chemical Society 2018 140 (8), 3061-3067
DOI: 10.1021/jacs.7b13562

ExtractBonds(), IdentifyChangedBonds() and ReadGeometries() function of this scripts are taken and
modified from https://github.com/Goodman-lab/VRAI-selectivity/blob/master/VRAI_selectivity_v7.py


'''

def ExtractBonds(rdkitmol):
    '''
    Extract all bonds in the rdkit molecule and the atom indexes of atoms in the bond

    '''
    bondlist = []

    for bond in rdkitmol.GetBonds():
        bondlist.append([])
        bondlist[-1] = [bond.GetBeginAtomIdx()+1, bond.GetEndAtomIdx()+1]

    return bondlist

def IdentifyChangedBonds(bondlist1, bondlist2):
    '''
    Find the uncommon bonds in the two product bond lists

    '''
    changed_bonds1 = []
    changed_bonds2 = []

    for bond_1 in bondlist1:
        bond_1.sort()

    for bond_2 in bondlist2:
        bond_2.sort()

    for bond_1 in bondlist1:
        if bond_1 in bondlist2:
            pass
        else:
            changed_bonds1 += [bond_1]

    for bond_2 in bondlist2:
        if bond_2 in bondlist1:
            pass
        else:
            changed_bonds2 += [bond_2]

    return changed_bonds1, changed_bonds2

def ReadGeometries(GMolFile):
    '''
    Finds optimisation steps in the Gaussian output file and extracts the coordinates for each step.
    If the output file is a frequency calculation only one geometry exists.
    Returns the coordinates as a list
    '''
    gausfile = open(GMolFile, 'r')
    GOutp = gausfile.readlines()
    gausfile.close()

    index = 0
    atoms = []
    coords = []

    for index in range(len(GOutp)):
        if index > 3:
            if len(GOutp[index].split()) < 8:
                break
            else:
                data = GOutp[index].split()
                atoms.append(data[3])
                coords.append([float(x) for x in data[0:3]])

    return atoms, coords

def calculate_bond_length(coords, bond_pairs):
    def distance(p1, p2):
        return math.sqrt(sum((c1 - c2) ** 2 for c1, c2 in zip(p1, p2)))
    
    return {tuple(pair): distance(coords[pair[0]-1], coords[pair[1]-1]) for pair in bond_pairs}

def calculate_B_A_ratio(bond3, bond2):
    exponent = -9.4 * (bond3 - bond2)
    B_A_ratio = math.exp(exponent)
    B = B_A_ratio / (1 + B_A_ratio)
    A = 1 - B
    return B,A

def process_molecule_files(filename1, filename3, filename4):
    
    rdkitmolTS1 = Chem.MolFromMolFile(filename1, removeHs=False, strictParsing=False, sanitize=False)
    rdkitmolP1 = Chem.MolFromMolFile(filename3, removeHs=False, strictParsing=False, sanitize=False)
    rdkitmolP2 = Chem.MolFromMolFile(filename4, removeHs=False, strictParsing=False, sanitize=False)

    atoms1, geometries1 = ReadGeometries(filename1)
    atoms3, geometries3 = ReadGeometries(filename3)
    atoms4, geometries4 = ReadGeometries(filename4)

    P1bonds = ExtractBonds(rdkitmolP1)
    P2bonds = ExtractBonds(rdkitmolP2)
    TS1bonds = ExtractBonds(rdkitmolTS1)

    unidenticalP1P2, unidenticalP2P1 = IdentifyChangedBonds(P1bonds, P2bonds)
    unidenticalTS1P1, unidenticalP1TS1 = IdentifyChangedBonds(TS1bonds, P1bonds)
    unidenticalTS1P2, unidenticalP2TS1 = IdentifyChangedBonds(TS1bonds, P2bonds)

    all_unbonds = [unidenticalTS1P1, unidenticalP1TS1, unidenticalTS1P2, unidenticalP2TS1]
    all_unbonds2 = [bond for sublist in all_unbonds for bond in sublist]
    
    sublist_counts = Counter(tuple(sublist) for sublist in all_unbonds2)
    all_unbonds3 = [list(sublist) for sublist, count in sublist_counts.items() if count == 1]
    bond_lengths = calculate_bond_length(geometries1, all_unbonds3)

    P1_bond, P2_bond = None, None
    for b in all_unbonds3:
        p1list = tuple(unidenticalTS1P1) + tuple(unidenticalP1TS1)
        p2list = tuple(unidenticalTS1P2) + tuple(unidenticalP2TS1)
        if b in p1list:
            P1_bond = b
        if b in p2list:
            P2_bond = b

    # Calculate the product ratio and round it to 2 decimal places
    P1r,P2r = calculate_B_A_ratio(bond_lengths.get(tuple(P1_bond)), bond_lengths.get(tuple(P2_bond)))
    
    # Print the product ratio (rounded to 2 decimal places)
    print(f"Product ratio (rounded to 2 decimal places): {P1r:.2f} / {P2r:.2f}")
    
    #print("Changed bonds (atom index starts from 1):", all_unbonds3)
    
    return calculate_B_A_ratio(bond_lengths.get(tuple(P1_bond)), bond_lengths.get(tuple(P2_bond)))



class dBondAnalyzer:
    def __init__(self, folder):
        self.folder = folder
        self.ts_file = None
        self.prod_files = []
        self.rmsd_values = {}
        self.ratios = {}
    
    def find_files(self):
        """Finds the TS1 and product files in the folder."""
        for file in glob.glob(os.path.join(self.folder, "*.mol")):
            if "_TS1.mol" in file:
                self.ts_file = file
            elif "_prod.mol" in file:
                self.prod_files.append(file)
    
    def calculate_ratio(self):
        self.find_files()
        
        print('TS1: ', self.ts_file.split('/')[-1])
        print('P1: ', self.prod_files[0].split('/')[-1])
        print('P2: ', self.prod_files[1].split('/')[-1])

        filename1 = self.ts_file
        filename3 = self.prod_files[0]
        filename4 = self.prod_files[1]
        
        self.P1r,self.P2r=process_molecule_files(filename1, filename3, filename4)


def print_usage():
    usage_msg = """Usage: python dBond.py <folder_path>
Usage: <folder_path> should contain _TS.mol and _prod.mol files
Usage: Please rename your structural file with the following suffix: '_TS1.mol', '_prod.mol'
Usage: this script only accepts .mol files"""
    print(usage_msg)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("folder_path", nargs="?", help="Path to the folder containing .mol files")
    parser.add_argument("--help", action="store_true", help="Show usage information")
    
    args = parser.parse_args()

    if args.help:
        print_usage()
        sys.exit(0)
    
    if not args.folder_path:
        print("Error: Missing folder path.\n")
        print_usage()
        sys.exit(1)
    
    folder_path = args.folder_path

    try:
        analyzer = dBondAnalyzer(folder_path)
        analyzer.calculate_ratio()
    except Exception as e:
        print("An error occurred while running the script.\n")
        print_usage()
        print(f"\nError details: {e}")
        sys.exit(1)
    






