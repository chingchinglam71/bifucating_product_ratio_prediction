from rdkit import Chem
from collections import Counter
import math
import os
import glob
import re
import argparse


'''

This script analyse TS/P1/P2 structure and extract KED value 
from RMCF_*.LOG from rmcf.py (Chem. Sci., 2021,12, 12682-12694)
to predict product ratio from bifucating reactions


Some parts of this scripts are taken from 
https://github.com/Goodman-lab/VRAI-selectivity/blob/master/VRAI_selectivity_v7.py


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

    all_unbonds = [unidenticalTS1P1 + unidenticalP1TS1, unidenticalTS1P2+ unidenticalP2TS1]
    all_unbonds2 = [bond for sublist in all_unbonds for bond in sublist]
    sublist_counts = Counter(tuple(sublist) for sublist in all_unbonds2)
    
    
    all_unbonds3 = [list(sublist) for sublist, count in sublist_counts.items() if count == 1]

    all_unbonds4 = [[item for item in sublist if item in all_unbonds3] for sublist in all_unbonds]

    P1_key_bonds = all_unbonds4[0]
    P2_key_bonds = all_unbonds4[1]
    

    return P1_key_bonds, P2_key_bonds, all_unbonds3


class RMCFAnalyzer:
    def __init__(self, folder):
        self.folder = folder
        self.ts_file = None
        self.prod_files = []
        self.rmsd_values = {}
        self.ratios = {}
        self.RMCF_log= None
    
    def find_files(self):
        """Finds the TS1 and product files in the folder."""
        for file in glob.glob(os.path.join(self.folder, "*.mol")):
            if "_TS1.mol" in file:
                self.ts_file = file
            elif "_prod.mol" in file:
                self.prod_files.append(file)

        for file in glob.glob(os.path.join(self.folder, "*.LOG")):
            if "RMCF_" in file:
                self.RMCF_log  = file
            
    
    def calculate_ratio(self, print_result=True):
        self.find_files()
        
        filename1 = self.ts_file
        filename3 = self.prod_files[0]
        filename4 = self.prod_files[1]
        
        P1_key_bonds, P2_key_bonds, all_unbonds3=process_molecule_files(filename1, filename3, filename4)

        # Define the target atom pairs
        target_pairs = {tuple(pair): None for pair in all_unbonds3 }
        
        # Read the log file
        
        log_file = self.RMCF_log

        with open(log_file, "r") as file:
            for line in file:
                match = re.match(r"\s*(\d+)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)", line)
                if match:
                    atom1, atom2, ked1, ked2, sum_ked = map(float, match.groups())
                    atom_pair = (int(atom1), int(atom2))
                    if atom_pair in target_pairs:
                        target_pairs[atom_pair] = sum_ked
        
        
        P1_ked = [target_pairs[tuple(pair)] for pair in P1_key_bonds]
        P2_ked = [target_pairs[tuple(pair)] for pair in P2_key_bonds]
        
        sum_P1_ked = sum(P1_ked)
        sum_P2_ked = sum(P2_ked)
        
        sum_ked=sum_P1_ked+sum_P2_ked
        self.P1_ratio=round(sum_P1_ked/sum_ked,2)
        self.P2_ratio=round(sum_P2_ked/sum_ked,2)

        if print_result==True:

            print('TS1: ', self.ts_file.split('/')[-1])
            print('P1: ', self.prod_files[0].split('/')[-1])
            print('P2: ', self.prod_files[1].split('/')[-1])
            print('RMCF LOG file: ', self.RMCF_log.split('/')[-1])
            
            print()
            
            print('P1 key bonds: ',str(P1_key_bonds))
            print('P2 key bonds: ',str(P2_key_bonds))
            print()
            # Print the results
            for pair, sum_ked in target_pairs.items():
                print(f"SUM(KED) for {pair}: {sum_ked:.3f}")
            
            print()
            print('P1 to P2 ratio: '+str(self.P1_ratio)+'/'+str(self.P2_ratio))




def print_usage():
    usage_msg = """Usage: python RMCFAnalyzer.py <folder_path>
Usage: <folder_path> should contain RMCF_*.LOG, _TS.mol and _prod.mol files
Usage: Please rename your structural file with the following suffix: '_TS.mol', '_prod.mol'
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
    	test=RMCFAnalyzer(folder_path)
    	test.calculate_ratio()
    except Exception as e:
        print("An error occurred while running the script.\n")
        print_usage()
        print(f"\nError details: {e}")
        sys.exit(1)


