import os
import glob
import sys
import argparse
from pymol import cmd
import pandas as pd

class dRMSDAnalyzer:
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
    
    def calculate_rmsd(self):
        """Calculates the RMSD for each product file compared to the TS1 file."""
        if not self.ts_file or not self.prod_files:
            print('unable to find files')
            return
        
        for prod_file in self.prod_files:
            cmd.load(self.ts_file, "ts_mol")
            prod_name = os.path.basename(prod_file)
            cmd.load(prod_file, "prod_mol")
            rmsd = cmd.align("prod_mol", "ts_mol")[0]
            self.rmsd_values[prod_name] = rmsd
            cmd.delete("all")
    
    def compute_product_ratios(self):
        """Computes the product ratio based on RMSD values."""
        total_rmsd = sum(self.rmsd_values.values())
        if total_rmsd == 0:
            self.ratios = {key: 0 for key in self.rmsd_values}
        else:
            self.ratios = {key: 1-(rmsd / total_rmsd) for key, rmsd in self.rmsd_values.items()}
    
    def run_analysis(self):
        """Runs the complete RMSD analysis and prints results."""
        self.find_files()
        if not self.ts_file:
            print("No TS1 file found.")
            return
        if not self.prod_files:
            print("No product files found.")
            return
        
        self.calculate_rmsd()
        self.compute_product_ratios()

        print('TS1: '+self.ts_file)
        
        print("RMSD Values:")
        for name, rmsd in self.rmsd_values.items():
            print(f"{name}: {rmsd:.2f} Å")
        
        print("\nProduct Ratios:")
        for name, ratio in self.ratios.items():
           print(f"{name}: {ratio:.2f}")


        data = {"Product": list(self.rmsd_values.keys()), "RMSD (Å)": list(self.rmsd_values.values()), 
                "Ratio": list(self.ratios.values())}
        self.result = pd.DataFrame(data)
        

def print_usage():
    usage_msg = """Usage: python dRMSD.py <folder_path>
Usage: <folder_path> should contain _TS.mol and _prod.mol files
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
        analyzer = dRMSDAnalyzer(folder_path)
        analyzer.run_analysis()
    except Exception as e:
        print("An error occurred while running the script.\n")
        print_usage()
        print(f"\nError details: {e}")
        sys.exit(1)



