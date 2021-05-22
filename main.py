import rdkit
from rdkit import Chem

tokens = ('@<TRIPOS>ATOM', '@<TRIPOS>BOND')





def main():
    ''''''
   # a = Chem.SDMolSupplier('mol.mol2')
    a = Chem.MolFromMol2File('mol.mol2')



# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()