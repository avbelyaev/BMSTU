from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import rdkit.Chem.rdMolAlign as ChemAlign
from rdkit.Chem import PyMol


def main():
    # mol1 = Chem.MolFromMolFile('epinephrine.mol')
    # mol2 = Chem.MolFromMolFile('clonidine.mol')

    # mols = [mol1, mol2]
    file = open('superposition.txt')
    smiles = file.read().split('\n')
    mols = [Chem.MolFromSmiles(m) for m in smiles][:2]
    mols = [Chem.AddHs(m) for m in mols]

    i = 0
    for m in mols:
        print(m)
        AllChem.EmbedMolecule(m)
        AllChem.MMFFOptimizeMolecule(m, maxIters=200)
        Chem.MolToMolFile(m, '{0}.mol'.format(i))

        Draw.MolToFile(m, 'm{0}.png'.format(i))
        i += 1

    refMol1 = AllChem.MMFFGetMoleculeProperties(mols[0])
    refMol2 = AllChem.MMFFGetMoleculeProperties(mols[1])

    pyO3A = AllChem.GetO3A(mols[0], mols[1], refMol1, refMol2)
    print('align')
    print(pyO3A.Align())
    print(pyO3A.Matches())

    Chem.MolToMolFile(mols[0], '0_.mol')
    Chem.MolToMolFile(mols[1], '1_.mol')
    # for m in mols:
    # ff = AllChem.MMFFGetMoleculeForceField(mols[0])
    # ff.Initialize()
    # ff.Minimize()
    ff = AllChem.UFFGetMoleculeForceField(mols[0])
    ff.Initialize()
    ff.Minimize(maxIts=200)

    rmsd = ChemAlign.AlignMol(mols[0], mols[1], atomMap=pyO3A.Matches())
    print('rmsd')
    print(rmsd)

    print('energy')
    print(ff.CalcEnergy())


    #  launch pymol in server mode. it lays in anaconda/bin. 'pymol -R'
    v = PyMol.MolViewer()
    v.ShowMol(mols[0])
    # v.ShowMol(mols[1])
    v.GetPNG(h=200)


if __name__ == '__main__':
    main()