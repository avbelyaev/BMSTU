from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import rdkit.Chem.rdMolAlign as ChemAlign
from rdkit.Chem import PyMol

# all the shit from above can be installed with anaconda unfortunately

EPINEPHRINE = 'C1CN=C(N1)NC2=C(C=CC=C2Cl)Cl'
CLONIPIDINE = 'CNCC(C1=CC(=C(C=C1)O)O)O'


def main():
    smiles = [EPINEPHRINE, CLONIPIDINE]

    mols = [Chem.MolFromSmiles(m) for m in smiles]
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

    # somehow this shit is the key
    pyO3A = AllChem.GetO3A(mols[0], mols[1], refMol1, refMol2)
    print('align')
    print(pyO3A.Align())
    print(pyO3A.Matches())

    Chem.MolToMolFile(mols[0], '0_.mol')
    Chem.MolToMolFile(mols[1], '1_.mol')

    ff = AllChem.UFFGetMoleculeForceField(mols[0])
    ff.Initialize()
    ff.Minimize(maxIts=200)

    rmsd = ChemAlign.AlignMol(mols[0], mols[1], atomMap=pyO3A.Matches())
    print('rmsd')
    print(rmsd)

    print('energy')
    print(ff.CalcEnergy())

    # launch pymol in server mode: `./pymol -R`
    # it resides in ~/anaconda/bin
    v = PyMol.MolViewer()
    v.ShowMol(mols[0])
    v.GetPNG(h=200)


if __name__ == '__main__':
    print('make sure pymol is launched in server mode')
    main()
