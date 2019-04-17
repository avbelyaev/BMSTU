from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors


METH = "CC(CC1=CC=CC=C1)NC"
BENZENE = "C1=CC=CC=C1"
NAFTALENE = "C1=CC=C2C=CC=CC2=C1"
ANTRCENE = "C1=CC=C2C=C3C=CC=CC3=CC2=C1"
TETRACENE = "C1=CC=C2C=C3C=C4C=CC=CC4=CC3=CC2=C1"
PENTACENE = "C1=CC=C2C=C3C=C4C=C5C=CC=CC5=CC4=CC3=CC2=C1"
HEXACENE = "C1=CC=C2C=C3C=C4C=C5C=C6C=CC=CC6=CC5=CC4=CC3=CC2=C1"
OCTACENE = "C1=CC=C2C=C3C=C4C=C5C=C6C=C7C=C8C=CC=CC8=CC7=CC6=CC5=CC4=CC3=CC2=C1"

THINGS = [METH, BENZENE, NAFTALENE, ANTRCENE, TETRACENE, PENTACENE, HEXACENE, OCTACENE]


for thing in THINGS:
    m = Chem.MolFromSmiles(BENZENE)
    # add hydrogen
    m_hydrogenized = Chem.AddHs(m)

    Draw.MolToFile(m_hydrogenized, 'mol.png')

    for desc in Descriptors._descList:
        desc_name = desc[0]
        desc_value = desc[1](m_hydrogenized)
        print(str(desc_value) + "\t" + desc_name)

# http://www.qsar4u.com/files/rdkit_tutorial/rdkit.html
