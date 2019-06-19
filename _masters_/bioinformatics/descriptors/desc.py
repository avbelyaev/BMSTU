# -*- coding: utf-8 -*-

from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors


# в кач-ве значения дескриптора выбрано количество бенольных колец
# см картинки http://www.xumuk.ru/encyklopedia/2/5111.html


class ChemDescr:
    def __init__(self, descriptor_value, lambda_value, smiles_formula):
        self.smiles = smiles_formula
        self.descr = descriptor_value
        self.light_len = lambda_value


BENZENE = ChemDescr(1, 255, "C1=CC=CC=C1")
NAFTALENE = ChemDescr(2, 275, "C1=CC=C2C=CC=CC2=C1")
ANTRCENE = ChemDescr(3, 370, "C1=CC=C2C=C3C=CC=CC3=CC2=C1")
TETRACENE = ChemDescr(4, 460, "C1=CC=C2C=C3C=C4C=CC=CC4=CC3=CC2=C1")
PENTACENE = ChemDescr(5, 580, "C1=CC=C2C=C3C=C4C=C5C=CC=CC5=CC4=CC3=CC2=C1")
HEXACENE = ChemDescr(6, 693, "C1=CC=C2C=C3C=C4C=C5C=C6C=CC=CC6=CC5=CC4=CC3=CC2=C1")
OCTACENE = ChemDescr(8, 782, "C1=CC=C2C=C3C=C4C=C5C=C6C=C7C=C8C=CC=CC8=CC7=CC6=CC5=CC4=CC3=CC2=C1")

THINGS = [BENZENE, NAFTALENE, ANTRCENE, TETRACENE, PENTACENE, HEXACENE, OCTACENE]

for t in THINGS:
    m = Chem.MolFromSmiles(t.smiles)
    # add hydrogen
    m_hydrogenized = Chem.AddHs(m)

    # Draw.MolToFile(m_hydrogenized, 'mol.png')

    for desc in Descriptors._descList:
        desc_name = desc[0]
        desc_value = desc[1](m_hydrogenized)
        print(str(desc_value) + "\t" + desc_name)

# http://www.qsar4u.com/files/rdkit_tutorial/rdkit.html
