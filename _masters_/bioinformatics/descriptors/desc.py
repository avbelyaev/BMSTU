# -*- coding: utf-8 -*-

from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors

# в кач-ве значения дескриптора выбрано количество кратных связей
# см картинки http://www.xumuk.ru/encyklopedia/2/5111.html
from scipy.stats import linregress


class ChemDescr:
    def __init__(self, descriptor_value, lambda_value, smiles_formula):
        self.descr = descriptor_value
        self.lambda_len = lambda_value
        self.smiles = smiles_formula


BUTADIENE = ChemDescr(2, 217, "C=CC=C")
OCTATERTAENE = ChemDescr(4, 312, "C=CC=CC=CC=C")
NAFTALENE = ChemDescr(2, 295, "C1=CC=C2C=CC=CC2=C1")
TETRACENE = ChemDescr(4, 460, "C1=CC=C2C=C3C=C4C=CC=CC4=CC3=CC2=C1")
PERILENE = ChemDescr(10, 432, "C1=CC2=C3C(=C1)C4=CC=CC5=C4C(=CC=C5)C3=CC=C2")
PENTACENE = ChemDescr(5, 540, "C1=CC=C2C=C3C=C4C=C5C=CC=CC5=CC4=CC3=CC2=C1")
CORONENE = ChemDescr(12, 421, "C1=CC2=C3C4=C1C=CC5=C4C6=C(C=C5)C=CC7=C6C3=C(C=C2)C=C7")
TERILENE = ChemDescr(15, 518, "C1=CC8=C7C(=C1)C2=CCC5C3=C2C(=CC=C3C4=CC=CC6=C4C5=CC=C6)C7=CC=C8")

TRAIN = [BUTADIENE, OCTATERTAENE, NAFTALENE, TETRACENE,
         PERILENE, PENTACENE, CORONENE, TERILENE]

BENZENE = ChemDescr(3, 255, "C1=CC=CC=C1")
ANTRACENE = ChemDescr(5, 370, "C1=CC=C2C=C3C=CC=CC3=CC2=C1")
HEXACENE = ChemDescr(6, 506, "C1=CC=C2C=C3C=C4C=C5C=C6C=CC=CC6=CC5=CC4=CC3=CC2=C1")
LICOPENE = ChemDescr(11, 693, "CC(=CCCC(=CC=CC(=CC=CC(=CC=CC=C(C)C=CC=C(C)C=CC=C(C)CCC=C(C)C)C)C)C)C")

TEST = [BENZENE, ANTRACENE, HEXACENE, LICOPENE]


# значение линейной регрессии вычислено автоматически в bio1.numbers
def linear_lambda_func(x):
    return 39 * x + 219


def main():
    expectation = []
    reality = []
    for t in TEST:
        expectation.append(t.lambda_len)

        predicted_lambda_len = linear_lambda_func(t.descr)
        reality.append(predicted_lambda_len)
        print(predicted_lambda_len)

    slope, intercept, r_value, p_value, std_err = linregress(expectation, reality)
    print(slope, intercept, r_value, p_value, std_err)
    # for t in TRAIN:
    #     m = Chem.MolFromSmiles(t.smiles)
    # add hydrogen
    # m_hydrogenized = Chem.AddHs(m)

    # Draw.MolToFile(m_hydrogenized, 'mol.png')

    # for desc in Descriptors._descList:
    #     desc_name = desc[0]
    #     desc_value = desc[1](m_hydrogenized)
    #     print(str(desc_value) + "\t" + desc_name)


if __name__ == '__main__':
    main()
