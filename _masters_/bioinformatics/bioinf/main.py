from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import  GraphDescriptors
import numpy as np


data = open('data.txt', 'r').readlines()

m = [Chem.MolFromSmiles(data[i])for i in range(len(data))]
for i in range(len(m)):
    m[i] = Chem.AddHs(m[i])
descriptors = Descriptors._descList

desc_names = []
for d in descriptors:
    desc_names.append(d[0])

output = np.zeros(shape=(len(m), len(desc_names)))
bal_descriptors = np.zeros(shape=(len(m), 1))

for i in range(len(m)):
    for j in range(len(descriptors)):
        output[i][j] = (descriptors[j][1](m[i]))
    bal_descriptors[i] = GraphDescriptors.BalabanJ(m[i])

np.savetxt('bal_desc.txt', bal_descriptors, delimiter=',')
np.savetxt('data_desc.txt', output, delimiter=',')