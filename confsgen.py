import os
import sys
import re

import numpy as np
import pandas as pd
import random

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw


def gen_conformers(mol, num_conformers,numThreads):
    """
    生成分子的多个构象
    """
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=num_conformers, numThreads=numThreads, randomSeed=42)
    for cid in cids:
        AllChem.MMFFOptimizeMolecule(mol, confId=cid)
    print(len(cids))
    return mol

def get_rmsd(mol):
    rmslist = []
    AllChem.AlignMolConformers(mol, RMSlist=rmslist)
    
    return rmslist

if __name__=='__main__':
    smi = sys.argv[1]
    m3d = Chem.AddHs(Chem.MolFromSmiles(smi))
    cids = AllChem.EmbedMultipleConfs(m3d, numConfs=500, numThreads=0, randomSeed=42)
    res = AllChem.MMFFOptimizeMoleculeConfs(m3d,numThreads=0)
    energies = [conf[1] for conf in res]
    min_energy_index = energies.index(min(energies))
    min_energy_conf = m3d.GetConformer(min_energy_index)
    min_energy_conf_position = min_energy_conf.GetPositions()
    rmsd = get_rmsd(m3d)
    # print(min_energy_conf_position)
    print(Chem.MolToMolBlock(m3d, confId=min_energy_index))
    writer = Chem.SDWriter(smi+'.sdf')
    writer.write(m3d, confId=min_energy_index)
    writer.flush()
    file = open(smi + 'conf.txt','w')
    file.write(str(res))
    file.write('\n\n')
    file.write(str(energies))
    file.write('\n\n')
    file.write(str(res[min_energy_index][0]))
    file.write('\n\n')
    file.write(str(rmsd))
    file.close() 

    #xyz文件生成
    nums=m3d.GetNumAtoms()
    file_in = open(smi+'.sdf', mode="r")
    content = file_in.read()
    file_in.close()
    match = re.search(r"^ {3,}-?\d.*(?:\r?\n {3,}-?\d.*)*", content, re.M)
    datasplit = []
    xyzs=[]

    if match:
        for line in match.group().splitlines():
            datasplit.append([part for part in line.split()][:4])
    datasplit=np.array(datasplit)
    temp = np.roll(datasplit,1,axis=1)
    xyzs=pd.DataFrame(temp)

    filename =smi + '.xyz'
    file_out = open(filename,'w')
    print(nums,file=file_out)
    file_out.write(smi)
    file_out.write('\n')
    file_out.close()

    xyzs.to_csv(filename,sep=' ',index=0,header=0,mode='a')
