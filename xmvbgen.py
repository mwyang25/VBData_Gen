import os
import re 
import sys

import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors

import pandas as pd 
import numpy as np

PTDic = {'1': 'H', '2': 'He', '3': 'Li', '4': 'Be', '5': 'B', '6': 'C', '7': 'N', '8': 'O', '9': 'F', '10': 'Ne',
         '11': 'Na', '12': 'Mg', '13': 'Al', '14': 'Si', '15': 'P', '16': 'S', '17': 'Cl', '18': 'Ar', '19': 'K', '20': 'Ca',
         '21': 'Sc', '22': 'Ti', '23': 'V', '24': 'Cr', '25': 'Mn', '26': 'Fe', '27': 'Co', '28': 'Ni', '29': 'Cu', '30': 'Zn',
         '31': 'Ga', '32': 'Ge', '33': 'As', '34': 'Se', '35': 'Br', '36': 'Kr', '37': 'Rb', '38': 'Sr', '39': 'Y', '40': 'Zr',
         '41': 'Nb', '42': 'Mo', '43': 'Tc', '44': 'Ru', '45': 'Rh', '46': 'Pd', '47': 'Ag', '48': 'Cd', '49': 'In', '50': 'Sn',
         '51': 'Sb', '52': 'Te', '53': 'I', '54': 'Xe', '55': 'Cs', '56': 'Ba',
         '74': 'W', '76': 'Os', '77': 'Ir', '78': 'Pt', '79': 'Au', '80': 'Hg', '81': 'Tl', '82': 'Pb', '83': 'Bi'}

def travel_atom(atom,another_atom_idx,fraglist):
    for neighbor in atom.GetNeighbors():
        if neighbor.GetIdx() != another_atom_idx and neighbor.GetIdx() not in fraglist:
            fraglist.add(neighbor.GetIdx())
            travel_atom(neighbor,another_atom_idx,fraglist)
    return fraglist

def get_elenums(mol,fraglist):

    ele_nums=0
    for idx in fraglist:
        atom_ele = mol.GetAtomWithIdx(idx).GetAtomicNum()
        ele_nums=ele_nums+atom_ele
    return ele_nums


def get_info(smile,nei_atom_id,sec_atom_id):

    mol = Chem.AddHs(Chem.MolFromSmiles(smile))
    atom1=mol.GetAtomWithIdx(nei_atom_id)
    atom2=mol.GetAtomWithIdx(sec_atom_id)

    frag1list,frag2list = set(),set()
    frag1list.add(nei_atom_id)
    frag2list.add(sec_atom_id)
    frag1_m1=list(travel_atom(atom1,sec_atom_id,frag1list))
    frag2_m1=list(travel_atom(atom2,nei_atom_id,frag2list))
    frag1 = [i+1 for i in frag1_m1]
    frag2 = [i+1 for i in frag2_m1]
    frag1len=len(frag1)
    frag2len=len(frag2)

    orb1 = int((get_elenums(mol,frag1_m1)-1)/2)
    orb2 = int((get_elenums(mol,frag2_m1)-1)/2)
    return frag1,frag2,frag1len,frag2len,orb1,orb2


smi=sys.argv[1]
bond_index=sys.argv[2]
fileName = smi + '.log'
if not os.path.isfile(fileName):
    print(f'No %s file exists'%fileName)
else:
    f = open(fileName, 'r')
    s = f.read()
    f.close()
    centNumList = []
    coordObj = re.findall(r'^\s+(\d+)\s+(\d+)\s+\d\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)$', s, flags=re.M)
    mol = Chem.AddHs(Chem.MolFromSmiles(smi))
    bond=mol.GetBondWithIdx(int(bond_index))
    begin_index=bond.GetBeginAtomIdx()
    end_index=bond.GetEndAtomIdx()

    frag_st, frag_nd, frag1len,frag2len, orb1,orb2 = get_info(smi, begin_index, end_index)
    
    file_out_path = smi + bond_index
    file_out = file_out_path +'cov.xmi'
    f = open(file_out,'w')
    print( smi + '\t'+ bond_index,file=f)
    print('$CTRL \nstr=cov nao=2 nae=2 \norbtyp=hao frgtyp=atom \nint=libcint basis=def2-svp \niscf=5 iprint=3 ITMAX=2000 \n$END',file=f)
    print('$FRAG',file=f)
    print('%d %d'%(frag1len,frag2len),file=f)
    print(' '.join([str(a) for a in frag_st]),file=f)
    print(' '.join([str(a) for a in frag_nd]),file=f)
    print('$END',file=f)
    print('$ORB',file=f)
    print('1*%d  1*%d  1  1'%(orb1,orb2),file=f)
    for i in range(orb1):
        print('1',file=f)
    for i in range(orb2):
        print('2',file=f)
    print('1',file=f)
    print('2',file=f)
    print('$END',file=f)
    print('$GEO',file=f)
    if coordObj is not None:
        for centNum in coordObj:
            centNumList.append(eval(centNum[0]))
        for coord in coordObj[-max(centNumList):]:
            print(PTDic[coord[1]] + coord[2].rjust(12) + coord[3].rjust(12) + coord[4].rjust(12),file=f)
    else: print('GEO ERROR')
    print('$END',file=f)
    f.close()
    file = file_out.replace('cov.xmi','.xmi')
    with open(file_out, 'r') as f1,open(file, 'w') as f2:
            for line in f1:
                f2.write(re.sub('str=cov','str=full',line))
    f1.close()
    f2.close()