import re
import os
import json
import pandas as pd
import numpy as np

def get_weights(string):
    pattern = r'WEIGHTS OF STRUCTURES\s+\**\s*\n\s*((?:\s*\d+\s+([-]?\d+\.\d+)\s+\**\s+(?:\S+\s*)+\n*)+)'
    matches = re.findall(pattern, string)

    weights = ['%.2f' % float(match[1]) for match in re.findall(r'(\d+)\s+([-]?\d+\.\d+)', matches[0][0])[:3]]
    return weights

def get_RE(string1,string2):
    corrCom = re.compile('TOTAL VB ENERGY :' + '\s+(-?\d+\.\d+)')
    match = corrCom.search(string1)
    if match:
        energy_tt = corrCom.search(string1).group(1)
        energy_cov = corrCom.search(string2).group(1)
        RE = float(energy_tt) - float(energy_cov)
    
    return RE

def getfolder_list(path):
    dirs1=[]
    for item in os.scandir(path):
        if item.is_dir():
            dirs1.append(item.path)
        else: pass
    return dirs1


obj=json.load(open('./smi_file.json', 'r'))
keys = list (obj.keys()) 
values= list (obj.values())
df = pd.read_csv("./H.csv")
df['RE'] = 0

dirs = getfolder_list('./')

for dir in dirs:
    order = dir.split('/')[-1]
    smi = keys[values.index(order)]
    print(f"Processing: '{dir}' '{smi}'")
    infos = df[df['molecule'] == smi]
    path = './'+ order + '/xmvb/'
    dirs1 = getfolder_list(path)
    for dir1 in dirs1:
        index = dir1.rsplit('/', 1)[-1]
        print(index)
        path1 =dir1 + '/'
        fileName1 = path1 + order+index+'.xmo'
        fileName2 = path1 + order+index+'cov.xmo'
        if not os.path.exists(fileName1):
            print('Error not exist 1')
        elif not os.path.exists(fileName2):
            print('Error not exist 2')
        else:
            s1 = open(fileName1, 'r').read()
            s2 = open(fileName2, 'r').read()

            corrCom = re.compile('TOTAL VB ENERGY :' + '\s+(-?\d+\.\d+)')
            match1 = corrCom.search(s1)
            match2 = corrCom.search(s2)
            if match1 and match2:
                energy_tt = corrCom.search(s1).group(1)
                energy_cov = corrCom.search(s2).group(1)
                RE = float(energy_cov) - float(energy_tt)
                RE_kcal = round(RE*627.5094, 2)
                # print(RE_kcal)
                df['RE'] = np.where((df['molecule'] == smi) & (df['end_idx'] == int(index)), RE_kcal, df['RE'])
                df.loc[(df['molecule'] == smi) & (df['end_idx'] == int(index)), ['Weight_1', 'Weight_2', 'Weight_3']] = get_weights(s1)
            elif match2: 
                print('Error not matching 1')
            else:
                print('Error not matching 2')
        
df.to_csv('results.csv', index=False)
