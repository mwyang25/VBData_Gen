import sys
import pandas as pd

#从分子的XYZ文件中提取坐标信息，生成半经验优化的MNDO的输入文件
PTDic = {'1': 'H', '2': 'He', '3': 'Li', '4': 'Be', '5': 'B', '6': 'C', '7': 'N', '8': 'O', '9': 'F', '10': 'Ne',
         '11': 'Na', '12': 'Mg', '13': 'Al', '14': 'Si', '15': 'P', '16': 'S', '17': 'Cl', '18': 'Ar', '19': 'K', '20': 'Ca',
         '21': 'Sc', '22': 'Ti', '23': 'V', '24': 'Cr', '25': 'Mn', '26': 'Fe', '27': 'Co', '28': 'Ni', '29': 'Cu', '30': 'Zn',
         '31': 'Ga', '32': 'Ge', '33': 'As', '34': 'Se', '35': 'Br', '36': 'Kr', '37': 'Rb', '38': 'Sr', '39': 'Y', '40': 'Zr',
         '41': 'Nb', '42': 'Mo', '43': 'Tc', '44': 'Ru', '45': 'Rh', '46': 'Pd', '47': 'Ag', '48': 'Cd', '49': 'In', '50': 'Sn',
         '51': 'Sb', '52': 'Te', '53': 'I', '54': 'Xe', '55': 'Cs', '56': 'Ba',
         '74': 'W', '76': 'Os', '77': 'Ir', '78': 'Pt', '79': 'Au', '80': 'Hg', '81': 'Tl', '82': 'Pb', '83': 'Bi'}


if __name__=='__main__':
    smi=sys.argv[1]
    data = pd.read_table(smi + '.xyz',skiprows=2,sep=' ', header=None)
    symbols = data.iloc[:,0]
   
    atom_order=[]
    keys = list (PTDic.keys()) 
    values= list (PTDic.values())
    for symbol in symbols:
        atom_order.append(keys[values.index(symbol)])

    data[0] = atom_order
    ones=[1 for index in range(len(data))]
    data['x'], data['y'], data['z']= ones, ones, ones
    data=data.reindex(columns=[0,1,'x',2,'y',3,'z'])
    data.iloc[0,2], data.iloc[0,4], data.iloc[0,6] = 0, 0, 0 
    data.loc[len(data.index)] = [0, 0.000000, 0, 0.000000, 0, 0.000000, 0]
    data[[0,'x','y','z']] = data[[0,'x','y','z']].astype('int')

    filename = smi + '.inp'
    file_out = open(filename,'w')
    control = 'iform=1  igeom=1 iop=-23 jop=3 nsav7=4 nsav13=2'
    print(control,file=file_out)
    file_out.write('\n')
    file_out.write('OM2')
    file_out.write('\n')
    file_out.close()

    data.to_csv(filename,sep=' ',index=0,header=0,mode='a')
