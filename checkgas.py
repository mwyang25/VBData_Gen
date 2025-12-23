import os
import re

def check_opt(string,i,file_name): 
    pattern = r"(?<=Converged\?)(?:\n[^\n]+){4}"

    matches = re.findall(pattern, string)
    if matches:
        # cont = matches[i].strip()
        cont = matches[i]
    else:
        print("No match found.",file_name)
    return cont

def match_freq(string,filename):
    freqObj = re.findall(r'^\s+Frequencies\s--\s+(-?\d+\.\d+.+)$', string, flags=re.M)
    freqs = []
    if freqObj is not None:
        for freq in freqObj:
            split_freqs = freq.split()
            numbers = [float(num) for num in split_freqs]
            freqs.extend(numbers)
        for frequency in freqs:
            if frequency < 0:
                print("Vibrational frequencies:",filename)
                # print(frequency)
            
    else:
        print("Match Frequencies Error") 

 
for item in os.listdir(os.getcwd()):
    if os.path.isdir(item):
        filename = item + '/gaussian/' +item+'.log'
        f = open(filename, 'r')
        s = f.read()
        f.close()
        try:
            content = check_opt(s,-1,item)
            content1 = check_opt(s,-2,item)
            # print(f'Content:\n {content}')
            match_freq(s,item)
            if content.count('YES') != 4 and content1.count('YES') != 4:
                print(item)

        except:
            print('erro:',item)