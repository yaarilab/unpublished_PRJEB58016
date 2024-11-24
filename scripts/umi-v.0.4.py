#!/usr/bin/env python3
import pandas as pd
import time as T
import numpy as np
import scipy.cluster.hierarchy as hc
import argparse
import sys
import os
from copy import deepcopy
import pkg_resources
required = {'Levenshtein'}
installed = {pkg.key for pkg in pkg_resources.working_set}
missing = required - installed
if missing:
    os.system("python3 -m pip install Levenshtein")
    sys.path.insert(1, '/usr/local/lib64/python3.9/site-packages')
import Levenshtein as lv


parser = argparse.ArgumentParser()
parser.add_argument('-i', type=str) # input umi 
parser.add_argument('-o', type=str) # output clustered umi table 
args = parser.parse_args()

time_s = T.time()
df1 = pd.read_table(args.i, header=None, names=['umi_origin','seq'],sep=' ')  # reading the file containing umi's only
df = pd.DataFrame({'umi_origin': df1.umi_origin.unique()})
df['Tless_umi'] = df.umi_origin.apply(lambda x: x[1:5] + x[6:10] + x[11:15])          # cutting the T's out of the umi (at indexes 0,6,11,16)
df = df.sort_values('Tless_umi')
df['Serial'] = range(0, len(df), 1)
df['newSerial'] = np.nan
df['finalgroup'] = np.nan
j = 0
while j <12:     #initial clustering of similar neighbors in the data frame sorted alphabeticly from the j index to the end of the string
    df['sort_umi'] = df.Tless_umi.apply(lambda x: x[j:12])
    df = df.sort_values('sort_umi')
    df = df.reset_index(drop=True)
    newSerial = df.newSerial.tolist()
    Tlessumi = df.Tless_umi.tolist()
    Serial = df.Serial.tolist()
    for i in df.index:
        if i < len(df) - 1:
            score1 = lv.distance(Tlessumi[i], Tlessumi[i + 1])
            if score1 <= 1:
                if not np.isnan(newSerial[i]):
                    newSerial[i + 1] = newSerial[i]
                else:
                    newSerial[i] = Serial[i]
                    newSerial[i + 1] = newSerial[i]
    df.newSerial = newSerial
    j += 1
for i in df.index:
    if np.isnan(newSerial[i]):
        newSerial[i] = df.Serial[i]
df.newSerial = newSerial
df = df.sort_values('newSerial')
df = df.reset_index(drop=True)

def internalClust(l): # internal clustering function according to levenshtein distance - threshold- 1 nucleotide, method- complete linkage
    if len(l) == 1:
        return [0]
    Target = list(l)
    cond=[]
    for i in range(0, len(Target)):
        for j in range(i+1, len(Target)):
            cond = cond + [lv.distance(Target[i], Target[j])]
    #print (cond)
    hc_cut= hc.cut_tree(hc.complete(cond), height=3).flatten()
    #print(hc_cut,i,j)
    return hc_cut


newSerial = df.newSerial.tolist()
uniq_id = [(newSerial[0], 0)]
for i in df.index:  #storing the index of each initial cluster in tuples- tuple format: (original serial number, start index, end index)
    if i != 0:
        if newSerial[i] != newSerial[i - 1]:
            uniq_id[-1] += i,
            if(uniq_id[-1][2]-uniq_id[-1][1] + uniq_id[-1][0]> newSerial[i]):  #overlapping number prevention
                uniq_id.append((uniq_id[-1][2]-uniq_id[-1][1] + uniq_id[-1][0]+1, i))
            else:
                uniq_id.append((newSerial[i], i))
uniq_id[-1] += len(df),
finalgroup = []
Tlessumi = df['Tless_umi'].tolist()


for i in range(0, len(uniq_id)):
    temp = uniq_id[i]
    toClust= deepcopy(Tlessumi[temp[1]:temp[2] ])
    orig_index= [b[0] for b in sorted(enumerate(toClust), key=lambda i: i[1])]
    temp1 = internalClust(sorted(toClust))     #subset sent to clustering function
    j=0
    temp2= deepcopy(temp1)
    for i in orig_index:
        temp2[i]=temp1[j]
        j+=1
   # if len(temp1)>1: print(temp1)
    finalgroup.extend(np.array(temp2) + temp[0])    #assignment of new unique serial numbers to each cluster
df['finalgroup'] = finalgroup

df = df.sort_values('finalgroup')
df = df.reset_index(drop=True)
finalgroup = df.finalgroup.tolist()
uniq_id = [(finalgroup[0], 0)]
for i in df.index:
    if i != 0:
        if finalgroup[i] != finalgroup[i - 1]:
            uniq_id[-1] += i - 1,
            uniq_id.append((finalgroup[i], i))
uniq_id[-1] += len(df),

newUmi = []
umiOrigin = df['umi_origin'].tolist()
for i in range(0, len(uniq_id)):
    temp = uniq_id[i]
    if temp[2] + 1 > len(umiOrigin):
        newUmi.extend([umiOrigin[temp[1]]] * (temp[2] - temp[1]))
    else:
        newUmi.extend([umiOrigin[temp[1]]] * (temp[2] + 1 - temp[1]))
df['new_umi'] = newUmi
df.to_csv(args.o, columns=['umi_origin', 'new_umi', 'finalgroup'], sep='\t', index=False)
print(T.time() - time_s, T.process_time())
