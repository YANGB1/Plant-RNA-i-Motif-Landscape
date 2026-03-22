import os
import time
import sys
import codecs
import argparse
import gc
import math
import numpy as np
import pandas as pd
import scipy
from scipy import stats
from scipy.stats.mstats import kruskalwallis
from scipy.stats import spearmanr
from scipy.stats import pearsonr

environ={}
sig_dic={}
with open("environment_variable.txt","r") as f:
    sig=0
    for line in f:
        q=line.strip().split("\t")
        if sig==0:
            for i in range(len(q)):
                environ[q[i]]=[]
                sig_dic[str(i)]=q[i]
                sig=sig+1
        else:
            for i in range(len(q)):
                environ[sig_dic[str(i)]].append(q[i])

density={}
sig_dic={}
with open("iMdensity_statistic.txt","r") as f:
    sig=0
    for line in f:
        q=line.strip().split("\t")
        if sig==0:
            for i in range(len(q)):
                density[q[i]]=[]
                sig_dic[str(i)]=q[i]
                sig=sig+1
        else:
            for i in range(len(q)):
                density[sig_dic[str(i)]].append(q[i])

density_new={}
for i in density.keys():
    if i !="Species":
        density_new[i]=[]

environ_new={}
for i in environ.keys():
    if i !="Label":
        environ_new[i]=[float(j) for j in environ[i]]

for i in environ["Label"]:
    for j in range(len(density["Species"])):
        if density["Species"][j]==i:
            for n in density_new.keys():
                density_new[n].append(float(density[n][j]))

out=open("PCC_density_enviromental.txt","w")
out.write("environmental factors"+"\t")
for j in density_new.keys():
    out.write(j+"\t"+"p"+"\t")
out.write("\n")
for i in environ_new.keys():
    out.write(i+"\t")
    for j in density_new.keys():
        pcc, p =pearsonr(environ_new[i], density_new[j])
        #pcc, p =spearmanr(environ_new[i], density_new[j])
        out.write(str(pcc)+"\t"+str(p)+"\t")
    out.write("\n")
out.close()
        

