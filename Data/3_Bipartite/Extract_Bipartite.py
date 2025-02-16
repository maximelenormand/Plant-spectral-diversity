# Imports
import os 
import random
import shutil
import copy
import operator
import sys
import numpy

# Parameter
wd='Data/'

# Read HSP
hspfile=open(wd + '2_Spectral/CoMat_Spectral.csv')
col=hspfile.readline().rstrip('\n\r').split(';')
hsp={}
for line in hspfile:
    attr=line.rstrip('\n\r').split(';')
    hsp[attr[0]]=tuple(float(s) for s in attr[1:])

# Read Vege
vegefile=open(wd + '1_Plant/CoMat_Plant.csv')
col=vegefile.readline().rstrip('\n\r').split(';')
vege={}
for line in vegefile:
    attr=line.rstrip('\n\r').split(';')
    vege[attr[0]]=tuple(float(s) for s in attr[1:])

# Loop
bip=open(wd + '3_Bipartite/Bipartite.csv','w')    
f=bip.write("Plant;Spectral;a;b;c;A;B;C\n")  

i=0
for key1, val1 in vege.items():
    
    i=i+1
    print ('Nb plants: ' + str(i))
    
    for key2, val2 in hsp.items():
        
        X=numpy.array(val1)
        Y=numpy.array(val2)
        
        x=(X>0).astype(float)
        y=(Y>0).astype(float)
        
        a=sum(tuple(0.5*(x+y-abs(x-y))))
        b=sum(x)-a
        c=sum(y)-a
        
        A=sum(tuple(0.5*(X+Y-abs(X-Y))))
        B=sum(X)-A
        C=sum(Y)-A
        
        f=bip.write(str(key1))
        f=bip.write(';')
        f=bip.write(str(key2))
        f=bip.write(';')
        f=bip.write(str(a))
        f=bip.write(';')
        f=bip.write(str(b))
        f=bip.write(';')
        f=bip.write(str(c))
        f=bip.write(';')
        f=bip.write(str(A))
        f=bip.write(';')
        f=bip.write(str(B))
        f=bip.write(';')
        f=bip.write(str(C))
        f=bip.write('\n')

bip.close()









