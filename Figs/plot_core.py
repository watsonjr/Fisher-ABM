import numpy as np
import scipy.special as SS
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from math import *
import os
from datatools import *
import cmath
import scipy.optimize as opt

path="../Data/"

#=========== Extraction of constants from julia code ==================


constants={}
def set_constants(filename):
    if not ".dat" in filename:
        filename=filename+".dat"
    fcst=open(os.path.join(path,filename ),"r")
    #constants={}
    global constants
    for l in fcst:
        line=l.split()
        if line and '#' in line[0]:
            continue
        if len(line)>1:
            exec("{} = {}".format(line[0],line[1]))
            exec("constants['{}'] = {}".format(line[0],line[1]))
    fcst.close()

def load_data(filename):
    if not ".npy" in filename:
        filename=filename+".npy"
    data= np.load(os.path.join(path,filename ))
    #get headers
    headers=[]
    fcst=open(os.path.join(path,filename.replace("npy","dat" )),"r")
    for l in fcst:
        line=l.split()
        if line and '#' in line[0]:
            line=[x.strip().replace('#','') for x in line ]
            headers=[ x for x in line if x]
            break
    datadict={}
    dim=len(data.shape)
    data=np.transpose(data,[dim-1]+range(dim)[:-1] ) #transpose to get dataset as first index
    for i,h in enumerate(headers):
        try:
            datadict[h]=data[i]
        except:
            print "ERROR in load_data: ", headers, data.shape[0]
    return datadict

def get_constants(**kwargs):
    args=[ kwargs[i] if i in kwargs else constants[i] for i in 
        ("PC_rp","PC_f","PC_q","PF_sig" ,"PF_n" ,"GRD_nx","GRD_dx","PC_v","PS_p","PS_n","PC_n")]
    return args
