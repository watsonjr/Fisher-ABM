# -*- coding: utf-8 -*-
import numpy as np, scipy, networkx as nx, scipy.linalg as la,scipy.stats as stats
from math import *
from numpy import array
import bisect as bisct,random
from itertools import product as iprod
import itertools
import matplotlib.pyplot as plt,matplotlib.axes as axes
from matplotlib import rc
from mpltools import special as pltspecial, style as pltstyle
from collections import defaultdict


###======================== DATA ANALYSIS TOOLSET ======================###
### (Mostly for time series)


#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#                    INTERACTIVITY
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------



def wplot(typ,*args,**kwargs):
    save=kwargs.pop('save',0)
    hold=kwargs.pop('hold',False)
    title=kwargs.pop('title',False)
    legs=kwargs.pop('legends',False)
    invert=kwargs.pop('invert',0)
    if not legs:
        legs=kwargs.pop('legend',False)
    if not legs:
        legs=kwargs.pop('leg',False)
    lkw=kwargs.pop('legopt',{})
    if 'xlabel' in kwargs:
        plt.xlabel(kwargs.pop('xlabel'))
    if 'ylabel' in kwargs:
        plt.ylabel(kwargs.pop('ylabel'))
    if 'labels' in kwargs:
        labs=kwargs.pop('xlabel')
        plt.xlabel(labs[0])
        plt.ylabel(labs[1])
    if 'log' in kwargs:
        lg=kwargs.pop('log')
        if not isinstance(lg,str):
            if typ!= 'hist':
                kwargs['log']='xy'
            else:
                kwargs['log']=lg
        else:
            if 'x' in lg:
                plt.xscale('log')
            if typ!= 'hist':
                if 'y' in lg :
                    plt.yscale('log')
            else:
                if 'y' in lg:
                    kwargs['log']=1
                if 'x' in lg:
                    if min(args[0])<=0:
                        args=([ar for ar in args[0] if ar>0],)+args[1:]
                    if not 'bins' in kwargs or isinstance(kwargs['bins'],int):
                        kwargs['bins']=np.logspace(log10(min(args[0])),
                            log10(max(args[0])),kwargs.get('bins',100),10)
    elif typ=='hist' and min(args[0])>0 and max(args[0])/min(args[0])>10**3:
        kwargs['log']=1
    if 'xs' in kwargs :
        xs=kwargs.pop('xs')
        args=list(itertools.chain.from_iterable([[xs[:len(a)],a] for a in args]))
    if not args:
        return
    try:
        if typ=='plot':
            handle=plt.plot(*args,**kwargs)
        if typ =='hist':
            if not 'bins' in kwargs:
                kwargs['bins']=100
            if kwargs.pop('cleverbins',False):
                if isinstance(kwargs['bins'],int):
                    kwargs['bins']=min(kwargs['bins'],len(args[0])/20)
            kwargs['hold']=1
            handle=[plt.hist(a,**kwargs)[2][0] for a in args]
        if typ=='scatter':
            handle=plt.scatter(*args,**kwargs)
        if typ=='error':
            handle=plt.errorbar(*args,**kwargs)
        if typ=='errorfill':
            handle=pltspecial.errorfill(*args,**kwargs)
        if typ=='contour':
            handle=plt.contour(*args,**kwargs)
        if typ=='contourf':
            handle=plt.contourf(*args,**kwargs)
        if typ=='bar':
            handle=plt.bar(*args,**kwargs)
    except e:
        print "Cannot plot:", sys.exc_info()[0],e
        return
    if title:
        plt.title(title)
    if legs:
        plt.legend(handle,legs,**lkw)
    if invert:
        if 'x' in invert:
            plt.xlim(plt.xlim()[::-1])
        elif 'y' in invert:
            plt.ylim(plt.ylim()[::-1])
    if save:
        plt.savefig(save)
        plt.clf()
        plt.close()
    elif not hold:
        plt.show()
        return handle

def errorbar(*args,**kwargs):
    return wplot('error',*args,**kwargs)
def errorfill(*args,**kwargs):
    return wplot('errorfill',*args,**kwargs)
def plot(*args,**kwargs):
    return wplot('plot',*args,**kwargs)
def hist(*args,**kwargs):
    return wplot('hist',*args,**kwargs)
def scatter(*args,**kwargs):
    return wplot('scatter',*args,**kwargs)
def contour(*args,**kwargs):
    return wplot('contour',*args,**kwargs)
def contourf(*args,**kwargs):
    return wplot('contourf',*args,**kwargs)
def bar(*args,**kwargs):
    return wplot('bar',*args,**kwargs)

def mhist(*data,**kwargs):
    if len(data)>1:
        return [mhist(x,**kwargs) for x in data]
    else:
        data=data[0]
    kwargs.setdefault('log',0)
    if kwargs['log']:
        if min(data)<=0:
            data=[ar for ar in data if ar>0]
        if not 'bins' in kwargs or isinstance(kwargs['bins'],int):
            rg=kwargs.get('range', ( max(10**-50,min(data)) ,max(10**-50,max(data))) )
            kwargs['bins']=np.logspace(log10(rg[0] ),
                log10(rg[1]),kwargs.get('bins',100),10)
    bhist,bbins=np.histogram(data,bins=kwargs.pop('bins',200),density=kwargs.pop('normed',1),
        range=kwargs.get('range',None))
    wbins=(bbins[1:]-bbins[:-1])
    bbins=(bbins[:-1]+bbins[1:])/2
    if kwargs.get('cumulative',0):
        if kwargs['cumulative']>0:
            bhist=np.cumsum(bhist)
        else:
            bhist=np.sum(bhist)-np.cumsum(bhist)
    if kwargs.get('typ',False):
        wplot(kwargs['typ'],bbins,bhist,**kwargs)
    if kwargs.get('width',0):
        return bbins,wbins,bhist
    return bbins,bhist

def export(*args,**kwargs):
    pass

def launch_cmd():
    mycmd=MyCmd()
    mycmd.cmdloop()



#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#                    MATH
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------

def powlaw(expo,xmin,xmax):
    #distribution proportional to x**expo
    xx=expo+1
    if not xx:
        return exp(log(xmin) + log(xmax/xmin)*random.random())
    return (xmin**xx + (xmax**xx-xmin**xx)*random.random())**(1./xx)


def rint(nb):
    return int(round(nb))

def make_matrix(matrix):
    #if matrix is under dictionary form
    if hasattr(matrix,"__iter__") and not hasattr(matrix,"keys"):
        return matrix
    if hasattr(matrix.keys()[0],'__iter__'):
        keys=sorted(matrix)
        if hasattr(matrix.values()[0],'__iter__'):
            elsh=array(matrix.values()[0]).shape
        else:
            elsh=()
        sh=(sorted(matrix)[-1][0]+1,sorted(matrix,key=lambda e:e[1])[-1][1]+1)+elsh
        sh=tuple(array(sh))
        mat=np.zeros(sh)
        for i,j in sorted(matrix):
            mat[i,j]=matrix[(i,j)]
        return mat
    return [[matrix[i][j] for j in sorted(matrix[i])] for i in sorted(matrix)]

def aggregate(data,**kwargs):
    dic=kwargs.get('dic',False)
    method=kwargs.get('method','random')
    level=kwargs.get('level',4)
    matrix=array(make_matrix(data)).T
    try:
        keys=sorted(data.values()[0].keys())
        if method=='id':
            kbunches=[[keys[0]]]
            bunches=[[0]]
            for i,k in enumerate(keys):
                if not i : #jump 0
                    continue
                if k[:level]==kbunches[-1][-1][:level]:
                    kbunches[-1].append(k)
                    bunches[-1].append(i)
                else :
                    kbunches.append([k])
                    bunches.append([i])

    except:
        keys=range(0,len(data))
    #print bunches,kbunches
    if method =='random' or 'neigh' in method:
        nmat=np.zeros((matrix.shape[0]/level,matrix.shape[1]))
        for i in xrange(matrix.shape[0]/level):
            for z in xrange(level):
                if 'neigh' in method:
                    idx =i*level+z
                else :
                    idx= random.randint(0,matrix.shape[0]-1)
                nmat[i]+=matrix[idx]/level
                if method=='random':
                    np.delete(matrix,idx,0)
    if method=='id':
        nmat=np.zeros((len(bunches),matrix.shape[1]))
        for i,b in enumerate(bunches):
            for idx in b:
                nmat[i]+=matrix[idx]
    if dic:
        nmat={j:{k:nmat[i,l] for l,k in enumerate(sorted(matrix[j])) } for i, j in enumerate(sorted(matrix)) if i<len(nmat.shape[0])}

    return nmat.T



def make_returns(data,dic=True,tfirst=True):
    #dic if return result as dictionary
    #tfirst if time is the first variable in the data
    matrix=array(make_matrix(data))
    if not tfirst:
        matrix=np.transpose(matrix)

    rets=np.empty(matrix.shape-array([1,0]))
    it = np.nditer(matrix, flags=['multi_index'])
    while not it.finished:
        idx=it.multi_index
        try :
            rets[idx]
        except:
            it.iternext()
            continue
        #print "%d <%s>" % (it[0], it.multi_index),matrix[tuple(it.multi_index+array([1,0]))]
        d=matrix[tuple(idx+array([1,0]))]

        if it[0]>0 and d>0:
            rets[idx] = log(d/it[0])
        else :
            rets[idx]=0
        it.iternext()


    if not tfirst :
        rets = np.transpose(rets)
    if dic:
        rets={j:{k:rets[i,l] for l,k in enumerate(sorted(matrix[j])) } for i, j in enumerate(sorted(matrix)) if i<len(rets.shape[0])}

    #rets[~np.all(rets == 0, axis=1)]   remove all zero lines
    return rets

def sample_weight(tws):
    #array of weights
    ws=1.*np.copy(array(tws))
    totw=sum(ws)
    nb=len(ws)
    ws/=totw/nb
    ordr=np.argsort(ws)
    boxes={i:() for i in range(nb)}
    lastb=0
    poor=[i for i in range(nb) if ws[i]<1.]
    rich=[i for i in range(nb) if ws[i]>1.]
    while poor :
        i = poor.pop(0)
        if not rich:
            if not poor:
                break
            print poor
            raise Exception('Sample weight error')
        j=rich.pop(0)
        ws[j]-=(1-ws[i])
        boxes[i]=(ws[i],j)
        if ws[j]<1:
            poor.append(j)
        else:
            rich.append(j)
    return boxes

def sample_box(boxes):
    i = random.randint(0,len(boxes)-1)
    if boxes[i] and random.random()>boxes[i][0]:
            return boxes[i][1]
    return i

    #for i,o in enumerate(ordr):
        #if ws[o]>1:
            #boxes[i]=o
            #ws[o]-=1
            #remains.append(o)
        #else:
            #boxes[i]=(o,ws[o])
            #ws[o]=0
    #while remains:
        #if len(boxes[lastb])==1:
            #lastb+=1
            #continue
        #if ws[remains[0]]<boxes[lastb][1]:
            #remains.append(remains.pop(0))
            #continue
        #ws[remains[0]]-=boxes[lastb][1]
        #boxes[lastb]+=(remains[0],)
        #lastb +=1
    #return boxes


def make_density(data,**kwargs):
    mat=array(data)
    if mat.shape[1]==2 :
        #data is couples (x,y)
        mat=mat.T
    normed=kwargs.get('normed',False)
    nbin=kwargs.get('bins',100)
    logy=kwargs.get('logy',False)
    remove0=1-kwargs.get('include_zeroes',1)
    bintype='linear'
    xmin,xmax=np.min(mat[0]),np.max(mat[0])
    if kwargs.get('logx',False):
        numat=mat[0][mat[0]>10**(-100)]
        xmin,xmax=min(nozero(nonan(numat))),max(nonan(numat))
        bintype='log'
        bins=np.logspace(log10(xmin),log10(xmax),nbin,False)
        binw=[bins[i+1]-bins[i] for i in xrange(nbin-1)]
        binw.append(xmax-bins[-1])
    xspan=xmax-xmin
    try:
        bins
    except: #linear spacing case
        binw=xspan/nbin
        bins=np.linspace(xmin,xmax,nbin,False)
    binnage=[[] for i in xrange(nbin)]
    for x,y in mat.T :
        if x<xmin or x>xmax:
            continue
        if bintype=='linear':
            xbin=int(floor(float(x-xmin)/binw))
        else :
            xbin=bisct.bisect_left(bins,x)
        if xbin ==nbin: #maxvalue
            xbin=nbin-1
        if remove0 and abs(y)>10**(-40):
            binnage[xbin].append(y)
    res=array([stats.describe(i)[2:]+(min(i),max(i)) if i else stats.describe([0])[2:]+(0,0) for i in binnage])
    sspercen=scipy.stats.scoreatpercentile
    if kwargs.get('relative',1):
        quantile=array([array([sspercen(i,50),sspercen(i,50)-sspercen(i,5),sspercen(i,95)-sspercen(i,50)]) if i else array([0,0,0]) for i in binnage])
        res2=array([-res[:,-2]+res[:,0],res[:,-1]-res[:,0]])
    else:
        quantile=array([array([sspercen(i,50),sspercen(i,5),sspercen(i,95)]) if i else array([0,0,0]) for i in binnage])
        res2=array([res[:,-2],res[:,-1]])
    quantile=quantile.T
    if normed :
        if bintype=='linear':
            res[:,0]/=sum(res[:,0])*binw
        else :
            res[:,0]/=sum(np.dot(res[:,0],binw))
    return bins,res[:,0],res[:,1],res2,quantile[0],quantile[1:],array([len(i) for i in binnage])


def prune(data,**kwargs):
    #todo
    return



def standardize(orig,tfirst=0):
    if not hasattr(orig,'any') or not orig.any():
        return orig
    if tfirst:
        dat=np.copy(orig).T
    else :
        dat=np.copy(orig)
    for i in range(dat.shape[0]):
        std=np.std(dat[i][dat[i]!=0])
        avg=np.mean(dat[i][dat[i]!=0])
        if std:
            dat[i]=(dat[i]-avg)/std
    if tfirst:
        dat=dat.T
    return dat

#======================== FILTERS AND TOOLS ==========================

#======================== FILTERS AND TOOLS ==========================

#======================== FILTERS AND TOOLS ==========================

#======================== FILTERS AND TOOLS ==========================
def rang(l,lst=False):
    if lst:
        return range(len(l))
    return xrange(len(l))
def ranges(*ls):
    return iprod(*[rang(l) for l in ls])


def dicmap(dic,lst,pos=None):
    if pos == None:
        return [dic[i] for i in lst]
    else :
        return [dic[i[pos]] for i in lst]


def jfilter(matrix,rep=np.nan,sigmas=1.5):
    #remove abrupt jumps
    mat=np.copy(matrix)
    mas=np.ma.masked_array(matrix,np.isnan(matrix))
    test=0
    for p in xrange(mat.shape[0]):
        if not len(mas[p].compressed()):
            continue
        desc=scipy.stats.describe(mas[p].compressed())[2:4]
        #print desc
        for x in np.nditer(mat[p], op_flags=['readwrite']):
            if abs(x-desc[0])>sigmas*sqrt(desc[1]):
                x[...]=rep
                test+=1
    print 'Filtered out', float(test)/mat.size
    return mat

def rankfilter(matrix,per=(10,90),axis=0):
    pass

def nonan(mat):
    return mat[np.isnan(mat)==False]
def nozero(mat):
    return mat[mat!=0]
def noinf(mat):
    return mat[mat!=np.inf]

def demean(mat,axis=1):
    try:
        mat.shape[1]
        return mat-mat.mean(axis)
    except :
        return mat-mat.mean()


def mlog(matrix):
    mat=np.copy(matrix)
    for x in np.nditer(mat, op_flags=['readwrite']):
        if x>0:
            x[...]=np.log10(x)
        else :
            x[...]=np.nan
    return mat

def ranks (v):
    t = np.argsort(v)
    r = np.empty(len(v),int)
    r[t] = np.arange(len(v))
    for i in xrange(1, len(r)):
        if v[t[i]] <= v[t[i-1]]: r[t[i]] = r[t[i-1]]
    return r

def linreg(data):
    if len(data)==2:
        x,y=data
    else :
        x,y=array(list(enumerate(data))).T
    A = np.vstack([x, np.ones(len(x))]).T
    m, c = np.linalg.lstsq(A, y)[0]
    return m,c,[m*xx + c for xx in x]

#==============================================================

#==============================================================

#==============================================================

#==============================================================

def difs(matrix,axis=1):
    if axis==0:
        return matrix[1:]-matrix[:-1]
    if axis==1:
        return matrix[:,1:]-matrix[:,:-1]

def reldifs(matrix,axis=1):
    ds=difs(matrix,axis)
    cop=matrix[:,:-1]
    ds[cop!=0]/=cop[cop!=0]
    return ds


def rets(matrix,axis=1):
    mat = np.copy(matrix)
    whereAre0s = mat==0;
    mat[whereAre0s] = np.nan
    if axis==0:
        return np.log10(mat[1:]/mat[:-1])
    if axis==1:
        return np.log10(mat[:,1:]/mat[:,:-1])


def returnmat(matl='RCA',typ='prod',**kwargs):
    if not isinstance(matl,str):
        dat=matl
    else :
        if 'prod' in typ :
            dat=tradeprods[matl]
        if 'cnt' in typ or 'count' in typ:
            dat=tradecnts[matl]
        if 'filter' in kwargs:
            dat=proj_filter(dat,kwargs.get('filter'),typ)
    if not dat.any():
        print 'No data for return matrix with filter', kwargs.get('filter')
        return np.zeros(dat.shape-array([0,1])),np.zeros(dat.shape[0])
    rets=make_returns(dat,0,0)
    avgrets=np.zeros(rets.shape[0])
    for i in range(rets.shape[0]):
        avgrets[i]=np.mean(rets[i][rets[i]!=0])
    return rets, avgrets


#======================== MOMENTS ==========================


def lrgfrac(matrix,axis=0):
    #Fraction represented by the largest contribution
    return np.amax(matrix,axis)/np.sum(matrix,axis)

def mom(matrix,exp=1,axis=0):
    mat=np.ma.masked_array(matrix,np.isnan(matrix))
    #print mat[:,0],sum(mat[:,0]**exp)/len(mat[:,0])
    if len(matrix.shape)==1:
        if len(mat.compressed()):
            return sum(mat.compressed()**exp)/len(mat.compressed())
        else:
            return 0
    if axis==1: #moment over years
        return array([sum(mat[:,y].compressed()**exp)/len(mat[:,y].compressed())
            if len(mat[:,y].compressed())>0 else 0 for y in xrange(mat.shape[1])])
    if axis==0: #moment over prod/cnt
        return array([sum(mat[c,:].compressed()**exp)/len(mat[c,:].compressed())
            if len(mat[c,:].compressed())>0 else 0 for c in xrange(mat.shape[0])])

def quant(matrix,per=50,axis=0):
    mat=np.ma.masked_array(matrix,np.isnan(matrix))
    if axis==1:
        return array([scipy.stats.scoreatpercentile(mat[:,y].compressed(), per) for y in xrange(mat.shape[1]) ])
    if axis==0:
        return array([scipy.stats.scoreatpercentile(mat[c,:].compressed(), per) for c in xrange(mat.shape[0]) ])

def rankwid(matrix,per=(16,84),axis=0):
    return (quant(matrix,per[1])-quant(matrix,per[0]))/2


def binomial(n, k):
    if not 0 <= k <= n:
        return 0
    if k == 0 or k == n:
        return 1
    P = k+1
    for i in xrange(k+2, n+1):
        P *= i
    return P//factorial(n-k)

def cumul(matrix,rank,axis=0,reduced=0):
    moments=array([mom(matrix,exp+1,axis) for exp in xrange(rank)]).T
    cumulants=np.empty_like(moments)
    binom={}
    for z in xrange(rank):
        for x in xrange(z):
            binom[(z,x)]=int(binomial(z,x))
    for yr in xrange(matrix.shape[axis]):
        cumulants[yr]=moments[yr]
        for z in xrange(rank):
            for x in xrange(z):
                cumulants[yr,z]-=cumulants[yr,x]*moments[yr,z-1-x]*binom[(z,x)]
        if reduced:
            cumulants[yr] = array([cumulants[yr,z]/moments[yr,0]**(z+1) for z in range(rank)])
    return cumulants.T

#=============== SELF CORRELATION

def selfcorr(mat,**kwargs):
    npd,tim=mat.shape
    maxdelay=kwargs.get('delay',tim/2)
    mean,var=cumul(mat,2)
    mean=mean.reshape(len(mean),1)
    var=var.reshape(len(var),1)
    newmat=(mat-mean)/np.sqrt(var)
    newmat[np.isnan(newmat)]=0
    scor=np.zeros((npd,maxdelay))
    rep=np.zeros((npd,maxdelay))
    for delay in xrange(maxdelay):
        for yr in xrange(tim-delay):
            nyr=yr+delay
            scor[:,delay]+=newmat[:,nyr]*newmat[:,yr]
            rep[:,delay]+=np.logical_and(newmat[:,nyr], newmat[:,yr])
    scor[scor!=0]/=rep[rep!=0]
    return scor
