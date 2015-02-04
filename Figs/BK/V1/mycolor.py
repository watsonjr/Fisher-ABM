
import pylab as pl
import numpy as np
import matplotlib.cm as cm
from matplotlib.colors import LogNorm, Normalize

def freecbar(rect,labelvec,cmap=cm.jet,label_pos='bottom',norm=None):
    """ Freefloating colorbar defined by rect=[lower x,lower y,len x, len y]
    Four or five labels.
    """
    if len(labelvec) == 5:
        tickvec = [0,25,50,75,100]
    else:
        #tickvec = [20,40,60,80]
        tickvec = [0,33,67,100]
    if label_pos=='left' or label_pos=='right': 
        cv = np.outer(np.ones(10),np.arange(0,1,0.01))
    else:
        cv = np.outer(np.arange(0,1,0.01),np.ones(10))

    cax = pl.axes(rect)
    if norm is None: norm = Normalize()
    pl.imshow(cv.transpose(),aspect='auto',origin='lower',cmap=cmap,norm=norm)
    if label_pos=='left' or label_pos=='right':
        cax.yaxis.set_ticks_position(label_pos)
        cax.set_xticks([])
        cax.set_yticks(tickvec)
        cax.set_yticklabels(labelvec)
    elif label_pos=='top':
        cax.xaxis.set_label_position(label_pos)
        cax.set_yticks([])
        cax.set_xticks(tickvec)
        cax.set_xticklabels(labelvec)
    else:
        cax.xaxis.set_label_position(label_pos)
        cax.set_yticks([])
        cax.set_xticks(tickvec)
        cax.set_xticklabels(labelvec)
    return cax


def WRY():
    import matplotlib as mpl
    cdict = {'red': (  (0.00, 1.0, 1.0),
                       (0.20, 1.0, 1.0),
                       (0.60, 0.3, 0.3),
                       (1.00, 1.0, 1.0)),
             'green': ((0.00, 1.0, 1.0),
                       (0.20, 0.5, 0.5),
                       (0.60, 0.0, 0.0),
                       (1.00, 1.0, 1.0)),
             'blue':  ((0.00, 1.0, 1.0),
                       (0.20, 0.4, 0.4),
                       (0.60, 0.0, 0.0),
                       (1.00, 0.4, 0.4))}
    return mpl.colors.LinearSegmentedColormap('GreenBlueRedYellow',cdict,256)

def GBRY():
    import matplotlib as mpl
    cdict = {'red': (  (0.00, 0.6, 0.6),
                       (0.20, 0.0, 0.0),
                       (0.40, 0.4, 0.4),
                       (0.50, 1.0, 1.0),
                       (0.60, 1.0, 1.0),
                       (0.80, 0.3, 0.3),
                       (1.00, 1.0, 1.0)),
             'green': ((0.00, 1.0, 1.0),
                       (0.20, 0.0, 0.0),
                       (0.40, 1.0, 1.0),
                       (0.50, 1.0, 1.0),
                       (0.60, 0.5, 0.5),
                       (0.80, 0.0, 0.0),
                       (1.00, 1.0, 1.0)),
             'blue':  ((0.00, 0.6, 0.6),
                       (0.20, 0.4, 0.4),
                       (0.40, 1.0, 1.0),
                       (0.50, 1.0, 1.0),
                       (0.60, 0.4, 0.4),
                       (0.80, 0.0, 0.0),
                       (1.00, 0.4, 0.4))}
    return mpl.colors.LinearSegmentedColormap('GreenBlueRedYellow',cdict,256)

def GBRY9010():
    import matplotlib as mpl
    cdict = {'red': (  (0.00, 0.6, 0.6),
                       (0.40, 0.0, 0.0),
                       (0.80, 0.4, 0.4),
                       (0.91, 1.0, 1.0),
                       (0.95, 1.0, 1.0),
                       (0.97, 0.3, 0.3),
                       (1.00, 1.0, 1.0)),
             'green': ((0.00, 1.0, 1.0),
                       (0.40, 0.0, 0.0),
                       (0.80, 1.0, 1.0),
                       (0.91, 1.0, 1.0),
                       (0.95, 0.5, 0.5),
                       (0.97, 0.0, 0.0),
                       (1.00, 1.0, 1.0)),
             'blue':  ((0.00, 0.6, 0.6),
                       (0.40, 0.4, 0.4),
                       (0.80, 1.0, 1.0),
                       (0.91, 1.0, 1.0),
                       (0.95, 0.4, 0.4),
                       (0.97, 0.0, 0.0),
                       (1.00, 0.4, 0.4))}
    return mpl.colors.LinearSegmentedColormap('GreenBlue90RedYellow10',cdict,256)

def GrGr():
    import matplotlib as mpl
    cdict = {'red': (  (0.00, 0.5, 0.5),
                       (1.00, 0.5, 0.5) ),
             'green': ((0.00, 0.5, 0.5),
                       (1.00, 0.5, 0.5) ),
             'blue':  ((0.00, 0.5, 0.5),
                       (1.00, 0.5, 0.5) ) }
    return mpl.colors.LinearSegmentedColormap('AllGrey',cdict,8)
