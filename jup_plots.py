import scipy.interpolate
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt, cm, mlab
from scipy.stats import beta

def contour_Q_ns(xi,yi,zi,lev,col=None,hexgrid=None,orient=None,boundlev=None,logclabel=True,clabel=r'$\log_{10}$ Stochastic Fraction',outline=True):
    # Make contour plot
    cmap = cm.get_cmap('magma')
    if col is None:
        col = [cmap(a) for a in np.linspace(0,1,len(lev))]
    if not boundlev is None:
        xib,yib,zib,levb = boundlev
        plt.contour(xib,yib,zib,levels=[levb],colors='white',linestyles='-.',linewidths=1.5)
    if outline:
        plt.contour(xi,yi,zi,levels=lev,colors='gray',linewidths=0.5)
    plt.contourf(xi,yi,zi,levels=lev,colors=col,alpha=1)
    if not hexgrid is None and (len(hexgrid) < 5 or hexgrid[5] == True):
        x,y,z,gsize = hexgrid
        plt.hexbin(x,y,np.log10(z),gridsize=gsize,cmap=cm.magma,bins=None)
    plt.xlabel(r"$\log_{10} \, Q$")
    plt.ylabel(r"$n_s$")

    # Add colorbar
    if orient is None:
        cb = plt.colorbar()
    else:
        cb = plt.colorbar(orientation=orient)
    cb.set_label(clabel)
    cbticks = cb.get_ticks()
    cblim = cb.get_clim()
    cb.set_ticks(cbticks)

    # Handle tick labels on colorbar
    if logclabel == True:
        cb.set_ticklabels([r"${%.1f}$" % x for x in np.log10(cbticks)])
    elif logclabel == False:
        cb.set_ticklabels([r"${%.2f}$" % x for x in (cbticks)])
    elif logclabel == -1:
        cb.set_ticklabels([r"${%.1f}$" % x for x in np.log10(1-cbticks)])
    else:
        cb.set_ticklabels([r"$_{10}{%.1f}$" % (x) for x in np.log10(cbticks[:logclabel+1])]+[r"${%.2f}$" % x for x in (cbticks[logclabel+1:])])

    return plt.gcf()

def binomialUpperBound(frac,fieldname,invflag=False,jeffreys=True):
    if invflag:
        z = list(1-frac[fieldname])
    else:
        z = list(frac[fieldname])
    counts = list(frac['len'])
    nstoch = list(frac['len']*z)

    if jeffreys:
        z = [beta.ppf(0.95,0.5+ns,0.5+(c-ns)) for c,ns in zip(counts,nstoch)]
    else:
        z = [beta.ppf(0.95,1+ns,1+(c-ns)) for c,ns in zip(counts,nstoch)]

    return z

def binomialLowerBound(frac,fieldname,invflag=False,jeffreys=True):
    if invflag:
        z = list(1-frac[fieldname])
    else:
        z = list(frac[fieldname])
    counts = list(frac['len'])
    nstoch = list(frac['len']*z)

    if jeffreys:
        z = [beta.ppf(0.05,0.5+ns,0.5+(c-ns)) for c,ns in zip(counts,nstoch)]
    else:
        z = [beta.ppf(0.05,1+ns,1+(c-ns)) for c,ns in zip(counts,nstoch)]

    return z

def poissonUpperBound(frac,fieldname,invflag=False):
    if invflag:
        z = list(1-frac[fieldname])
    else:
        z = list(frac[fieldname])
    length = list(frac['len'])
    counts = list(frac['len'])
    nstoch = list(frac['len']*z)

    # Number by which to multiply counts [-,1,2,3,...] in order to get 95% upper bound
    # integral(@(lam) pois(lam,0),0,offset[0]) = 0.95
    # integral(@(lam) pois(lam,k),0,k*offset[k]) = 0.95 for k > 0
    offset=np.power(10.0,[0.47,0.67,0.49,0.41,0.36,0.32,0.29,0.27,0.26,0.24,0.23])

    # If there are zero events, just take the offset value
    z = [(offset[0]/counts[j] if nstoch[j] < 0.00005 else z[j]) for j in range(0,len(z))]

    # If there are nonzero events, multiply by the offset
    for i in range(1,len(offset)):
        z = [(z[j]*offset[i] if (nstoch[j] > (i-0.5) and nstoch[j] < (i+0.5)) else z[j]) for j in range(len(z))]

    # If events is greater than 10, use an approximation to get the 95% upper bound
    # For large N, Poisson becomes approximately Gaussian, and 95% of probability is
    # contained within 2*sigma of mean
    z = [(z[j]*(1+2/np.sqrt(nstoch[j]))) if nstoch[j] > 10.5 else z[j] for j in range(len(z))]

    return z

def getBinCenters(xedge,yedge):
    x = np.pad((xedge[:-1]+xedge[1:])/2,(1,1),'constant')
    y = np.pad((yedge[:-1]+yedge[1:])/2,(1,1),'constant')
    x,y = np.meshgrid(x,y)
    x = np.transpose(x).reshape((1,-1)).tolist()[0]
    y = np.transpose(y).reshape((1,-1)).tolist()[0]
    return x,y
