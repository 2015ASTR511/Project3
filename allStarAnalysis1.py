import matplotlib.pyplot as plt
import numpy as np
import pyfits
from astroML.stats import binned_statistic_2d


def funcStdDev(vector):
    return np.sqrt(np.var(vector))


def plotPanel(axes, xVec, yVec, zVec, CCstat="", CCtype="", xMin="", xMax="", yMin="", yMax="", xLabel="", yLabel="", title=""): 

    ## axes, etc.
    # axes limits
    if (xMin==""): xMin = np.min(xVec)
    if (xMax==""): xMax = np.max(xVec)
    axes.set_xlim(xMin, xMax)
    if (yMin==""): yMin = np.min(yVec)
    if (yMax==""): yMax = np.max(yVec)
    axes.set_ylim(yMin, yMax)
    # axes labels
    if (xLabel==""): xLabel = "X"
    if (yLabel==""): yLabel = "Y" 
    axes.set_xlabel(xLabel, fontsize=14)
    axes.set_ylabel(yLabel, fontsize=14)
    # title 
    axes.set_title(title)
    # make bigger ticmark labels
    for tick in axes.xaxis.get_major_ticks():
        tick.label.set_fontsize(12) 
    for tick in axes.yaxis.get_major_ticks():
        tick.label.set_fontsize(12) 

    ## bin
    N, xedges, yedges = binned_statistic_2d(xVec, yVec, zVec,'count', bins=100)
    if (CCstat==""):
        Zmean, xedges, yedges = binned_statistic_2d(xVec, yVec, zVec,'mean', bins=100)
    else:
        Zmean, xedges, yedges = binned_statistic_2d(xVec, yVec, zVec,funcStdDev, bins=100)

    ## plot
    cmap_multicolor = plt.cm.jet
    cmap_multicolor.set_bad('w', 1.)
    plt.imshow(Zmean.T, origin='lower',
           extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
           aspect='auto', interpolation='nearest', cmap=cmap_multicolor)
    # special cases
    if (CCtype=="FeH"):
        cb = plt.colorbar(ticks=np.arange(-2.5, 1, 1), pad=0.16,
                  format=r'$%.1f$', orientation='horizontal')
        cb.set_label(r'$\mathrm{mean\ [Fe/H]\ in\ pixel}$', fontsize=12)
        plt.clim(-2.5, 0.5)
    if (CCtype=="logg"):
        cb = plt.colorbar(ticks=np.arange(0, 4, 1), pad=0.16,
                  format=r'$%.1f$', orientation='horizontal')
        cb.set_label(r'$\mathrm{mean\ log(g)\ in\ pixel}$', fontsize=12)
        plt.clim(0, 4)
    if (CCtype=="vRad"):
        cb = plt.colorbar(ticks=np.arange(0, 150, 25), pad=0.16,
                  format=r'$%.1f$', orientation='horizontal')
        cb.set_label(r'$\mathrm{\ vel. disp. (km/s) \ in\ pixel}$', fontsize=12)
        plt.clim(0, 150)

    # density contours over the colors
    levels = np.linspace(0, np.log10(N.max()), 7)[2:]
    plt.contour(np.log10(N.T), levels, colors='k',
                  extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])

 


##################################################################### 
### read APOGEE parameters
# from SDSS Data Release 12 
data = pyfits.open('allStar-v603.fits')[1].data[::1]  # 163,278

# some ad hoc quality cuts 
OKcond = ((data['LOGG'] > -10) & (data['PARAM_M_H'] > -10))
OKcond = (OKcond & (data['J'] > -10) & (data['J'] < 30))
OKcond = (OKcond & (data['H'] > -10) & (data['H'] < 30))
OKcond = (OKcond & (data['K'] > -10) & (data['K'] < 30))
# and some more severe cuts
# we don't want very bright objects as they seem outliers (~500) 
OKcond = (OKcond & (data['K'] > 6.0))
### ZI: but the above is so sloppy and unbecoming of a good survey 
###     astronomer - we need to (also) use FLAGS!!!

# full sky; 101,399
title = "APOGEE, full sky"

# for investigating variations with the sky position
if (0):
    # high latitudes
    OKcond = (OKcond & (data['GLAT'] > 20.0))  # 19,249
    title = "APOGEE, b>20"

if (0):
    # low latitudes
    OKcond = (OKcond & (data['GLAT'] < 5.0) & (data['GLAT'] > -5.0))  # 43,002
    title = "APOGEE, |b|<5 deg."

## select good data and make vectors for convenience 
dataOK = data[OKcond]   
J = dataOK['J']
H = dataOK['H']
K = dataOK['K']
JK = J-K
logg = dataOK['LOGG']
Teff = dataOK['TEFF']
FeH = dataOK['PARAM_M_H']
alphaFe = dataOK['PARAM_ALPHA_M']
gLat = dataOK['GLAT']
W2 = dataOK['WISE_4_5']
KW2 = K-W2
vRad = dataOK['VHELIO_AVG']



### plot

# Create figure and subplots
fig = plt.figure(figsize=(8, 8))
# this work well in *.py version but not so well in ipython notebook
fig.subplots_adjust(wspace=0.25, left=0.1, right=0.95, bottom=0.12, top=0.95)


### all these color-coded by the mean FeH 
# K vs. J-K color-magnitude diagram 
axes = plt.subplot(221)
plotPanel(axes, JK, K, logg, CCtype="logg", xMin=0.0, xMax=4.0, yMin=14.01, yMax=6.0,
          xLabel="J-K", yLabel="K", title=title)

# logg vs. Teff, color-coded by mean FeH 
axes = plt.subplot(222, xticks=[3000, 3500, 4000, 4500, 5000, 5500, 6000])
plotPanel(axes, Teff, logg, FeH, CCtype="FeH", xMin=5500, xMax=3500, yMin=4.2, yMax=-0.5,
          xLabel=r'$\mathrm{T_{eff}}$', yLabel="log(g)")

# [alpha/Fe] vs. [Fe/H], color-coded by logg
axes = plt.subplot(223)
plotPanel(axes, FeH, alphaFe, logg, CCtype="logg", xMin=-1.8, xMax=0.65, yMin=-0.2, yMax=0.5,
          xLabel=r'$\mathrm{[Fe/H]}}$', yLabel=r'$\mathrm{[\alpha/Fe]}}$')


# [alpha/Fe] vs. [Fe/H], color-coded by radial velocity 
axes = plt.subplot(224)
plotPanel(axes, FeH, alphaFe, vRad, CCstat="stdev", CCtype="vRad", xMin=-1.8, xMax=0.65, yMin=-0.2, yMax=0.5,
          xLabel=r'$\mathrm{[Fe/H]}}$', yLabel=r'$\mathrm{[\alpha/Fe]}}$')

plt.savefig('./apogee1.png')
plt.show() 

