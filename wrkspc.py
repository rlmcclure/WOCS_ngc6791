#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 13:00:15 2020

@author: rlmcclure
"""

import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from astropy.io import fits
import scipy.signal as sig
from scipy import stats
import lightkurve as lk
import pandas as pd
#import eleanor.Visualization as vis
#RUN functions.py 
#%% k2mosaic steps http://k2mosaic.geert.io/
'''
For each quarter in each of the four years we have to do 
$ k2mosaic tpflist Q[#] [CH#] > tpflist.txt
$ k2mosaic mosaic tpflist.txt 

if you get stopped add --cadence #####..##### between mosaic and tpf 


e.g. $ k2mosaic mosaic 30771..33935 tpflist.txt 

Channels are 81, 53, 29, 1
so this has to occur for 
Q4, Q8, Q12, Q16 for channel 81
Q3, Q7, Q11, Q15 for channel 53 or 1
Q2, Q6, Q10, Q14, Q17 for channel 29 
Q1, Q5, Q9,  Q13 for channel 1 or 53

We start with just channel 81 for our initial round of analysis 

to re
'''
#%%
#filesPath = '../../../../../../Volumes/rmathieu/NGC6791/'#/k2mosaic-master/'
#%
#targFile = 'kepler/kplr002297232-2010078095331_lpd-targ.fits'
#tagFilePath = filesPath+targFile
#hdul = fits.open(targFilePath)
#%
#hdul.info()
#%
#primary = hdul['PRIMARY'].data
#targettables = hdul['TARGETTABLES'].data
#aperture = hdul['APERTURE'].data
#%%
#filesPath = '/Volumes/rmathieu/NGC6791/k2mosaic-master/mosaic-ch81'
filesPath = '/Users/kookoo2052/HOLD/mosaic-ch81'
#%%
filesList = [f for f in os.listdir(filesPath) if (f.endswith(".fits") and os.path.isfile(os.path.join(filesPath, f)))]
filesList.sort()
#%%
fistarList = [f for f in os.listdir(filesPath) if (f.endswith(".fistar") and os.path.isfile(os.path.join(filesPath, f)))]
fistarList.sort()
#%%Extrac the stamp
filesIDs = np.empty([len(filesList)])
stampBox = np.empty([len(filesList),200,200])

#%%Stamp loc [723:923,20:220]
for idx in np.arange(len(filesList)):
    
    #Create list of file IDs
    filesIDs[idx] = int(filesList[idx][-10:-5])

    #Create the full path string
    fileStr = filesList[idx]
    filePath = filesPath+'/'+fileStr
    
    #Pull imag for each file
    [hdr,img] = getvar(filePath,'IMAGE',verbose=0)
    stampBox[idx]=img[723:923,20:220]

#%%
ids = [72.0, 41.0, 81.0, 24.0, 39.0, 26.0, 10.0, 69.0, 85.0, 82.0]
for id in ids:
    Xs = []
    Ys = []
#%% FISTAR output tracking y:23:43,x:150:170
xmin=20+156
xmax=20+161
ymin=723+30
ymax=723+35
siz = len(fistarList)#np.size(stampBoxNonNan,0)
ID = np.empty(siz)
X = np.empty(siz)
Y = np.empty(siz)

for fidx in np.arange(len(fistarList)):
    df = fistar2df(filesPath+'/'+fistarList[fidx])
    if df.empty != 1:
        dfSS = df[(df['X'] >xmin) & (df['X'] <xmax) & (df['Y'] >ymin) & (df['Y'] <ymax)]
        
        ID[fidx] = int(dfSS.index[0])
        X[fidx]  = float(dfSS['X'][dfSS.index[0]])
        Y[fidx]  = float(dfSS['Y'][dfSS.index[0]])
#%%
plt.title('fistar X Drift')
plt.scatter(np.arange(len(X)),X,marker='.',c='k',s=.5)
plt.ylim(177.7,178.1)
plt.ylabel('fistar Centroid Location')
plt.xlabel('Frame Number, Cadence')
plt.savefig('fistarXDrift.png',dpi=200,bbox_inches='tight')
plt.show()
#%
plt.title('fistar Y Drift')
plt.scatter(np.arange(len(Y)),Y,marker='.',c='k',s=.5)
plt.ylim(755.1,755.5)
plt.ylabel('fistar Centroid Location')
plt.xlabel('Frame Number, Cadence')
plt.savefig('fistarYDrift.png',dpi=200,bbox_inches='tight')
plt.show()
#%% m soares code from 5/13/2020
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)

from astropy.io import fits 
from astropy.table import Table
from astroquery.mast import Catalogs
from astroquery.vizier import Vizier

v = Vizier()
from astropy import coordinates
from astropy import units as u

#Live dangerously
import warnings
warnings.filterwarnings("ignore")


# interface with astrometry package instead of opening this fits file
tbl = fits.open('/Users/msoares/Desktop/corr.fits')[1].data





#%%            Xs.append(dfSS.loc[id]['X'])
#            Ys.append(dfSS.loc[id]['Y'])
#            
    #%plot some shit to figure out what this is showing
    idN = int(id)
    xBest = np.median(Xs)
    yBest = np.median(Ys)
#    xBest = stats.mode(np.round(Xs[np.argwhere(Xs>1)]))[0][0][0]
#    yBest = stats.mode(np.round(Ys[np.argwhere(Ys>1)]))[0][0][0]
    #%
    plt.title('Returned Locations for Object %s'%(idN))
    plt.scatter(Xs,Ys)
    #plt.scatter(Xs[np.argwhere(Xs>1)],Ys[np.argwhere(Ys>1)])
    plt.xlabel('X Location Values')
    plt.ylabel('Y Location Values')
    plt.savefig('Locations%s.png'%(idN),dpi=200,bbox_inches='tight')
    plt.show()
    #%
    plt.title('Y Drift for Object %s, Selected Mode Value'%(idN))
#    plt.scatter(np.arange(len(np.argwhere(Ys>1))),Ys[np.argwhere(Ys>1)])
    plt.scatter(np.arange(len(Ys)),Ys)
    plt.xlabel('Cadence Frame')
    plt.ylabel('Y Location Value')
    plt.ylim(yBest-.5,yBest+.5)
    plt.savefig('YDrift%sZoom.png'%(idN),dpi=200,bbox_inches='tight')
    plt.show()
    
    plt.title('Y Drift for Object %s'%(idN))
#    plt.scatter(np.arange(len(np.argwhere(Ys>1))),Ys[np.argwhere(Ys>1)])
    plt.scatter(np.arange(len(Ys)),Ys)
    plt.xlabel('Cadence Frame')
    plt.ylabel('Y Location Value')
    plt.savefig('YDrift%s.png'%(idN),dpi=200,bbox_inches='tight')
    plt.show()
        #%
    plt.title('X Drift for Object %s, Selected Mode Value'%(idN))
    plt.xlabel('Cadence Frame')
    plt.scatter(np.arange(len(Xs)),Xs)
    plt.ylabel('X Location Value')
    plt.ylim(xBest-.5,xBest+.5)
    plt.savefig('XDrift%sZoom.png'%(idN),dpi=200,bbox_inches='tight')
    plt.show()
    
    plt.title('X Drift for Object %s'%(idN))
#    plt.scatter(np.arange(len(np.argwhere(Xs>1))),Xs[np.argwhere(Xs>1)])
    plt.scatter(np.arange(len(Xs)),Xs)
    plt.xlabel('Cadence Frame')
    plt.ylabel('X Location Value')
    plt.savefig('XDrift%s.png'%(idN),dpi=200,bbox_inches='tight')
    plt.show()
#%% remove nan from field NEED TO TRACK INDEXES OF THIS
non_nan_ind = ~np.isnan(stampBox)
stampBoxNonNan = np.resize(stampBox[non_nan_ind],
                           [np.shape(stampBox)[0]-len(np.where(non_nan_ind==0)),
                           np.shape(stampBox)[1],np.shape(stampBox)[2]])
#%% get medians
meanField = np.mean(stampBoxNonNan,0)
medianField = np.mean(stampBoxNonNan,0)
yMax,xMax = np.where(medianField == np.amax(medianField))
#%% median subtracted frames
stampMedSub = stampBoxNonNan-medianField
#%%
for nn in np.arange(0,len(stampBoxNonNan),100):
#    quick_imshow_savepngmv(stampBoxNonNan[nn],'Cadence '+str(cadList[nn]),'cadStamp',nn,'terrain',0)
    quick_imshow_savepngmv(np.log(np.abs(stampBoxNonNan[nn,0:150,50:220])),'Field '+str(nn),'cadStampSUBFIELDf',nn,'binary_r',0)
    #%%
y,x = [100,30]
pltimg = stampMedSub
plt.plot(pltimg[:,x,y],label='data')
plt.legend()
plt.show()

fig = plt.figure()
ax = plt.subplot(221)
ax.plot(sig.detrend(pltimg[:,x,y]),label='sig.detrend(type=linear)')
ax = plt.subplot(222)
ax.plot(sig.detrend(pltimg[:,x,y],type='constant'),label='sig.detrend(type=constant)')
plt.show()
#%%
nn=80
quick_imshow_savepngmv(stampBoxNonNan[nn],'Cadence '+str(filesIDs[nn]),'cadStamp',nn,'gist_heat_r',1)
quick_imshow_savepngmv(medianField,'Median Field','medianField',0,'gist_heat_r',1)
#%%
quick_imshow_savepngmv(stampMedSub[nn],'Median Subtracted Cadence '+str(filesIDs[nn]),'cadStampSub',nn,'BrBG',1)
#%%

quick_imshow_savepngmv(img,'Full Channel 81 Field','fullField',0,'gist_heat_r',1)
#%% pixel by pixel plot from eleanor.Visualization tess package
FOI = stampMedSub      
def pixel_by_pixel(box, colrange=None, rowrange=None, cmap='viridis', mask=None,
                    xlim=None,ylim=None, color_by_pixel=True, 
                    freq_range=[1/20., 1/0.1]):
    
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    """
    Creates a pixel-by-pixel light curve using the corrected flux.
    Contribution from Oliver Hall.
    Parameters
    ----------
    colrange : np.array, optional
         A list of start column and end column you're interested in
         zooming in on.
    rowrange : np.array, optional
         A list of start row and end row you're interested in zooming
         in on.
    cmap : str, optional
         Name of a matplotlib colormap. Default is 'viridis'.
    #data_type : str, optional
         The type of flux used. Either: 'raw', 'corrected', 'amplitude',
         or 'periodogram'. If not, default set to 'corrected'.
    mask : np.array, optional
         Specifies the cadences used in the light curve. If not, default
         set to good quality cadences.
    xlim : np.array, optional
         Specifies the xlim on the subplots. If not, default is set to 
         the entire light curve.
    ylim : np.array, optional
         Specifies the ylim on the subplots, If not, default is set to 
         the entire light curve flux range.
    color_by_pixel : bool, optional
         Colors the light curve given the color of the pixel. If not,
         default is set to True.
    #freq_range : list, optional
         List of minimum and maximum frequency to search in Lomb Scargle
         periodogram. Only used if data_type = 'periodogram'. If None,
         default = [1/20., 1/0.1].
    """
    if colrange is None:
        colrange = [0, np.shape(box)[-1]]
    
    if rowrange is None:
        rowrange = [0, np.shape(box)[-2]]
        
    nrows = int(np.round(colrange[1]-colrange[0]))
    ncols = int(np.round(rowrange[1]-rowrange[0]))
    
      
    fig = plt.figure(figsize=(20,8))
    outer = gridspec.GridSpec(1,2, width_ratios=[1,4])
    
    inner = gridspec.GridSpecFromSubplotSpec(ncols, nrows, hspace=0.1, wspace=0.1,
                                             subplot_spec=outer[1])
    
    i, j = rowrange[0], colrange[0]
    
    #    if mask is None:
    #        q = self.obj.quality == 0
    #    else:
    #        q = mask == 0
    
    
    ## PLOTS TARGET PIXEL FILE ##
    ax = plt.subplot(outer[0])
    
    c = ax.imshow(np.nanmedian(box[:,rowrange[0]:rowrange[1],colrange[0]:colrange[1]],0),cmap=cmap)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.15)
    plt.colorbar(c, cax=cax, orientation='vertical')
    
    ## PLOTS PIXEL LIGHT CURVES ##
    for ind in range( int(nrows * ncols) ):
        print(ind)
        ax = plt.Subplot(fig, inner[ind])
    
        flux = box[:,i,j]
        y = flux
        x = np.arange(len(y))
    #        time = self.obj.time
    #        corr_flux = self.obj.corrected_flux(flux=flux)
    
    #        if data_type.lower() == 'corrected':
    #            y = corr_flux[q]/np.nanmedian(corr_flux[q])
    #            x = time[q]
    
    #        elif data_type.lower() == 'amplitude':
    #            lc = lk.LightCurve(time=time, flux=corr_flux)
    #            pg = lc.normalize().to_periodogram()
    #            x = pg.frequency.value
    #            y = pg.power.value
    #    
    #        elif data_type.lower() == 'raw':
    #            y = flux[q]/np.nanmedian(flux[q])
    #            x = time[q]
    #        
    #        elif data_type.lower() == 'periodogram':
    #            freq, power = LombScargle(time, corr_flux).autopower(minimum_frequency=freq_range[0],
    #                                                                 maximum_frequency=freq_range[1],
    #                                                                 method='fast')
    #            y = power
    #            x = 1/freq
    
        if color_by_pixel is False:
            color = 'k'
        else:
            rgb = c.cmap(c.norm(np.nanmedian(box,0)[i,j]))
            color = mpl.colors.rgb2hex(rgb)
    
        ax.plot(x, y, c=color)
    
        j += 1
        if j == colrange[1]:
            i += 1
            j  = colrange[0]
    
        if ylim is None:
            ax.set_ylim(np.percentile(y, 1), np.percentile(y, 99))
        else:
            ax.set_ylim(ylim[0], ylim[1])
    
        if xlim is None:
            ax.set_xlim(np.min(x)-0.1, np.max(x)+0.1)
        else:
            ax.set_xlim(xlim[0], xlim[1])
    
#        if data_type.lower() == 'amplitude':
#            ax.set_yscale('log')
#            ax.set_xscale('log')
#            ax.set_ylim(y.min(), y.max())
#            ax.set_xlim(np.min(x),
#                        np.max(x))
    
        ax.set_xticks([])
        ax.set_yticks([])
    
        fig.add_subplot(ax)
    plt.savefig('pixbypix_X%iby%i_Y%iby%i.png'%(rowrange[0],rowrange[1],colrange[0],colrange[1]),dpi=400,bbox_inches='tight')
    return fig
        
#%%
pixel_by_pixel(stampMedSub,[1,3],[1,3],ylim=[500,500])
#DUMP
#%%

#%% build in sliding detrender MATLAB CODE -- don't need?

%remove artifact
data_artifact_rmvd = zeros(size(raw_data));
current_section = [];

%step size
s = sr*param.ArtifactWidth;

for ii = 1:s:(length(raw_data)-s)
    if range(dtrn(ii:ii+s-1))<(3.0*param.ArtifactTh) && ii<(length(raw_data)-s*2)%2x up/down +1x for wiggle
        current_section = cat(2,current_section,dtrn(ii:ii+s-1));
    else
        if length(current_section)>(param.MinSecLength*sr)
            
            i = ii-length(current_section);
            data_artifact_rmvd(i:ii-1)= detrend(current_section);
            
            %either this or a loop for each sw_new
            current_section = [];
        else
            current_section = [];
            continue
        end
        continue
    end
    
end

chan_data.ProcessedData = data_artifact_rmvd;
chan_data.ArtifactTh = param.ArtifactTh;
%update channel data
data.Channels{chan} = chan_data;




l#%%
#DUMP
#%
#os.chdir('/Users/kookoo2025/Dropbox/2019-20/Blue Binaries/k2mosaic-master/')
#filesPath = '../../../../../../Volumes/rmathieu/NGC6791/'#/k2mosaic-master/'
#%
#exFile = 'kepler/kplr002297232-2010078095331_lpd-targ.fits'
#exFile = 'cadence-fits/k2mosaic-q04-ch80-cad13433.fits'
##%%
#idx = 1519 #0 through 1519
#for idx in [0,1,3,200,250,300,375,376]:
#    cadStr = str(cadList[idx])
#    fileStr = 'k2mosaic-q04-ch81-cad'+cadStr+'.fits'
#    filePath = filesPath+fileStr
##    infoprint(filePath,1)
#    [hdr,img] = getvar(fileStr,'IMAGE',0)
#    quick_imshow_savepngmv(img[704:943,10:239],'Cadence '+cadStr,'cadStamp',cadStr,'terrain',1)
#    

#%%

yMax,xMax = np.where(medianStamp == np.amax(medianStamp))
y =  yMax-9
x =  xMax+9
plt.plot(stampMedSub[:,x,y],label='Median Subtracted')
plt.plot(stampMeanSub[:,x,y],label='Mean Subtracted')
plt.ylim(-10,10)
plt.legend()
plt.show()




#%%