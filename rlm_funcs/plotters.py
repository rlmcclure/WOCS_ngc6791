#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 15:42:49 2020

@author: rlm
"""

import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import lightkurve as lk

#%%
def quick_imshow(twod_field,field_title='Quick Plot',cmapp='viridis',cmin=None,cmax=None,rtnfig=0):
    '''
    plots a two-dimensional field using imshow with option to return the figure object
    ---------
    params:
        twod_field = a two dimensional array to plot
        field_title = plot title
        cmapp = color map, defaults to viridis
        cmin = color range minimum, defaults to the field minimum
        cmax = color range maximum, defaults to the field maximum
        rtnfig = boolean to return the figure as an object, defaults to False
        
    returns:
        fig (optional) = object of the figure, for itteration over
        figure 
    '''
    if cmin==None:
        cmin=np.amin(twod_field)
        cmax=np.amax(twod_field)
    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.imshow(twod_field,origin='l',cmap=cmapp,clim=(cmin,cmax))
    fig.suptitle(field_title)
    sm = plt.cm.ScalarMappable(cmap=cmapp,norm=plt.Normalize(cmin,cmax))
    sm._A=[]
    cb = plt.colorbar(sm,fraction=0.046, pad=0.04)
    cb.ax.tick_params(labelsize=8)
    if rtnfig==1:
        return(fig)
    else:
        plt.show()
#%%
def quick_imshow_savepngmv(twod_field,field_title='Quick Plot',save_str='quickplot',n=0,cmapp='viridis',show=0,save=1,cmin=None,cmax=None,rtnfig=0):
    '''
    plots a two-dimensional field using imshow with option to return the figure object
    ---------
    params:
        twod_field = a two dimensional array to plot
        field_title = plot title
        save_str = string for the file save name
        n = for itteration over many files, defaults to 0
        cmapp = color map, defaults to viridis
        show = boolean to show the figure, defaults to False
        save = boolean to save the figure, defaults to True
        cmin = color range minimum, defaults to the field minimum
        cmax = color range maximum, defaults to the field maximum
        rtnfig = boolean to return the figure as an object, defaults to False
        
    returns:
        fig (optional) = object of the figure, for itteration over
        figure 
    '''
    if cmin==None:
        cmin=np.amin(twod_field)
        cmax=np.amax(twod_field)
    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.imshow(twod_field,origin='l',cmap=cmapp,interpolation='none',clim=(cmin,cmax))
    fig.suptitle(field_title)
    sm = plt.cm.ScalarMappable(cmap=cmapp,norm=plt.Normalize(cmin,cmax))
    sm._A=[]
    cb = plt.colorbar(sm,fraction=0.046, pad=0.04)
    cb.ax.tick_params(labelsize=8)
    if save==1:
        plt.savefig('%s%s.png'%(save_str,n),dpi=200,bbox_inches='tight')
    if rtnfig==1:
        return(fig)
    elif show ==1:
        plt.show()
    else:
        plt.close()
        
#%% pixel by pixel plot from eleanor.Visualization tess package
def pixel_by_pixel(box, colrange=None, rowrange=None, cmap='BrBG', mask=None,
                    xlim=None,ylim=None, color_by_pixel=True, 
                    save=1,savestr=None):
    
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
    
    ## PLOTS TARGET PIXEL FILE ##
    ax = plt.subplot(outer[0])
    
    c = ax.imshow(np.nanmedian(box[:,rowrange[0]:rowrange[1],colrange[0]:colrange[1]],0),origin='l',cmap=cmap)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.15)
    plt.colorbar(c, cax=cax, orientation='vertical')
    
    ## PLOTS PIXEL LIGHT CURVES ##
    for ind in range( int(nrows * ncols) ):
        ax = plt.Subplot(fig, inner[ind])
    
        flux = box[:,i,j]
        y = flux
        x = np.arange(len(y))
        
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
    
        ax.set_xticks([])
        ax.set_yticks([])
    
        fig.add_subplot(ax)
    
    if save == 1:
        plt.savefig('%s_pixbypix_X%iby%i_Y%iby%i.png'%(savestr,rowrange[0],rowrange[1],colrange[0],colrange[1]),dpi=400,bbox_inches='tight')
    return fig