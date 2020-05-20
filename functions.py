#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 12:13:38 2020

@author: kookoo2052

Helper Functions
"""
#%% Astrometry query
def run_ast(outfile):
    try_again = True
    submission_id = None
    
    while try_again:
        try:
            if not submission_id:
                wcs_header = ast.solve_from_image(outfile,
                                                  submission_id=submission_id, 
                                                  detect_threshold = 3)
                
            else:
                wcs_header = ast.monitor_submission(submission_id,
                                                    solve_timeout=420)
        except TimeoutError as e:
            submission_id = e.args[1]
        else:
            # got a result, so terminate
            try_again = False
    
    if wcs_header:
        # Code to execute when solve succeeds
        print('\nAstrometry Complete\n')
        return(wcs_header)
    else:
        # Code to execute when solve fails
        print('\nAstrometry Failed\n')
#%% pixel by pixel plot from eleanor.Visualization tess package
def pixel_by_pixel(box, colrange=None, rowrange=None, cmap='BrBG', mask=None,
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
    plt.savefig('pixbypix_X%iby%i_Y%iby%i.png'%(rowrange[0],rowrange[1],colrange[0],colrange[1]),dpi=400,bbox_inches='tight')
    return fig
#%%
def fistar2df(filesstr):
    '''
    Takes a string file location and outputs a pandas data frame with X&Y
    Index the data frame to get the X(Y) loc by df.loc[TARGETID]['X(Y)']
    _________________
    Params:(string) file path
    
    Returns:(pandas.DataFrame) with X, Y, & FWHM col with index from star ID in fistar
    _________________
    '''
    with open(filesstr) as f:    
        # first we allocate memory
        siz = 0
        for l in f:
            if l[0] != '#':
                siz+=1
        ID = np.empty(siz)
        X = np.empty(siz)
        Y = np.empty(siz)
        FWHM = np.empty(siz)
        
    with open(filesstr) as f:  
        ii = 0
        for l in f:
            if l[0] != '#':
                ID[ii] = l[:7]
                
                X[ii] = float(l[7:16])
                
                Y[ii] = float(l[16:25])
                
                FWHM[ii] = float(l[25:31])
                
                ii+=1
                    
    df = pd.DataFrame(data={'X':X,'Y':Y,'FWHM':FWHM},index=ID)
        
#            etc cols:  
#            NPix 
#            sigma  
#            delta  
#            Flux      
#            SN    
#            Bg      
#            Amp      
#            S      
#            D      
#            K    
#            CMax   

    return(df)
#%%
def infoprint(cubestr,verbose=None,varindx=None):
    '''
    a function to print out the header info for the stuff in these fits files
    --------
    params:
        cubestr = filestring for cube
        verbose = bool, defaults to true for printing the info for all headers
        
    returns:
        fileinfo: list of string names of variables in file 
    --------
    '''
    #open file
    with fits.open(cubestr) as hdul:
    
        #if there's a specified variable, print only that
        if varindx != None:
            fileinfo = hdul[varindx].header
            
            #if output on then print header info for each variable
            if verbose != 0:
                print(' ')
                print(fileinfo)
                print(' ')
                        
        #otherwise do all
        else:
            #get the string names of variables in HDUList
            fileinfo = [n[1] for n in hdul.info(output=False)]
            
            #step thorugh and print headers
            for ii in np.arange(np.size(hdul)):
                    hdr = hdul[ii].header
                    #if output on then print header info for each variable
                    if verbose != 0:
                        print(' ')
                        print(fileinfo[ii])
                        print(hdr)
                        print(' ')
                    
#    #close file
#    hdul.close()
    return(fileinfo)
    
#%%
def getvar(cubestr,varindx,verbose=None,neststr=None):
    '''
    a function to snag one variable's cube from a fits file
    --------
    params:
        cubestr = filestring for cube
        varindx = variable's number or string index
        verbose = bool, defaults to true for printing the info for all headers
        neststr = nested variable's string index (optional)
        
    returns:
        varstr: list of string names of variables in file 
        vardata: the variable file data, hdul[var_indx].data
    --------
    '''
    #open file
    with fits.open(cubestr) as hdul:
    
        #get the string names of variables in HDUList
        fileinfo = [n[1] for n in hdul.info(output=False)]
        
        #get variable data
        vardata = hdul[varindx].data
        hdr = hdul[varindx].header
        
        #if nested index is present
        if neststr != None:
    #        parent_vardata = np.copy(vardata)
            vardata = hdr[neststr]
            
            
        #if output on then print variable string
        if verbose != 0:
            print(' ')
            if type(varindx) != str:
                print(fileinfo[varindx])
            print(hdr)
            print(' ')
                
#    #close file
#    hdul.close()
    return(hdr,vardata)
#%%
def quick_imshow(twod_field,field_str,cmapp,cmin=None,cmax=None,rtnfig=0):
    if cmin==None:
        cmin=np.amin(twod_field)
        cmax=np.amax(twod_field)
    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.imshow(twod_field,origin='l',cmap=cmapp,clim=(cmin,cmax))
    fig.suptitle(field_str)
    sm = plt.cm.ScalarMappable(cmap=cmapp,norm=plt.Normalize(cmin,cmax))
    sm._A=[]
    cb = plt.colorbar(sm,fraction=0.046, pad=0.04)
    cb.ax.tick_params(labelsize=8)
    if rtnfig==1:
        return(fig)
    else:
        plt.show()
    #%%
def quick_imshow_savepngmv(twod_field,field_title,save_str,n,cmapp,show=0,save=1,cmin=None,cmax=None,rtnfig=0):
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

