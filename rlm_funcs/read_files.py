#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 16:01:10 2020

@author: rlm
"""
import os
import numpy as np
from astropy.io import fits
import pandas as pd

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
        if verbose != None:
            print(' ')
            if type(varindx) != str:
                print(fileinfo[varindx])
            print(hdr)
            print(' ')
                
    return(hdr,vardata)
    #%%=============================================================
#=====open different type of files, from M. Soares =============
#===============================================================

def readcolumn(var,col,infile,datformat='float',div=None,fill=False):
	fin=open(infile,mode='r')
	data=fin.readline().split(div)
	if (datformat=='float'):
		while((not len(data)==0) and (not data == [''])):
			if(not (data[0].startswith('#'))):
				try:
					var.append(float(data[col-1]))
				except ValueError:
					if(fill):
						var.append(np.nan)
					else:
						raise
				except IndexError:
					if(fill):
						var.append(np.nan)
					else:
						raise
			data=fin.readline().split(div)

	if (datformat=='str'):
		while((not len(data)==0) and (not data == [''])):
			if(not (data[0].startswith('#'))):
				try:
					var.append(data[col-1])
				except IndexError:
					if(fill):
						var.append(None)
					else:
						raise IndexError

			data=fin.readline().split(div)
	fin.close()
	return

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