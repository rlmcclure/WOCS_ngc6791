#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 12 16:33:52 2020

@author: rlmcclure

Procedure document for all NGC 6791 Blue Lurkers Search
"""
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
Q3, Q7, Q11, Q15 for channel 53 
Q2, Q6, Q10, Q14 for channel 1
Q1, Q5, Q9,  Q13, Q17 for channel 29 
(Q1 is only 34 days, Q0 is 9 days, the rest are 93 days, see Haas et al. 2010)

We start with just channel 81 for our initial round of analysis 
Would be nice to figure out how to extract just the super stamp from this
'''
#%% FITSH steps https://fitsh.net/wiki/fistar
'''
fistar command over the whole directory
    for i in *.fits; do fistar "$i" -o $(echo $i | sed s/fits/fistar/)  -s flux 
    --comment --model elliptic --flux-threshold 500 --algorithm uplink 
    --iterations symmetric=2,general=1 
    --format id,x,y,fwhm,npix,sigma,delta,flux,s/n,bg,amp,s,d,k,cmax
    --mag-flux 16.62,30; done
    
    --format id,x,y,bg,amp,s,d,k,flux,s/n,cmax,fwhm,npix,sigma,delta--mag-flux 16.62,30


if we have to automatically identify the stamp within a field,
    do a sliding ID by when np.sum(np.isnan(img[723:923,20:220]))==0

lowest FWHM of the fistar output files copy it as astroref.fistar for all 
    frames to get shifted to transform the files from the each quarter  
    using the grmatch and fitrans commands
    
    grmatch will match points between the astroref.fistar file and each of the 
    cadence files. It will spit out a file that ends in .itrans for each of the 
    cadences, which contains the astrometric transformation information
    
    fitrans step actually performs the transformation and spits out a file that 
    ends in .xtrns for each cadence. The .xtrns files are the rotated cadence 
    files that we will use to create a photometric reference frame and to 
    perform image subtraction.
'''
#%%
import os
#import matplotlib as mpl
#import matplotlib.pyplot as plt
#import matplotlib.gridspec as gridspec
import numpy as np
from astropy.io import fits
import scipy.signal as sig
from scipy import stats
#import lightkurve as lk
import pandas as pd

#from astroquery.astrometry_net import AstrometryNet
#ast = AstrometryNet()
#ast.api_key = 'dnzncgdyolxsukyy'

#from astroquery.mast import Catalogs
#from astroquery.vizier import Vizier
#v = Vizier()

from astropy import coordinates
from astropy import units as u
from astropy.table import Table

import kplr
client = kplr.API()

#from rlm_funcs.plotters import *
from rlm_funcs.read_files import *

#%% run functions
#runfile('/Users/kookoo2052/Dropbox/2019-20/WOCS/WOCS_ngc6791/functions.py', wdir='/Users/kookoo2052/Dropbox/2019-20/WOCS/WOCS_ngc6791')
#%%
toggle_runcmd = 0
osfailstr = '\nheyhi so that command failed... plz investigate.\n'
#%% set quarter and channel
qu = 4
ch = 81
#%% until astrometry is inline, params here
raC = 290.5004083333
decC = 38.339662222222
#%% params
fluxthld = 597
#%%
filesPath = '/Users/kookoo2052/HOLD/mosaic-ch81/q04c81'#q4
#filesPath = '/Volumes/rmathieu/NGC6791/q04c81'
filesPathfits = '/Volumes/GoogleDrive/Shared drives/WOCS Mathieu/NGC 6791/q04c81'
filesPathfits = '/Volumes/GoogleDrive/Shared\ drives/WOCS\ Mathieu/NGC\ 6791/q04c81'
#%%
#filesPath = filesPathr

filesPathfits = '/Volumes/GoogleDrive/Shared drives/WOCS Mathieu/NGC 6791/q04c81'
filesList = [f for f in os.listdir(filesPathfits) if (f.startswith("k2mosaic") and f.endswith(".fits") and os.path.isfile(os.path.join(filesPathfits, f)) and not f.endswith("ref.fits") and not f.endswith("-xtrns.fits") and not f.endswith("-xtrns-m.fits") and not f.endswith("-conv.fits") and not f.endswith("-subtr.fits"))]
filesList.sort()

#%% double check for a flag of bad cadences 
#flaged vs cadence by cadence drift anamolies 
#%% loop for fistar
toggle_runcmd=0
#filesPath=filesPathl
filesPathfits = '/Volumes/GoogleDrive/Shared\ drives/WOCS\ Mathieu/NGC\ 6791/q04c81'

fistarList = []
for ii, f in enumerate(filesList[:]):
    infile = filesPathfits+'/'+f
    fistarfile = f[:-5]+'.fistar'
    fistarList.append(fistarfile)
    outfile = filesPath+'/'+fistarfile
    
    if toggle_runcmd == 1:
        cmdline = 'fistar ' + infile + ' -o ' + outfile + ' -s flux --comment --model elliptic --flux-threshold ' + str(fluxthld) +' --algorithm uplink --iterations symmetric=2,general=1 --format id,x,y,fwhm,npix,sigma,delta,flux,s/n,bg,amp,s,d,k,cmax --mag-flux 16.62,30'
        out=os.system(cmdline)
        if out != 0:
            print(osfailstr)
            toggle_runcmd = 0
    

#    for i in *.fits; do fistar "$i" -o $(echo $i | sed s/fits/fistar/)  -s flux 
#    --comment --model elliptic --flux-threshold 500 --algorithm uplink 
#    --iterations symmetric=2,general=1 
#    --format id,x,y,fwhm,npix,sigma,delta,flux,s/n,bg,amp,s,d,k,cmax
#    --mag-flux 16.62,30; done
##%%
#fistarList = [f for f in os.listdir(filesPath) if (f.endswith(".fistar") and os.path.isfile(os.path.join(filesPath, f)) and not f.endswith("-xtrns.fistar") and not f.startswith("astroref"))]
#fistarList.sort()



#%% find minimum FWHM fistar file for astroref.
fidFistarSuccess = []
fidxFistarSuccess = []
fidCut = []
fidxCut = []
nCut = 0

#filesPath = filesPathr

for x in np.arange(len(fistarList)):
    df = fistar2df(filesPath+'/'+fistarList[x])
    if df.empty != 1:
        fidxFistarSuccess.append(x)
        fidFistarSuccess.append(fistarList[x][:-7])
    else:
        fidxCut.append(x)
        fidCut.append(fistarList[x][:-7])
        nCut += 1
        
#%% pull all bad cadence quality flags, and then all quality flags
'''
KEPLER_QUALITY_FLAGS = {
    "1": "Attitude tweak",
    "2": "Safe mode",
    "4": "Coarse point",
    "8": "Earth point",
    "16": "Zero crossing",
    "32": "Desaturation event",
    "64": "Argabrightening",
    "128": "Cosmic ray",
    "256": "Manual exclude",
    "1024": "Sudden sensitivity dropout",
    "2048": "Impulsive outlier",
    "4096": "Argabrightening",
    "8192": "Cosmic ray",
    "16384": "Detector anomaly",
    "32768": "No fine point",
    "65536": "No data",
    "131072": "Rolling band",
    "262144": "Rolling band",
    "524288": "Possible thruster firing",
    "1048576": "Thruster firing"
}
'''
filesPathfits = '/Volumes/GoogleDrive/Shared drives/WOCS Mathieu/NGC 6791/q04c81'

q = []
for ii, f in enumerate(fidCut):
    filen = filesPathfits+'/'+f+'.fits'
    hd,val = getvar(filen,'IMAGE',neststr='QUALITY')
    q.append(val)
    
qKEEP = []
for ii, f in enumerate(fidFistarSuccess):
    filen = filesPathfits+'/'+f+'.fits'
    hd,val = getvar(filen,'IMAGE',neststr='QUALITY')
    qKEEP.append(val)

print('\n%2f percent had bad flags in QUARTER %2i CHANNEL %2i.\n'%(100*len(np.where(np.asarray(qKEEP))[0])/len(fidFistarSuccess),qu,ch)) 

fid = np.asarray(fidFistarSuccess)[np.where(np.asarray(qKEEP)==0)[0]]
fidx = np.asarray(fidxFistarSuccess)[np.where(np.asarray(qKEEP)==0)[0]]

#%%
siz=len(fid) #fileIDs of good fields
fwhmMed = np.empty(siz)
fwhmMean = np.empty(siz)
fwhmMin = np.empty(siz)
fwhmMax = np.empty(siz)
fwhmStd = np.empty(siz)
#%WHAT IS GOING ON RIGHT HERE??????????????????????????????????????????????????
for ii, x in enumerate(fidx):
    df = fistar2df(filesPath+'/'+fistarList[x])
    fwhmMed[ii] = np.median(df['FWHM'])
    fwhmMean[ii] = np.mean(df['FWHM'])
    fwhmMin[ii] = np.amin(df['FWHM'])
    fwhmMax[ii] = np.amax(df['FWHM'])
    fwhmStd[ii] = np.std(df['FWHM'])
    
#%
dfFWHM = pd.DataFrame(data={'median':fwhmMed,'mean':fwhmMean,'std':fwhmStd,'min':fwhmMin,'max':fwhmMax},index=fid)
del(fwhmMed,fwhmMean,fwhmStd,fwhmMin,fwhmMax)
#%
print('\nMinimum FWHM:\n',dfFWHM.loc[dfFWHM['min'].idxmin])
print('\nMinimum Mean FWHM:\n',dfFWHM.loc[dfFWHM['mean'].idxmin])
#% Make copy of minimum FWHM file as astroref file
astroref = filesPath+'/'+'astroref.fistar'
astrorefits = filesPath+'/'+'astroref.fits'

#% toggle to run commands
if toggle_runcmd == 1:
    cmdline = 'cp ' + filesPath+'/'+dfFWHM.loc[dfFWHM['mean'].idxmin].name+'.fistar' + ' ' + astroref
    cmdline = 'cp ' + filesPath+'/'+dfFWHM.loc[dfFWHM['mean'].idxmin].name+'.fits' + ' ' + astrorefits
    out = os.system(cmdline)
    if out != 0:
        print(osfailstr)
        toggle_runcmd = 0

#%% Plot the median FWHM as a function of frame number
#dfFWHM.plot(y='mean')


#%% grmatch
'''
grmatch will match points between the astroref.fistar file and each of the 
    cadence files. It will spit out a file that ends in .itrans for each of the 
    cadences, which contains the astrometric transformation information
'''
trid = []
resd = np.empty(siz)
resdW = np.empty(siz)
utry = np.empty(siz)
rato = np.empty(siz)

#toggle_runcmd = 1
#filesPath = filesPathl
#loop for grmatch
for ii, f in enumerate(fid):
    infile = filesPath+'/'+f+'.fistar'
    transfile = infile+'.itrans'
    trid.append(transfile)
        
    #% toggle to run commands
    if toggle_runcmd == 1:
        cmdline = 'grmatch --match-points -r ' + astroref + ' --col-ref 2,3 --col-ref-ordering +8 -i ' + infile + ' --col-inp 2,3 --col-inp-ordering +8 --weight reference,column=12 --triangulation maxinp=1000,maxref=1000,auto,unitarity=0.01 --order 2 --max-distance 1 --comment --output-transformation ' + transfile + ' --output /dev/null'
        out = os.system(cmdline)
        if out != 0:
            print(osfailstr)
            toggle_runcmd = 0

    with open(transfile) as ff:  
        for l in ff:
            if l[2:10] == 'Residual':
                if l[-9:] == '(native)\n':
                    if 'nan' in l[11:-9]:
                        resd[ii] = float('Nan')
                    else:
                        resd[ii] = float(l[11:-9])
                else:
                    if 'nan' in l[11:-11]:
                        resdW[ii] = float('Nan')
                    else:
                        resdW[ii] = float(l[11:-11])
            elif l[2:11] == 'Unitarity':
                utry[ii] = float(l[13:-1])
            elif l[2:7] == 'Ratio':
                rato[ii] = float(l[8:-10])
                
dfGRM = pd.DataFrame(data={'Residuals':resd,'Weighted Residuals':resdW,'Unitarity':utry,'Ratio':rato},index=trid)
del(resd,resdW,utry,trid)
#%% fitrans do this in chunks https://fitsh.net/wiki/man/fistar
'''
fitrans step actually performs the transformation and spits out a file that 
    ends in -xtrns for each cadence. The -xtrns files are the rotated cadence 
    files that we will use to create a photometric reference frame and to 
    perform image subtraction.
'''
#loop for fitrans
for ii, f in enumerate(fid):
    infile = filesPath+'/'+f+'.fits'
    transfile = infile[:-5]+'.itrans'
    xtrnsfile = infile[:-5]+'-xtrns.fits'
    
    #% toggle to run commands
    if toggle_runcmd == 1:
        cmdline = 'fitrans ' + infile + ' -k --input-transformation ' + transfile + ' --reverse -o ' + xtrnsfile 
        out = os.system(cmdline)
        if out != 0:
            print(osfailstr)
            toggle_runcmd = 0
            
#%% mask the frames fiign do this https://fitsh.net/wiki/man/fiign
'''
Mask the files from the transformation onward 
    $ fiign -a -g --mask-value -999 -i files.xtrns -o files.xtrns

all .xtrns files, the PhotRef, and then the image subtracted frames will also be masked
Keep copies of the unmasked files too
'''

#toggle_runcmd=1
#loop for fitrans
for ii, f in enumerate(fid):
    xtrnsfilein = filesPath+'/'+f+'-xtrns.fits'
    xtrnsfileout = filesPath+'/'+f+'-xtrns-m.fits'
    
    #% toggle to run commands
    if toggle_runcmd == 1:
        cmdline = 'fiign -a -g --mask-value -999 -i ' + xtrnsfilein + ' -o ' + xtrnsfileout
        out = os.system(cmdline)
        if out != 0:
            print(osfailstr)
            toggle_runcmd = 0
            
toggle_runcmd=0
#%% photometric reference frame USE FICOMBINE
'''
using the sharpest 20 (40) FWHM as fcombine fails because there are too many files
and we want a sharp reference frame

easiest way may be to copy the 20 into another file and then preform the photoref in that folder

'''
photorefs = dfFWHM.nsmallest(20,'mean').index
prxtrnsm = list(photorefs)
for ii, pr in enumerate(photorefs):
    prxtrnsm[ii] = filesPath+'/'+str(pr+'-xtrns-m.fits')
#%
#toggle_runcmd = 1 
#infile = filesPath+'/'+'*xtrns.fits'
outfile = filesPath+'/'+'photoref.fits'
#% toggle to run commands
if toggle_runcmd == 1:
    cmdline = 'ficombine ' + str(prxtrnsm).replace('\'','').replace(',','')[1:-1] + ' --mode mean -o ' + outfile 
    out = os.system(cmdline)
    if out != 0:
        print(osfailstr)
        toggle_runcmd = 0
        
toggle_runcmd = 0
            
#%% run grtrans
'''
astrometry for the stamp as csv using the photoref.fits above
    get the RA & Dec from teh WCS header I assume or some part of the astrometry steps
take the ID RA Dec and Gmag and make into a space separated column .dat file 'gaia_ra_dec_gmag.dat'



grtrans --comment --input /gaia_ra_dec_gmag.dat --wcs arc,ra=raC,dec=decC --col-radec 2,3 --col-out 5,6 --output /xi_eta_grtrans_out.dat
'''
#toggle_runcmd = 1
#run for fistar
radec = filesPath+'/'+'gaia_ra_dec_gmag.dat'
xieta = filesPath+'/'+'xi_eta_grtrans_out.dat'
    
#% toggle to run commands
if toggle_runcmd == 1:
    cmdline = 'grtrans '+'--comment --input '+radec+' --wcs arc,ra='+str(raC)+',dec='+str(decC)+' --col-radec 2,3 --col-out 5,6 --output  '+xieta
    out = os.system(cmdline)
    if out != 0:
        print(osfailstr)
        toggle_runcmd = 0
            
toggle_runcmd = 0
#%% create .reg file form melinda_regslct.py
'''
the region file separates the astrometry output into X evenly spaced regions and then picks the brightest non saturated star in each to create an alignment grid

use an astroref.fistar, the columns in this astroref.fistar are: x coordinates, y coordinates, and flux.
grid the sources in the fistar into a 30x30 grid then find the brightest non-saturated star in each grid
the coordinates and flux for said star will be output into a .reg file
'''

astroref = filesPath+'/'+'astroref.fistar'
regfile = filesPath+'/'+'photoref.reg'
#toggle_runcmd = 1
#% toggle to run commands
if toggle_runcmd == 1:

    #column indexes in the fistar file
    colx = 2; coly=3; colflux = 5

    #read in from fistar
    # x locs arr
    xcoor = []; readcolumn(xcoor,colx,astroref)
    xcoor=np.array(xcoor)

    # y locs arr
    ycoor = []; readcolumn(ycoor,coly,astroref)
    ycoor = np.array(ycoor)

    # flux arr
    flux = []; readcolumn(flux,colflux,astroref)
    flux = np.array(flux)

    #the number of grid points in x and y direction
    grid_x = 30; grid_y = 30
    mx = np.zeros(grid_x*grid_y)-1
    my = np.zeros(grid_x*grid_y)-1
    ma = np.zeros(grid_x*grid_y)

    #the pixel numbers in x and y direction, might need to change
    hd,pix_x=getvar(filesPath+'/'+fid[0]+'.fits','IMAGE',neststr='NAXIS1')
    hd,pix_y=getvar(filesPath+'/'+fid[0]+'.fits','IMAGE',neststr='NAXIS2')
    
    saturation = 20000 #64000. #saturation limit, not sure if this is correct

    #taking the image coordinate value stored as xcoor/ycoor multiply this by grid number and divide this value by the total number of pixels in the x/y directions.

    # determine which grid a paticular star belongs to.
    bx = (xcoor*grid_x/pix_x).astype(int)
    by = (ycoor*grid_y/pix_y).astype(int)

    #give each grid the coordinates of the brightest, non-saturated star in the grid
    for i in range(len(xcoor)):
        j = by[i]*bx[i]+bx[i]
        if ma[j]<flux[i] and flux[i]<saturation:
            mx[j] = xcoor[i]
            my[j] = ycoor[i]
            ma[j] = flux[i]
    #           print 'mx[j] = ', mx[j], 'my[j]=', my[j]
    with open(regfile,mode='w') as fout:
        for i in range(int(grid_x*grid_y)):
            fout.write("%8.0f %8.0f %8.0f\n" % (mx[i],my[i],ma[i]))

toggle_runcmd = 0
#created the file region file with cols: x image coordinate, y image coordinate, flux


#%% Convolve the frames before we do the subtraction too. ficonv [options] [-r <reference>] [-i <input>] [-o <output>]

#%% loop for ficonv
'''
Takes a *.itrans, *xtrns.fits, ??region file, & stacked photref frame reffile in the script.
Results in the output of *.kernel and -sub.fits files. 

ficonv -i photref.fits -r [loop file]-xtrns.fits -it PhotRef.reg -k "b/0;i/0;d=[1]/0" -ok [loop file].kernel -oc [loop file]-conv.fits

ficonv [options] [-r <reference>] [-i <input>] [-o <output>]

The purpose of this program is twofold. First, using a reference and input FITS
 image, the program tries to figure out the best-fit convolution kernel which
 transforms the reference image to the input image. The convolution kernel
 solution is saved as a separate, human-readable file. Second, using an existing
 such kernel solution file, `ficonv` convolves the reference image and optionally
 co-adds an additional image to the result. The program figures out the desired
 mode (i.e. whether to fit a convolution kernel or use an existing kernel solution
 to convolve an image) from the presence or lack of the command line options.

-r, --input-reference <image file>
    Name of the reference FITS image file.
-i, --input <image file>
    Name of the input FITS image file (required only for kernel fitting).
?-k, --kernel <kernel set>
    List of kernel bases used for fitting convolution kernel.See also Kernel specifications below for the format of this <kernel set> argument.
?--input-kernel <file>
    Name of the file containing kernel bases. The kernel bases in this file should have no associated coefficients if convolution fitting is done, otherwise the kernel basis file must contain the convolution coefficients
--output-kernel <file>
    Name of the file where the coefficients for the kernel bases are stored after convolution kernel fitting
-o, -oc, --output, --output-convolved <image file>
    Name of the output file which is the reference image convolved with the kernel solution (which can either be a previously fitted and now read from a file or the result of the current fit)
--output-subtracted <image file>
    The difference between the input image and the convolved reference image.
'''
toggle_runcmd = 1

fluxthld = 5

d = 2#1

#fistarListShift = []
for ii, f in enumerate(fid):
    infile = filesPath+'/'+f+'-xtrns-m.fits'
    reffile = filesPath+'/'+'photoref.fits'
    regfile = filesPath+'/'+'photoref.reg'
    outsub = filesPath+'/'+f+'-subtr.fits'
    outfile = filesPath+'/'+f+'-conv.fits'
    outkern = filesPath+'/'+f+'.kernel'
    
#%   fistarListShift.append(fistarfile)
    
    if toggle_runcmd == 1:
#        cmdline = 'ficonv'+' -i '+reffile+' -r '+infile+'  -k "b/0;i/0;d='+str(d)+'/0"'+' --output-kernel-list '+outkern+' --output-subtracted  '+outsub+' -o '+outfile
        cmdline = 'ficonv'+' -i '+reffile+' -r '+infile+' --input-stamp '+regfile+'  -k "b/0;i/0;d='+str(d)+'/0"'+' --output-kernel-list '+outkern+' --output-subtracted  '+outsub+' -o '+outfile
        
#        cmdline = 'ficonv' + ' -output-kernel ' + outkern + ' -output-subtracted ' + outsub + ' -r ' + reffile + ' -i ' + infile + ' -o ' + outfile 
        
        out = os.system(cmdline)
#        print(cmdline)
        if out != 0:
            print(cmdline)
            print(osfailstr)
            toggle_runcmd = 0
            
            
#toggle_runcmd = 0

#%% fistar with updated parameters just on the astroref file to make the most precise xi_eta transformation
'''
fistar astroref.fits -o astroref.fistar -s flux --comment --model elliptic --flux-threshold 597 --algorithm uplink --iterations symmetric=2,general=3 --format id,x,y,bg,amp,s,d,k,flux,s/n,fwhm,ellip,pa --mag-flux 16.62,30 --sort flux
'''
#toggle_runcmd = 1
#run for fistar
fitsfile = filesPath+'/'+'astroref.fits'
fistarfile = filesPath+'/'+'astroref.fistar'
    
#% toggle to run commands
if toggle_runcmd == 1:
    cmdline = 'fistar ' + fitsfile + ' -o ' + fistarfile + ' -s flux --comment --model elliptic --flux-threshold '+str(fluxthld)+' --algorithm uplink --iterations symmetric=2,general=3 --format id,x,y,bg,amp,s,d,k,flux,s/n,fwhm,ellip,pa --mag-flux 16.62,30 --sort flux'
    out = os.system(cmdline)
    if out != 0:
        print(osfailstr)
        toggle_runcmd = 0
            
#toggle_runcmd = 0

#%% grmatch with the final optimized frame
'''
grmatch --input-reference xi_eta_grtrans_out.dat --input
 astroref.fistar --col-ref 5,6 --col-ref-ordering -4
 --col-inp 2,3 --col-inp-ordering 9 --order 1 --max-distance 2 --triangulation
 unitarity=0.01,maxref=119,maxinp=119,full,mixed --fit iterations=2
 --weight reference,column=4,magnitude,power=2 --comment --output-transformation
 reference_gaia_matched.trans --output-matched
 reference_gaia_matched.txt
'''
#toggle_runcmd = 1
#run for grmatch
xieta = filesPath+'/'+'xi_eta_grtrans_out.dat'
fistar = filesPath+'/'+'astroref.fistar'
outtrns = filesPath+'/'+'reference_gaia_matched.trns'
outtxt = filesPath+'/'+'reference_gaia_matched.txt'

#% toggle to run commands
if toggle_runcmd == 1:
    cmdline = 'grmatch '+'--input-reference '+xieta+' --input '+fistar+' --col-ref 5,6 --col-ref-ordering -4 --col-inp 2,3 --col-inp-ordering 9 --order 1 --max-distance 2 --triangulation unitarity=0.01,maxref=119,maxinp=119,full,mixed --fit iterations=2 --weight reference,column=4,magnitude,power=2 --comment --output-transformation '+outtrns+' --output-matched '+ outtxt
    out = os.system(cmdline)
    if out != 0:
        print(osfailstr)
        toggle_runcmd = 0
            
#toggle_runcmd = 0


#%% grtrans with the final optimized frame
'''
grtrans --input xi_eta_grtrans_out.dat --col-xy 5,6 --input-transformation reference_gaia_matched.trans --col-out 7,8 --output gaia_full_transformed_to_astrometricref.txt
'''
#toggle_runcmd = 1
#run for grtrans
xieta = filesPath+'/'+'xi_eta_grtrans_out.dat'
intrns = filesPath+'/'+'reference_gaia_matched.trns'
outtxt = filesPath+'/'+'gaia_full_transformed_to_astrometricref.txt'

#% toggle to run commands
if toggle_runcmd == 1:
    cmdline = 'grtrans '+'--input '+xieta+' --col-xy 5,6 --input-transformation '+intrns+' --col-out 7,8 --output '+ outtxt
    out = os.system(cmdline)
    if out != 0:
        print(osfailstr)
        toggle_runcmd = 0
            
#toggle_runcmd = 0

#%% grtrans for the gaia to kepler transformation
'''
grtrans --input reference_gaia_matched.txt --col-xy 5,6 --input-transformation reference_gaia_matched.trans --col-out 20,21 --output gaiaxrK2_transformed_to_astrometricref.txt
'''
toggle_runcmd = 1
#run for grtrans
intrns = filesPath+'/'+'reference_gaia_matched.trns'
intxt = filesPath+'/'+'reference_gaia_matched.txt'
outtxt = filesPath+'/'+'gaiaxrK2_transformed_to_astrometricref.txt'

#% toggle to run commands
if toggle_runcmd == 1:
    cmdline = 'grtrans '+'--input '+intxt+' --col-xy 5,6 --input-transformation '+intrns+' --col-out 20,21 --output '+ outtxt
    out = os.system(cmdline)
    if out != 0:
        print(osfailstr)
        toggle_runcmd = 0
            
toggle_runcmd = 0



#%% loop for fiphot
'''
performs aperture photometry on normal or subtracted/convolved images

First the script calls the fiphot program in the subtracted frame 
 photometry mode,then it transforms one sets of x,y coordinates back into the 
 original frame coordinate using the .itrans file, and finally it prints
 out other relevant information that the light curve requires. 

Doing normal aperture photometry on the master image (i.e. the photometry reference), requires: 
- fiphot program used for aperture photometry
- list of x,y coordinates for the stars we want to do photometry on 
- fits image (photoref, with regions properly masked)
- list of aperture to do the photometry on **start with a big list, making plots, and then choose the good few ones. 
- method to convert magnitude to flux, I had a guess for Kepler, which is not too bad 

fiphot [options] [<input>] [-o|--output <output>]

Melinda calls

gain %.3f
zeropoint %.2f
apertures has the quotes around the string value somehow like '%s'

gain,zeropoint,apertures,serial set in file

fiphot --input refframe --input-list xylist --col-id colid 
 --col-xy colx,coly --gain gain --mag-flux zeropoint,30
 --apertures 'apertures' --sky-fit 'median,sigma=5,iterations=2'
 --disjoint-radius 2 --serial serial --format 'IXY-----,sMm'
 --nan-string 'NaN' --aperture-mask-ignore 'saturated' --comment '--comment'
 --single-background 3 -op rawphot -k
 
 We are doing fiphot multiple times
 
first:
fiphot --input /Users/msoaresfurtado/Dropbox/Research_NGC6791/Notebooks/photref.fits --input-list /Users/msoaresfurtado/Dropbox/Research_NGC6791/Notebooks/gaiaxrK2_transformed_to_astrometricref.txt --col-id 1 --col-xy 20,21 --col-mag 4 --apertures 2.5:5.0:6.0,2.6:5.0:6.0,2.7:5.0:6.0,2.8:5.0:6.0,2.9:5.0:6.0,3.0:5.0:6.0,3.1:5.0:6.0,3.2:5.0:6.0,3.3:5.0:6.0,3.4:5.0:6.0,3.5:5.0:6.0,3.55:5.0:6.0,3.6:5.0:6.0,3.65:5.0:6.0,3.7:5.0:6.0,3.75:5.0:6.0,3.8:5.0:6.0,3.85:5.0:6.0,3.9:5.0:6.0,3.95:5.0:6.0 --sky-fit mode,iterations=4,sigma=3 --disjoint-radius 2 --format IXY-----,sMm --comment --output-raw-photometry output/fiphot_output_raw_photometry_full.txt --output output/fiphot_output_full.txt
'''
#%%fiphot 1
'''
 fiphot --input photoref.fits --input-list gaia_full_transformed_to_astrometricref.txt --col-id 1 --col-xy 7,8 --col-mag 4 --apertures 2.5:5.0:6.0,2.6:5.0:6.0,2.7:5.0:6.0,2.8:5.0:6.0,2.9:5.0:6.0,3.0:5.0:6.0,3.1:5.0:6.0,3.2:5.0:6.0,3.3:5.0:6.0,3.4:5.0:6.0,3.5:5.0:6.0,3.55:5.0:6.0,3.6:5.0:6.0,3.65:5.0:6.0,3.7:5.0:6.0,3.75:5.0:6.0,3.8:5.0:6.0,3.85:5.0:6.0,3.9:5.0:6.0,3.95:5.0:6.0 --sky-fit mode,iterations=4,sigma=3 --disjoint-radius 2 --format IXY-----,sMm --comment --output-raw-photometry output/fiphot_output_raw_photometry_full.txt --output output/fiphot_output_full.txt
'''
#toggle_runcmd = 1

#files
reffile = filesPath+'/'+'photoref.fits'
photout = filesPath+'/'+'fiphot_output_raw_photometry_full.txt'
outfile = filesPath+'/'+'fiphot_output_full.txt'
#stuff
xylist = filesPath+'/'+'gaia_full_transformed_to_astrometricref.txt'
colx = 7; coly = 8; colid = 1; colmag = 4;
apertures = "2.5:5.0:6.0,2.6:5.0:6.0,2.7:5.0:6.0,2.8:5.0:6.0,2.9:5.0:6.0,3.0:5.0:6.0,3.1:5.0:6.0,3.2:5.0:6.0,3.3:5.0:6.0,3.4:5.0:6.0,3.5:5.0:6.0,3.55:5.0:6.0,3.6:5.0:6.0,3.65:5.0:6.0,3.7:5.0:6.0,3.75:5.0:6.0,3.8:5.0:6.0,3.85:5.0:6.0,3.9:5.0:6.0,3.95:5.0:6.0"

if toggle_runcmd == 1:
    cmdline = 'fiphot' + ' --input '+reffile+' --input-list '+xylist+' --col-id '+str(colid)+' --col-xy '+str(colx)+','+str(coly)+' --col-mag '+str(colmag)+' --apertures '+apertures+' --sky-fit mode,sigma=3,iterations=4'+' --disjoint-radius 2 '+' --format IXY-----,sMm '+' --comment '+'--output-raw-photometry '+photout+' --output '+outfile
    
    out = os.system(cmdline)
    if out != 0:
        print(cmdline)
        print(osfailstr)
        toggle_runcmd = 0
            
            
#toggle_runcmd = 0

#%% 
'''
So now, fitting the appertrues to a gaussian function, this gives you the fitting  functions
Now going to do photometry (next two fiphots) with the appropriate Gaia magnitudes and the fitting functions on all the frames
'''
#%% 3:20pm started
'''
second:
fiphot --input-subtracted cad15100-subtr.fits --input-raw-photometry fiphot_output_raw_photometry_full_processed.txt 
--sky-fit median,iterations=4,sigma=3 --format IXY-----,sMm --gain 105.99 --input-kernel k2mosaic-q04-ch81-cad15100.kernel 
--comment --output subtraction_photometry_output/photometry_output_15100.out
'''
#toggle_runcmd = 1


for ii, f in enumerate(fid):
    
    #files
    insub = filesPath+'/'+f+'-subtr.fits'
    inkern = filesPath+'/'+f+'.kernel'
    photout = filesPath+'/'+f+'-subphot.out'
    inphot = filesPath+'/'+'fiphot_output_raw_photometry_full.txt'
    
    #stuff
    iterations = 4
    sigma = 3
    gain = 105.99
    

    if toggle_runcmd == 1:
        cmdline = 'fiphot' + ' --input-subtracted '+insub+' --input-raw-photometry '+inphot+' --sky-fit median,iterations='+str(iterations)+',sigma='+str(sigma)+' --format IXY-----,sMm'+' --gain '+str(gain)+' --input-kernel '+inkern+' --comment'+' --output '+photout
        
        out = os.system(cmdline)    
        if out != 0:
            print(cmdline)
            print(osfailstr)
            toggle_runcmd = 0
            
            
#toggle_runcmd = 0
#%%
'''
why? the only difference here is using the gaia to kepler instead of the gaia transformed???

third:
fiphot --input photref.fits --input-list /gaiaxrK2_transformed_to_astrometricref.txt 
--col-id 1 --col-xy 20,21 --col-mag 4 --apertures 2.5:5.0:6.0,2.6:5.0:6.0,2.7:5.0:6.0,2.8:5.0:6.0,2.9:5.0:6.0,3.0:5.0:6.0,3.1:5.0:6.0,3.2:5.0:6.0,3.3:5.0:6.0,3.4:5.0:6.0,3.5:5.0:6.0,3.55:5.0:6.0,3.6:5.0:6.0,3.65:5.0:6.0,3.7:5.0:6.0,3.75:5.0:6.0,3.8:5.0:6.0,3.85:5.0:6.0,3.9:5.0:6.0,3.95:5.0:6.0 
--sky-fit mode,iterations=4,sigma=3 --disjoint-radius 2 --format IXY-----,sMm --comment --output-raw-photometry output/fiphot_output_raw_photometry_full.txt --output output/fiphot_output_full.txt
'''

#toggle_runcmd = 1


#files
reffile = filesPath+'/'+'photoref.fits'
trnsin = filesPath+'/'+'gaiaxrK2_transformed_to_astrometricref.txt'
photout = filesPath+'/'+'fiphot_output_raw_photometry_full.txt'
outfile = filesPath+'/'+'fiphot_output_full.txt'
#stuff
xylist = filesPath+'/'+'gaia_full_transformed_to_astrometricref.txt'
colx = 20; coly = 21; colid = 1; colmag = 4;
#apertures = "2.5:5.0:6.0,2.6:5.0:6.0,2.7:5.0:6.0,2.8:5.0:6.0,2.9:5.0:6.0,3.0:5.0:6.0,3.1:5.0:6.0,3.2:5.0:6.0,3.3:5.0:6.0,3.4:5.0:6.0,3.5:5.0:6.0,3.55:5.0:6.0,3.6:5.0:6.0,3.65:5.0:6.0,3.7:5.0:6.0,3.75:5.0:6.0,3.8:5.0:6.0,3.85:5.0:6.0,3.9:5.0:6.0,3.95:5.0:6.0"

if toggle_runcmd == 1:
    cmdline = 'fiphot' + ' --input '+reffile+' --input-list '+xylist+' --col-id '+str(colid)+' --col-xy '+str(colx)+','+str(coly)+' --col-mag '+str(colmag)+' --apertures '+apertures+' --sky-fit mode,sigma=3,iterations=4'+' --disjoint-radius 2 '+' --format IXY-----,sMm '+' --comment '+'--output-raw-photometry '+photout+' --output '+outfile
    
    out = os.system(cmdline)
    if out != 0:
        print(cmdline)
        print(osfailstr)
        toggle_runcmd = 0
            
            
#toggle_runcmd = 0



#%%



#%%
#%% loop for fiarith subtraction -- obsolete
'''
fiarith "'fiign-xtrns.fits'-'fiign-subtract.fits'" -o diff.fits
'''
#toggle_runcmd = 1

fluxthld = 5

fistarListShift = []
for ii, f in enumerate(fid):
    infile = filesPath+'/'+f+'-subtr.fits'
    ficonvfile = f+'-xtrns.fits'
    fistarListShift.append(fistarfile)
    outfile = filesPath+'/'+fistarfile
    
    
    while toggle_runcmd == 1:
        cmdline = 'fiarith ' + infile + ' -o ' + outfile + ' -s flux --comment --model elliptic --flux-threshold ' + str(fluxthld) +' --algorithm uplink --iterations symmetric=2,general=1 --format id,x,y,fwhm,npix,sigma,delta,flux,s/n,bg,amp,s,d,k,cmax --mag-flux 16.62,30'
        out = os.system(cmdline)
        if out != 0:
            os.error(osfailstr)
            toggle_runcmd = 0
            
            
#toggle_runcmd = 0

#%%
'''
Below here is the fun checking stuff things
-
--
---
----
-----
------
-------
--------
-------
------
-----
----
---
--
-
-
--
---
----
-----
------
-------
--------
-------
------
-----
----
---
--
-
-
--
---
----
-----
------
-------
--------
-------
------
-----
----
---
--
-
-
--
---
----
-----
------
-------
--------
-------
------
-----
----
---
--
-
'''

##%%Extract the aligned stamp to look at it
#
#stampBox = np.empty([len(fid),206,206])#,200,200])
#
##%Stamp loc [723:923,20:220]
#for ii, f in enumerate(fid):
#    
#    #Create the full path string
#    fileStr = f+'-xtrns.fits'
#    filePath = filesPath+'/'+fileStr
#    
#    #Pull imag for each file
#    [hdr,img] = getvar(filePath,'PRIMARY')
#    stampBox[ii]=img[720:926,17:223]#[723:923,20:220]
#
##%% extract the masked aligned stamp
#stampBoxMasked = np.empty([len(fid),206,206])#,200,200])
#
##%Stamp loc [723:923,20:220]
#for ii, f in enumerate(fid):
#    
#    #Create the full path string
#    fileStr = f+'-xtrns-m.fits'
#    filePath = filesPath+'/'+fileStr
#    
#    #Pull imag for each file
#    [hdr,img] = getvar(filePath,'PRIMARY')
#    stampBoxMasked[ii]=img[720:926,17:223]#[723:923,20:220]
##%%
#for nn in np.arange(0,len(stampBoxMasked),100):
#    quick_imshow_savepngmv(stampBoxMasked[nn,:,:],'Masked-Aligned Output '+str(nn),'pngdump/mask-trns',nn,'BrBG',0)
##%%
#pixel_by_pixel(stampBoxMasked[:,33:38,158:163],cmap='BrBG',ylim=(0,10000),save=1,savestr='xtrns-m')
##%% extract the subtracted stamp
#stampBoxSubt = np.empty([len(fid),206,206])#,200,200])
#
##%Stamp loc [723:923,20:220]
#for ii, f in enumerate(fid):
#    
#    #Create the full path string
#    fileStr = f+'-subtr.fits'
#    filePath = filesPath+'/'+fileStr
#    
#    #Pull imag for each file
#    [hdr,img] = getvar(filePath,'PRIMARY')
#    stampBoxSubt[ii]=img[720:926,17:223]
##%%
#pixel_by_pixel(stampBoxSubt[:,33:38,158:163],cmap='BrBG',ylim=(0,3000),save=1,savestr='subt')
##%%
#for nn in np.arange(0,len(stampBoxSubt),100):
#    quick_imshow_savepngmv(stampBoxSubt[nn,:,:],'Subt Output '+str(nn),'pngdump/subt',nn,'BrBG',cmin=-3500,cmax=3500,rtnfig=0)
#    #%% extract the convovled stamp
#stampBoxConv = np.empty([len(fid),206,206])#,200,200])
#
##%Stamp loc [723:923,20:220]
#for ii, f in enumerate(fid):
#    
#    #Create the full path string
#    fileStr = f+'-conv.fits'
#    filePath = filesPath+'/'+fileStr
#    
#    #Pull imag for each file
#    [hdr,img] = getvar(filePath,'PRIMARY')
#    stampBoxConv[ii]=img[720:926,17:223]
##%%
#pixel_by_pixel(stampBoxConv[:,33:38,158:163],cmap='BrBG',ylim=(0,10000),save=1,savestr='conv')
##%%
#for nn in np.arange(0,len(stampBoxConv),100):
#    quick_imshow_savepngmv(stampBoxConv[nn,:,:],'Conv Output '+str(nn),'pngdump/conv',nn,'BrBG',0)
#    
##%% extract the OG stamp
#stampBoxOG = np.empty([len(fid),200,200])
#
##%Stamp loc [723:923,20:220]
#for ii, f in enumerate(fid):
#    
#    #Create the full path string
#    fileStr = f+'.fits'
#    filePath = filesPath+'/'+fileStr
#    
#    #Pull imag for each file
#    [hdr,img] = getvar(filePath,'IMAGE')
#    stampBoxOG[ii]=img[723:923,20:220]
##%%
#for nn in np.arange(0,len(stampBoxOG),100):
#    quick_imshow_savepngmv(stampBoxOG[nn,:,:],'Original Field '+str(nn),'pngdump/og',nn,'BrBG',0)
#    
##%% show that the drift if < 1 pixel
#for nn in np.arange(0,len(stampBox),100):
#    quick_imshow_savepngmv(np.log(np.abs(stampBoxMasked[nn,3:153,52:202])),'Field '+str(nn),'pngdump/cadStamp150bwlog',nn,'binary_r',0,cmin=0,cmax=np.log(10000))
#    
##%26:46,153:173
#for nn in np.arange(0,len(stampBox),100):
#    quick_imshow_savepngmv(stampBox[nn,33:38,158:163],'Field '+str(nn),'pngdump/cadStampSmoll',nn,'BrBG',0,cmin=0,cmax=10000)
##%%
#pixel_by_pixel(stampBox[:,33:38,158:163],cmap='BrBG',ylim=(0,10000),save=0)
##%% FISTAR output tracking changing this to plot the average centroid 
#'''
#Plotting the average centroid of the region within the stamp 
#(previously plotting one bright point in y:32:40,x:156:165)
#'''
#xmin=20+153#+25#
#xmax=20+162#+203#-25
#ymin=723+29#+25#
#ymax=723+37#+203#-25
#siz = len(fid)#np.size(stampBoxNonNan,0)
#ID = np.empty(siz)
#X = np.empty(siz)
#Y = np.empty(siz)
#
#
#for fidx in np.arange(len(fistarList)):
#    df = fistar2df(filesPath+'/'+fistarList[fidx])
#    if df.empty != 1:
#        dfSS = df[(df['X'] >xmin) & (df['X'] <xmax) & (df['Y'] >ymin) & (df['Y'] <ymax)]
#        
##        print(dfSS)
#        ID[fidx] = int(dfSS.index[0])
##        X[fidx]  = dfSS['X'].mean(axis=0)
#        X[fidx]  = float(dfSS['X'][dfSS.index[0]])
##        Y[fidx]  = dfSS['Y'].mean(axis=0)
#        Y[fidx]  = float(dfSS['Y'][dfSS.index[0]])
#
#
##%%
#cmapp = 'cool'
#cm = plt.cm.get_cmap()
#xy = range(len(Y))
#z = xy
#        
#plt.title('fistar One Star X Location, Drift')
#plt.scatter(np.arange(len(X)),X,marker='.',c=z,s=.5,cmap=cmapp)
#plt.ylim(177.7,178.1)#(120,122.5)
#plt.ylabel('fistar Centroid Location')
#plt.xlabel('Frame Number, Cadence')
##plt.savefig('fistarOneStarXDrift.png',dpi=200,bbox_inches='tight')
#plt.show()
##%
#plt.title('fistar One Star Y Location, Drift')
#plt.scatter(np.arange(len(Y)),Y,marker='.',c=z,s=.5,cmap=cmapp)
#plt.ylim(755.1,755.5)#(792.6,793.2)#(824.5,827)
#plt.ylabel('fistar Centroid Location')
#plt.xlabel('Frame Number, Cadence')
##plt.savefig('fistarOneStarYDrift.png',dpi=200,bbox_inches='tight')
#plt.show()
#
##%
#cmapp = 'cool'
#cm = plt.cm.get_cmap()
#xy = range(len(Y))
#z = xy
#
#plt.title('fistar One Star Location, Drift')
#plt.scatter(X,Y,c=z,marker='.',s=.5,cmap=cmapp)
#plt.colorbar()
#plt.ylim(755.1,755.5)#(792.6,793.2)#(824.5,827)
#plt.xlim(177.7,178.1)#(201,203.5)#120,122.5)
#plt.ylabel('fistar Y Location')
#plt.xlabel('fistar X Location')
##plt.savefig('fistarOneStarDrift.png',dpi=200,bbox_inches='tight')
#plt.show()
#
#'''
#I want to add a flag for each frame with a QFLAG and then a flag for each 0 point that isn't shown here
#'''
#
##%% Create a fistar file for each masked and transformed frame to show that change in the drift
##% loop for fistar
##toggle_runcmd = 1
#fistarListShift = []
#for ii, f in enumerate(fid):
#    infile = filesPath+'/'+f+'-xtrns-m.fits'
#    fistarfile = f+'-xtrns.fistar'
#    fistarListShift.append(fistarfile)
#    outfile = filesPath+'/'+fistarfile
#    
#    
#    if toggle_runcmd == 1:
#        cmdline = 'fistar ' + infile + ' -o ' + outfile + ' -s flux --comment --model elliptic --flux-threshold ' + str(fluxthld) +' --algorithm uplink --iterations symmetric=2,general=1 --format id,x,y,fwhm,npix,sigma,delta,flux,s/n,bg,amp,s,d,k,cmax --mag-flux 16.62,30'
#        out = os.system(cmdline)
#        if out != 0:
#            print(osfailstr)
#            
#toggle_runcmd = 0 
##%% box loc [720:926,17:223]
#
#xmin=17+156 
#xmax=17+165
#ymin=720+32
#ymax=720+40
#siz = len(fid)#np.size(stampBoxNonNan,0)
#IDShift = np.empty(siz)
#XShift = np.empty(siz)
#YShift = np.empty(siz)
#
#for fidx in np.arange(len(fistarListShift)):
#    df = fistar2df(filesPath+'/'+fistarListShift[fidx])
#    if df.empty != 1:
#        dfSSShift = df[(df['X'] >xmin) & (df['X'] <xmax) & (df['Y'] >ymin) & (df['Y'] <ymax)]
#        
##        print(dfSS)
#        IDShift[fidx] = int(dfSSShift.index[0])
##        XShift[fidx]  = dfSSShift['X'].mean(axis=0)#
#        XShift[fidx]  = float(dfSSShift['X'][dfSSShift.index[0]])
##        YShift[fidx]  = dfSSShift['Y'].mean(axis=0)
#        YShift[fidx]  = float(dfSSShift['Y'][dfSSShift.index[0]])
#
#
##%%
#cmapp = 'cool'
#cm = plt.cm.get_cmap()
#xy = range(len(YShift))
#z = xy
#        #%
#plt.title('fistar Single Star X Location, Drift after Allignment')
#plt.scatter(np.arange(len(XShift)),XShift,marker='.',c=z,s=.5,cmap=cmapp)
##plt.ylim(-.1,.1)#(118,120.5)
#plt.ylabel('fistar Centroid Location')
#plt.xlabel('Frame Number, Cadence')
##plt.savefig('fistarStarXDriftAlligned.png',dpi=200,bbox_inches='tight')
#plt.show()
##%
#plt.title('fistar Single Star Y Location, Drift after Allignment')
#plt.scatter(np.arange(len(YShift)),YShift,marker='.',c=z,s=.5,cmap=cmapp)
##plt.ylim(-.1,.1)#(826,828.5)
##plt.ylim(-.1,1)
#plt.ylabel('fistar Centroid Location')
#plt.xlabel('Frame Number, Cadence')
##plt.savefig('fistarStarYDriftAlligned.png',dpi=200,bbox_inches='tight')
#plt.show()
##%
#plt.title('fistar Single Star Central Location after Allignment')
#plt.scatter(XShift,YShift,c=z,marker='.',s=.5,cmap=cmapp)
#plt.colorbar()
##plt.ylim(-.1,.1)#(826,828.5)
##plt.xlim(-.1,.1)#(118,120.5)
#plt.ylabel('fistar Centroid Y Location')
#plt.xlabel('fistar Centroid X Location')
##plt.savefig('fistarStarDriftAlligned.png',dpi=200,bbox_inches='tight')
#plt.show()
#
#            
##%% get medians 
#medianStamp = np.median(stampBox,0)
###% median subtracted frames
#stampMedSub = stampBox-medianStamp
##%%
##for ii, f in enumerate(stampMedSub):
##    quick_imshow_savepngmv(f,'Median Subtracted Cadence '+str(fid[ii][-12:-7]),'../pngdump/medSub',ii,'BrBG',cmin=-35000,cmax=35000)
#for ii in np.arange(0,len(stampMedSub),100):
#    quick_imshow_savepngmv(stampMedSub[ii],'Median Subtracted Cadence '+str(fid[ii][-12:-7]),'../pngdump/medSub',ii,'BrBG',cmin=-35000,cmax=35000)
#quick_imshow_savepngmv(medianStamp,'Median Super Stamp','../pngdump/medStamp',0,'BrBG')
#
##%% get mean field 
#meanStamp = np.mean(stampBox,0)
##% mean subtracted frames
#stampMeanSub = stampBox-meanStamp
#
##%%
##for ii, f in enumerate(stampMeanSub):
##    quick_imshow_savepngmv(f,'Mean Subtracted Cadence '+str(fid[ii][-12:-7]),'../pngdump/meanSub',ii,'BrBG',cmin=-35000,cmax=35000)
#for ii in np.arange(0,len(stampMeanSub),100):
#    quick_imshow_savepngmv(stampMeanSub[ii],'Mean Subtracted Cadence '+str(fid[ii][-12:-7]),'../pngdump/meanSub',ii,'BrBG',cmin=-35000,cmax=35000)
#quick_imshow_savepngmv(meanStamp,'Mean Super Stamp','../pngdump/meanStamp',0,'BrBG')
#
##%% get mean field unaligned frame
#meanStampOG = np.mean(stampBoxOG,0)
##% mean subtracted frames
#stampMeanSubOG = stampBoxOG-meanStampOG
##%
##for ii, f in enumerate(stampMeanSubOG):
##    quick_imshow_savepngmv(f,'Unaligned Mean Subtracted Cadence '+str(fid[ii][-12:-7]),'../pngdump/meanSubOG',ii,'BrBG',cmin=-35000,cmax=35000)
#for ii in np.arange(0,len(stampMeanSubOG),100):
#    quick_imshow_savepngmv(stampMeanSubOG[ii],'Unaligned Median Subtracted Cadence '+str(fid[ii][-12:-7]),'../pngdump/meanSubOG',ii,'BrBG',cmin=-35000,cmax=35000)
#quick_imshow_savepngmv(meanStampOG,'Unaligned Mean Super Stamp','../pngdump/meanStampOG',0,'BrBG')


#%% 
'''
Check the magnitude limits by locating the objects 
http://nova.astrometry.net/upload

Then we go into the gaia catalog 
https://gea.esac.esa.int/archive/

'''

#%% Check limiting magnitude step
#'''
#for  each  fits file, load it, load it as a temp file of just the subfield?
#
#I think we want to move this up to the top eventually so the RA, Dec are pulled immediately 
#
#'''
#
##%Stamp loc [723:923,20:220]
#for ff in np.arange(len(filesList)):
#    
#    filesList[idx]
#    filePath = filesPath+'/'+fileStr
#    
#    #Pull imag for each file
#    [hdr,img] = getvar(filePath,'IMAGE',verbose=0)
##%%    stampBox[idx]=img[723:923,20:220]
#    
#    
##make a temp header for the temp file
#tmphdr = hdr
#tmphdr['NAXIS1']=200
#tmphdr['NAXIS2']=200
#outfile = 'current_cadence.fits'
#
##%
#hdu = fits.PrimaryHDU(data=medianField,header=tmphdr)
#hdu.writeto(outfile, overwrite=True)
##%
#wcs_header = run_ast(outfile)
##%%
#'''
#figure out how to query the corr.fits file automatically, probs gonna be os.()
#'''
#tbl = fits.open('corr.fits')[1].data
#tbl.columns
#'''
#[('field_x', '>f8'), ('field_y', '>f8'), ('field_ra', '>f8'), ('field_dec', '>f8'), 
# ('index_x', '>f8'), ('index_y', '>f8'), ('index_ra', '>f8'), ('index_dec', '>f8'),
# ('index_id', '>i4'), ('field_id', '>i4'), ('match_weight', '>f8')]
#'''
##%%
#t = Table()
#t = Table(names=('RA', 'DEC'), \
#              dtype=('f4', 'f4'))
#for i in range(0,len(tbl)):
#    t.add_row((tbl.index_ra[i],tbl.index_dec[i]))
#    
#output='ra_dec.dat'
#t.write(output,format='ascii',overwrite=True)
#
##%%
#df = t.to_pandas()
#df=df.assign(GaiaID=0,TIC=0,Vmag=0,Jmag=0,Kmag=0,kepmag=0,Imag=0,Rmag=0,Bmag=0,Rad=0,Mass=0,Lum=0,logg=0,LClass=0,\
#         Teff=0,EBPRP=0)
#for i in range(0,len(t)):
#    print('Working on source #',i)
#    coords=str(t[i][0])+ ' '+ str(t[i][1])
#    c = coordinates.SkyCoord(coords,unit=('deg','deg'),frame='icrs')
#    result = v.query_region(c, radius=0.0005*u.deg)
#    print(len(result),'total catalogs for source',i)
#    
#    #collect all the values and choose the one with the best precision
#    Vmags, Jmags, Kmags, Imags, Rmags, Bmags,Rads,Masses,Lums,loggs,RVs,Teffs,\
#    =[],[],[],[],[],[],[],[],[],[],[],[]
#    for ii in result:
#        if (df.GaiaID[i]==0) & ('GAIA' in (ii.columns)):
#            df.GaiaID[i]=ii['GAIA'][0]
#        if (df.GaiaID[i]==0) &('GaiaID' in (ii.columns)):
#            df.GaiaID[i]=ii['GAIA'][0]
#        if (df.TIC[i]==0) & ('TIC' in (ii.columns)):
#            df.TIC[i]=ii['TIC'][0]
#        if (df.LClass[i]==0) & ('LClass' in (ii.columns)):
#            df.LClass[i]=str(ii['LClass'][0])
#        if ('Vmag' in (ii.columns)) & (len(Vmags)<6):
#            Vmags.append(ii['Vmag'][0])
#        if ('Jmag' in (ii.columns)) & (len(Jmags)<6):
#            Jmags.append(ii['Jmag'][0])
#        if ('Kmag' in (ii.columns)) & (len(Kmags)<6):
#            Kmags.append(ii['Kmag'][0])
#        if (df.kepmag[0]==0) & ('kepmag' in (ii.columns)):
#            df.kepmag[i]=str(ii['kepmag'][0])
#        if ('Imag' in (ii.columns)) & (len(Imags)<4):
#            Imags.append(ii['Imag'][0])
#        if ('Rmag' in (ii.columns)) & (len(Rmags)<4):
#            Rmags.append(ii['Rmag'][0])
#        if ('Bmag' in (ii.columns)) & (len(Bmags)<4):
#            Bmags.append(ii['Bmag'][0])
#        ###
#        if ('Rad' in (ii.columns)) & (len(Rads)<4):
#            a=str(ii['Rad'][0])
#            if '--' in a:
#                continue
#            else:
#                Rads.append(ii['Rad'][0])
#        if ('Mass' in (ii.columns)) & (len(Masses)<4):
#            a=str(ii['Mass'][0])
#            if '--' in a:
#                continue
#            else:
#                Masses.append(ii['Mass'][0])
#        if ('Lum' in (ii.columns)) & (len(Lums)<4):
#            a=str(ii['Lum'][0])
#            if '--' in a:
#                continue
#            else:
#                Lums.append(ii['Lum'][0])
#        if ('logg' in (ii.columns)) & (len(loggs)<4):
#            a=str(ii['logg'][0])
#            if '--' in a:
#                continue
#            else:
#                loggs.append(ii['logg'][0])
#        if ('Teff' in (ii.columns)) & (len(Teffs)<4):
#            a=str(ii['Teff'][0])
#            if '--' in a:
#                continue
#            else:
#                Teffs.append(ii['Teff'][0])
#    df.Vmag[i]=str(round(np.nanmedian(Vmags),2))
#    df.Jmag[i]=str(round(np.nanmedian(Jmags),2))
#    df.Kmag[i]=str(round(np.nanmedian(Kmags),2))
#    df.Kmag[i]=str(round(np.nanmedian(Kmags),2))
#    df.Imag[i]=str(round(np.nanmedian(Imags),2))
#    df.Rmag[i]=str(round(np.nanmedian(Rmags),2))
#    df.Bmag[i]=str(round(np.nanmedian(Bmags),2))
#    df.Rad[i]=str(round(np.nanmedian(Rads),2))
#    df.Mass[i]=str(round(np.nanmedian(Masses),2))
#    df.Lum[i]=str(round(np.nanmedian(Lums),2))
#    df.logg[i]=str(round(np.nanmedian(loggs),2))
#    df.Teff[i]=str(round(np.nanmedian(Teffs),2))
#    print('--- ')
#df['Vmag'] = df['Vmag'].astype(float)
#df['Kmag'] = df['Kmag'].astype(float)
#df['VminK']=df['Vmag']-df['Kmag']
#df['Jmag'] = df['Jmag'].astype(float)
#df['JminK']=df['Jmag']-df['Kmag']
#t2 = Table.from_pandas(df)
