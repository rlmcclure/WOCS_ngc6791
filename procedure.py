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
Q3, Q7, Q11, Q15 for channel 53 or 1
Q2, Q6, Q10, Q14 for channel 29 
Q1, Q5, Q9,  Q13, Q17 for channel 1 or 53 (Q1 is only 34 days, Q0 is 9 days)

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
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from astropy.io import fits
import scipy.signal as sig
from scipy import stats
import lightkurve as lk
import pandas as pd

from astroquery.astrometry_net import AstrometryNet
ast = AstrometryNet()
ast.api_key = 'dnzncgdyolxsukyy'

from astroquery.mast import Catalogs
from astroquery.vizier import Vizier
v = Vizier()

from astropy import coordinates
from astropy import units as u
from astropy.table import Table
#%%
toggle_runcmd = 0
osfailstr = '\nheyhi so that command failed... plz investigate.\n'
#%%
filesPath = '/Users/kookoo2052/HOLD/mosaic-ch81/q04c81'#q4
#filesPath = '/Volumes/rmathieu/NGC6791/q8c81'

filesList = [f for f in os.listdir(filesPath) if (f.endswith(".fits") and os.path.isfile(os.path.join(filesPath, f)) and not f.endswith("-xtrns.fits"))]
filesList.sort()
#%% loop for fistar
fluxthld = 5

fistarList = []
for ii, f in enumerate(filesList[:]):
    infile = filesPath+'/'+f
    fistarfile = f[:-5]+'.fistar'
    fistarList.append(fistarfile)
    outfile = filesPath+'/'+fistarfile
    
    if toggle_runcmd == 1:
        cmdline = 'fistar ' + infile + ' -o ' + outfile + ' -s flux --comment --model elliptic --flux-threshold ' + str(fluxthld) +' --algorithm uplink --iterations symmetric=2,general=1 --format id,x,y,fwhm,npix,sigma,delta,flux,s/n,bg,amp,s,d,k,cmax --mag-flux 16.62,30'
        os.system(cmdline)
        if out != 0:
            print(osfailstr)
    

#    for i in *.fits; do fistar "$i" -o $(echo $i | sed s/fits/fistar/)  -s flux 
#    --comment --model elliptic --flux-threshold 500 --algorithm uplink 
#    --iterations symmetric=2,general=1 
#    --format id,x,y,fwhm,npix,sigma,delta,flux,s/n,bg,amp,s,d,k,cmax
#    --mag-flux 16.62,30; done
#%
#fistarList = [f for f in os.listdir(filesPath) if (f.endswith(".fistar") and os.path.isfile(os.path.join(filesPath, f)))]
#fistarList.sort()

#%% find minimum FWHM fistar file for astroref.
fid = []
fidx = []
fidCut = []
fidxCut = []
nCut = 0

for x in np.arange(len(fistarList)):
    df = fistar2df(filesPath+'/'+fistarList[x])
    if df.empty != 1:
        fidx.append(x)
        fid.append(fistarList[x])
    else:
        fidxCut.append(x)
        fidCut.append(fistarList[x])
        nCut += 1
#%%
siz=len(fid) #fileIDs of good fields
fwhmMean = np.empty(siz)
fwhmMin = np.empty(siz)
fwhmMax = np.empty(siz)
fwhmStd = np.empty(siz)
#%
for ii, x in enumerate(fidx):
    df = fistar2df(filesPath+'/'+fistarList[x])
    fwhmMean[ii] = np.mean(df['FWHM'])
    fwhmMin[ii] = np.amin(df['FWHM'])
    fwhmMax[ii] = np.amax(df['FWHM'])
    fwhmStd[ii] = np.std(df['FWHM'])
    
#%
dfFWHM = pd.DataFrame(data={'mean':fwhmMean,'std':fwhmStd,'min':fwhmMin,'max':fwhmMax},index=fid)
del(fwhmMean,fwhmStd,fwhmMin,fwhmMax)
#%
print('\nMinimum FWHM:\n',dfFWHM.loc[dfFWHM['min'].idxmin])
print('\nMinimum Mean FWHM:\n',dfFWHM.loc[dfFWHM['mean'].idxmin])
#% Make copy of minimum FWHM file as astroref file
astroref = filesPath+'/'+'astroref.fistar'

#% toggle to run commands
if toggle_runcmd == 1:
    cmdline = 'cp ' + filesPath+'/'+dfFWHM.loc[dfFWHM['mean'].idxmin].name + ' ' + astroref
    out = os.system(cmdline)
    if out != 0:
        print(osfailstr)
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

#loop for grmatch
for ii, f in enumerate(fid):
    infile = filesPath+'/'+f
    transfile = infile[:-7]+'.itrans'
    trid.append(transfile)
        
    #% toggle to run commands
    if toggle_runcmd == 1:
        cmdline = 'grmatch --match-points -r ' + astroref + ' --col-ref 2,3 --col-ref-ordering +8 -i ' + infile + ' --col-inp 2,3 --col-inp-ordering +8 --weight reference,column=12 --triangulation maxinp=1000,maxref=1000,auto,unitarity=0.01 --order 2 --max-distance 1 --comment --output-transformation ' + transfile + ' --output /dev/null'
        out = os.system(cmdline)
        if out != 0:
            print(osfailstr)

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
    infile = filesPath+'/'+f[:-7]+'.fits'
    transfile = infile[:-5]+'.itrans'
    xtrnsfile = infile[:-5]+'-xtrns.fits'
    
    #% toggle to run commands
    if toggle_runcmd == 1:
        cmdline = 'fitrans ' + infile + ' -k --input-transformation ' + transfile + ' --reverse -o ' + xtrnsfile 
        out = os.system(cmdline)
        if out != 0:
            print(osfailstr)

#%%Extract the aligned stamp

stampBox = np.empty([len(fid),206,206])#,200,200])

#%Stamp loc [723:923,20:220]
for ii, f in enumerate(fid):
    
    #Create the full path string
    fileStr = f[:-7]+'-xtrns.fits'
    filePath = filesPath+'/'+fileStr
    
    #Pull imag for each file
    [hdr,img] = getvar(filePath,'PRIMARY',verbose=0)
    stampBox[ii]=img[720:926,17:223]#[723:923,20:220]
    
#% extract the OG stamp
stampBoxOG = np.empty([len(fid),200,200])

#%Stamp loc [723:923,20:220]
for ii, f in enumerate(fid):
    
    #Create the full path string
    fileStr = f[:-7]+'.fits'
    filePath = filesPath+'/'+fileStr
    
    #Pull imag for each file
    [hdr,img] = getvar(filePath,'IMAGE',verbose=0)
    stampBoxOG[ii]=img[723:923,20:220]

#%% show that the drift if < 1 pixel
for nn in np.arange(0,len(stampBox),100):
    quick_imshow_savepngmv(np.log(np.abs(stampBox[nn,3:153,52:202])),'Field '+str(nn),'cadStamp150bwlog',nn,'binary_r',0,cmin=0,cmax=np.log(10000))
    
#%26:46,153:173
for nn in np.arange(0,len(stampBox),100):
    quick_imshow_savepngmv(stampBox[nn,33:38,158:163],'Field '+str(nn),'cadStampSmoll',nn,'BrBG',0,cmin=0,cmax=10000)
#%%
pixel_by_pixel(stampBox[:,33:38,158:163],cmap='BrBG',ylim=(0,10000))
#%% FISTAR output tracking y:23:43,x:150:170
xmin=20+156
xmax=20+165
ymin=723+32
ymax=723+40
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
#%
plt.title('fistar X Drift')
plt.scatter(np.arange(len(X)),X,marker='.',c='k',s=.5)
plt.ylim(177,179)
plt.ylabel('fistar Centroid Location')
plt.xlabel('Frame Number, Cadence')
#plt.savefig('fistarXDrift.png',dpi=200,bbox_inches='tight')
plt.show()
#%
plt.title('fistar Y Drift')
plt.scatter(np.arange(len(Y)),Y,marker='.',c='k',s=.5)
plt.ylim(755,756)
plt.ylabel('fistar Centroid Location')
plt.xlabel('Frame Number, Cadence')
#plt.savefig('fistarYDrift.png',dpi=200,bbox_inches='tight')
plt.show()

#%% photometric reference frame USE FICOMBINE
#%% get medians 
medianStamp = np.median(stampBox,0)
##% median subtracted frames
stampMedSub = stampBox-medianStamp
#%%
#for ii, f in enumerate(stampMedSub):
#    quick_imshow_savepngmv(f,'Median Subtracted Cadence '+str(fid[ii][-12:-7]),'pngdump/medSub',ii,'BrBG',cmin=-35000,cmax=35000)
for ii in np.arange(0,len(stampMedSub),100):
    quick_imshow_savepngmv(stampMedSub[ii],'Median Subtracted Cadence '+str(fid[ii][-12:-7]),'pngdump/medSub',ii,'BrBG',cmin=-35000,cmax=35000)
quick_imshow_savepngmv(medianStamp,'Median Super Stamp','pngdump/medStamp',0,'BrBG')

#%% get mean field 
meanStamp = np.mean(stampBox,0)
#% mean subtracted frames
stampMeanSub = stampBox-meanStamp

#%%
#for ii, f in enumerate(stampMeanSub):
#    quick_imshow_savepngmv(f,'Mean Subtracted Cadence '+str(fid[ii][-12:-7]),'pngdump/meanSub',ii,'BrBG',cmin=-35000,cmax=35000)
for ii in np.arange(0,len(stampMeanSub),100):
    quick_imshow_savepngmv(stampMeanSub[ii],'Mean Subtracted Cadence '+str(fid[ii][-12:-7]),'pngdump/meanSub',ii,'BrBG',cmin=-35000,cmax=35000)
quick_imshow_savepngmv(meanStamp,'Mean Super Stamp','pngdump/meanStamp',0,'BrBG')

#%% get mean field unaligned frame
meanStampOG = np.mean(stampBoxOG,0)
#% mean subtracted frames
stampMeanSubOG = stampBoxOG-meanStampOG
#%
#for ii, f in enumerate(stampMeanSubOG):
#    quick_imshow_savepngmv(f,'Unaligned Mean Subtracted Cadence '+str(fid[ii][-12:-7]),'pngdump/meanSubOG',ii,'BrBG',cmin=-35000,cmax=35000)
for ii in np.arange(0,len(stampMeanSubOG),100):
    quick_imshow_savepngmv(stampMeanSubOG[ii],'Unaligned Median Subtracted Cadence '+str(fid[ii][-12:-7]),'pngdump/meanSubOG',ii,'BrBG',cmin=-35000,cmax=35000)
quick_imshow_savepngmv(meanStampOG,'Unaligned Mean Super Stamp','pngdump/meanStampOG',0,'BrBG')


#%% fcombine fails because there are too many files?
toggle_runcmd = 1 
infile = filesPath+'/'+'*xtrns.fits'
outfile = filesPath+'/'+'photoref.fits'
#% toggle to run commands
if toggle_runcmd == 1:
    cmdline = 'ficombine ' + infile + ' --mode mean -o ' + outfile 
    out = os.system(cmdline)
    if out != 0:
        print(osfailstr)

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
