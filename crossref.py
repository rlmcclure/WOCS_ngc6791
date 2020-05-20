#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 15 23:52:22 2019

@author: rlmcclure
"""

from astropy.coordinates import SkyCoord
from astropy import units as u
import csv
import numpy as np
import sys

#Use astroquery to match instead of her table

    
#%%
# OUR DATA from WIYN
#manually set wiyn header since not in file
#%
def get_paired_indexes(filename):
    '''
    '''
    header_wiyn = ['File Specific Ind','Indx','Gaia g-Band','??','??','RA Hr','RA Min','RA Sec','Dec Hr','Dec Min','Dec Sec','??','Status']
    
    # create file
    ff = filename+'_matches'+'.txt'
    f_new = open(ff,'w+')
    
    with open(filename) as f:
        csv_reader = csv.reader(f, delimiter=' ')
        #count lines, ID where data begins
        line_count_wiyn = 0
        data_row_wiyn = []
        comment_row_wiyn = []
        for row_raw in csv_reader:
            #remove extra spaces
            if len(row_raw) > 2:
                row = []
                for r in np.arange(len(row_raw)):
                    if row_raw[r]!='':    
                        row.append(row_raw[r])
            else:
                row = row_raw
            #check if line is data
            try:
                float(row[0])
                float(row[1])
                data_row_wiyn = np.append(data_row_wiyn,line_count_wiyn)
            #line is header/text
            except ValueError:
                comment_row_wiyn = np.append(comment_row_wiyn,line_count_wiyn)
            line_count_wiyn += 1
    f.close
    
    #print the header
    tt = ('\nThe header values are '+str(header_wiyn)+'\n\n')
    f_new.write(tt)
    #ADD DIALOG TO ASSIGN COLLUMNS FOR RA DEC etc??
    #% columns of interest INCOMING 
    base_col_wiyn = np.zeros((len(data_row_wiyn),1))
    #grab IDs 
    ID_ind_wiyn = 1
    ID_col_wiyn = np.copy(base_col_wiyn)
    #Right Asenscion
    RA_hr_ind_wiyn = 5
    RA_hr_wiyn = np.copy(base_col_wiyn)
    RA_mn_ind_wiyn = 6
    RA_mn_wiyn = np.copy(base_col_wiyn)
    RA_sc_ind_wiyn = 7
    RA_sc_wiyn = np.copy(base_col_wiyn)
    #Declination
    DEC_hr_ind_wiyn = 8
    DEC_hr_wiyn = np.copy(base_col_wiyn)
    DEC_mn_ind_wiyn = 9
    DEC_mn_wiyn = np.copy(base_col_wiyn)
    DEC_sc_ind_wiyn = 10
    DEC_sc_wiyn = np.copy(base_col_wiyn)
    #Gaia G-Band
    g_ind_wiyn = 2
    g_col_wiyn = np.copy(base_col_wiyn)
    
    #% get the RA, DEC, g band, & ID separate
    with open(filename) as f:
        csv_reader = csv.reader(f, delimiter=' ')
        #now let's grab each collumn for each header element
        line_n = 0
        ind = 0
        for row_raw in csv_reader:
            #remove extra spaces
            if len(row_raw) > 13:
                row = []
                for r in np.arange(len(row_raw)):
                    if row_raw[r]!='':    
                        row.append(row_raw[r])
            else:
                row = row_raw
                
            #skip header data
            if any(line_n == comment_row_wiyn):
                line_n += 1
                continue
            else:
                #fill in data
                ID_col_wiyn[ind] = row[ID_ind_wiyn]
                RA_hr_wiyn[ind] = row[RA_hr_ind_wiyn]
                RA_mn_wiyn[ind] = row[RA_mn_ind_wiyn]
                RA_sc_wiyn[ind] = row[RA_sc_ind_wiyn]
                DEC_hr_wiyn[ind] = row[DEC_hr_ind_wiyn]
                DEC_mn_wiyn[ind] = row[DEC_mn_ind_wiyn]
                DEC_sc_wiyn[ind] = row[DEC_sc_ind_wiyn]
                g_col_wiyn[ind] = row[g_ind_wiyn]
                
                line_n +=1
                ind +=1
    f.close
                
    #% Compare and output the correct ID_ind_wiyn for matches
                
    #create lists with paired narrowed down indicies 
    wiyn_ID = []
    usyd_ID = []
    for s in RA_sc_usyd:
        if any(np.round(s,1)==np.round(RA_sc_wiyn,1)):
            #identify new values that have matches starting with RA
            RA_sc_match = np.unique(np.where(np.round(s,1)==np.round(RA_sc_usyd,1))[0])
            for ss in RA_sc_match:
                #extract wiyn indexes that are equal to identified matching new values with Dec
                RA_DEC_wiyn_ind = np.unique(np.where(np.round(DEC_sc_wiyn,1)==np.round(DEC_sc_usyd[ss],1)))
                #narrow again to matching RA values
                if len(RA_DEC_wiyn_ind)>0:
                    RA_DEC_wiyn_ind_match = RA_DEC_wiyn_ind[np.where(np.round(RA_sc_wiyn[RA_DEC_wiyn_ind],1)==np.round(RA_sc_usyd[ss],1))[0]]
                    #check hr and min values
                    if len(RA_DEC_wiyn_ind_match)>0:
                        if np.round(RA_mn_wiyn[RA_DEC_wiyn_ind_match],1)[0][0]!=np.round(RA_mn_usyd[ss],1)[0]:
                            print('RA Min mismatch at usyd ID=',ID_col_usyd[ss],'ind=',ss,' and wyin ID=',RA_DEC_wiyn_ind_match)
                        elif np.round(RA_hr_wiyn[RA_DEC_wiyn_ind_match],1)[0][0]!=np.round(RA_hr_usyd[ss],1)[0]:
                            print('RA Hr mismatch at usyd ID=',ID_col_usyd[ss],'ind=',ss,' and wyin ID=',RA_DEC_wiyn_ind_match)
                        elif np.round(DEC_mn_wiyn[RA_DEC_wiyn_ind_match],1)[0][0]!=np.round(DEC_mn_usyd[ss],1)[0]:
                            print('DEC Min mismatch at usyd ID=',ID_col_usyd[ss],'ind=',ss,' and wyin ID=',RA_DEC_wiyn_ind_match)
                        elif np.round(DEC_hr_wiyn[RA_DEC_wiyn_ind_match],1)[0][0]!=np.round(DEC_hr_usyd[ss],1)[0]:
                            print('DEC Hr mismatch at usyd ID=',ID_col_usyd[ss],'ind=',ss,' and wyin ID=',RA_DEC_wiyn_ind_match)
                        #append succesful matches
                        else:
                            tt = '\n'
                            f_new.write(tt)
                            tt = 'WIYN ID:'+str(ID_col_wiyn[RA_DEC_wiyn_ind_match][0][0])+'   USYD ID ['+str(cls_col_usyd[ss])+']:'+str(ID_col_usyd[ss][0])+'\n'
                            f_new.write(tt)
                            tt = 'Gaia g-Band:'+str(g_col_wiyn[RA_DEC_wiyn_ind_match][0][0])+' Sloan g-Band:'+str(g_col_usyd[ss][0])+'\n'
                            f_new.write(tt)
                            tt = 'RA Hr Values:  '+str(RA_hr_wiyn[RA_DEC_wiyn_ind_match][0][0])+str(RA_hr_usyd[ss][0])+'\n'
                            f_new.write(tt)
                            tt = 'RA Min Values: '+str(RA_mn_wiyn[RA_DEC_wiyn_ind_match][0][0])+str(RA_mn_usyd[ss][0])+'\n'
                            f_new.write(tt)
                            tt = 'RA Sec Values: '+str(RA_sc_wiyn[RA_DEC_wiyn_ind_match][0][0])+str(RA_sc_usyd[ss][0])+'\n'
                            f_new.write(tt)
                            tt = 'Dec Hr Values: '+str(DEC_hr_wiyn[RA_DEC_wiyn_ind_match][0][0])+str(DEC_hr_usyd[ss][0])+'\n'
                            f_new.write(tt)
                            tt = 'Dec Min Values:'+str(DEC_mn_wiyn[RA_DEC_wiyn_ind_match][0][0])+str(DEC_mn_usyd[ss][0])+'\n'
                            f_new.write(tt)
                            tt = 'Dec Sec Values:'+str(DEC_sc_wiyn[RA_DEC_wiyn_ind_match][0][0])+str(DEC_sc_usyd[ss][0])+'\n'
                            wiyn_ID.append(ID_col_wiyn[RA_DEC_wiyn_ind_match][0][0])
                            usyd_ID.append(ID_col_usyd[ss][0])
                            
                    else: #append 0 where no matching locations
                        wiyn_ID.append(0)
                        usyd_ID.append(ID_col_usyd[ss][0])
    f_new.close()
    return(wiyn_ID,usyd_ID)
    
#%%
# USYD Colman Data at https://docs.google.com/spreadsheets/d/1_stWkn0ypppBu8oBca9eZ6LRNLn9x_E0zzrahe2CYCQ/edit?pli=1#gid=1538639369 
#%%
# RUN FUNCTION HERE for 6791
#%%  
# Run This to Set-Up USYD Colman Data 
#%
filename_usyd = '6791 & 6819 possible binaries - 6791.csv'
with open(filename_usyd) as f:
    csv_reader = csv.reader(f, delimiter=',')
    #count lines, ID where data begins
    line_count_usyd = 0
    header_usyd = []
    for row in csv_reader:
        #check if line is data
        try:
            float(row[0])
        #line is header/text
        except ValueError:
            header_end_usyd = line_count_usyd
            if len(header_usyd) == 0:
                header_usyd = row
            data_start_usyd = line_count_usyd+1
        line_count_usyd += 1
f.close

#print the header
print('\nThe header values are',header_usyd)
#ADD DIALOG TO ASSIGN COLLUMNS FOR RA DEC etc??
#% columns of interest INCOMING 
base_col_usyd = np.zeros((line_count_usyd-data_start_usyd,1))
#grab IDs 
ID_ind_usyd = 0
ID_col_usyd = np.copy(base_col_usyd)
#Right Asenscion
RA_ind_usyd = 5
RA_col_usyd = np.copy(base_col_usyd)
#Declination
DEC_ind_usyd = 6
DEC_col_usyd = np.copy(base_col_usyd)
#Sloan G-Band
g_ind_usyd = 2
g_col_usyd = np.copy(base_col_usyd)
#Classification
cls_ind_usyd = 4
cls_col_usyd = []


#% get the RA, DEC, g band, & ID separate
with open(filename_usyd) as f:
    csv_reader = csv.reader(f, delimiter=',')
    #now let's grab each collumn for each header element
    line_n = 0
    for row in csv_reader:
        #skip header data
        if line_n == header_end_usyd:
            line_n += 1
            continue
        else:
            #fill in data
            ind = line_n-data_start_usyd
            ID_col_usyd[ind] = row[ID_ind_usyd]
            RA_col_usyd[ind] = row[RA_ind_usyd]
            DEC_col_usyd[ind] = row[DEC_ind_usyd]
            g_col_usyd[ind] = row[g_ind_usyd]
            cls_col_usyd.append(row[cls_ind_usyd])
            
            line_n +=1
            
            
#% set up columns to compare to 
#Right Asenscion
RA_hr_usyd = np.copy(base_col_usyd)
RA_mn_usyd = np.copy(base_col_usyd)
RA_sc_usyd = np.copy(base_col_usyd)
#Declination
DEC_hr_usyd = np.copy(base_col_usyd)
DEC_mn_usyd = np.copy(base_col_usyd)
DEC_sc_usyd = np.copy(base_col_usyd)

#fill in appropriate data values
for ind in np.arange(line_count_usyd-data_start_usyd):
    c = SkyCoord(RA_col_usyd[ind], DEC_col_usyd[ind], frame='icrs', unit='deg')
    RA_hr_usyd[ind],RA_mn_usyd[ind],RA_sc_usyd[ind] = c.ra.hms
    DEC_hr_usyd[ind],DEC_mn_usyd[ind],DEC_sc_usyd[ind] = c.dec.dms
    
    
(wiyn_ID_6791,usyd_ID_6791) = get_paired_indexes('ngc6791master.sort.sep19')
#%%
# RUN FUNCITON FOR 6819
#%%  
# Run This to Set-Up USYD Colman Data 
#%
filename_usyd = '6791 & 6819 possible binaries - 6819.csv'
with open(filename_usyd) as f:
    csv_reader = csv.reader(f, delimiter=',')
    #count lines, ID where data begins
    line_count_usyd = 0
    header_usyd = []
    for row in csv_reader:
        #check if line is data
        try:
            float(row[0])
        #line is header/text
        except ValueError:
            header_end_usyd = line_count_usyd
            if len(header_usyd) == 0:
                header_usyd = row
            data_start_usyd = line_count_usyd+1
        line_count_usyd += 1
f.close

#print the header
print('\nThe header values are',header_usyd)
#ADD DIALOG TO ASSIGN COLLUMNS FOR RA DEC etc??
#% columns of interest INCOMING 
base_col_usyd = np.zeros((line_count_usyd-data_start_usyd,1))
#grab IDs 
ID_ind_usyd = 0
ID_col_usyd = np.copy(base_col_usyd)
#Right Asenscion
RA_ind_usyd = 5
RA_col_usyd = np.copy(base_col_usyd)
#Declination
DEC_ind_usyd = 6
DEC_col_usyd = np.copy(base_col_usyd)
#Sloan G-Band
g_ind_usyd = 2
g_col_usyd = np.copy(base_col_usyd)
#Classification
cls_ind_usyd = 4
cls_col_usyd = []

#% get the RA, DEC, g band, & ID separate
with open(filename_usyd) as f:
    csv_reader = csv.reader(f, delimiter=',')
    #now let's grab each collumn for each header element
    line_n = 0
    for row in csv_reader:
        #skip header data
        if line_n == header_end_usyd:
            line_n += 1
            continue
        else:
            #fill in data
            ind = line_n-data_start_usyd
            ID_col_usyd[ind] = row[ID_ind_usyd]
            RA_col_usyd[ind] = row[RA_ind_usyd]
            DEC_col_usyd[ind] = row[DEC_ind_usyd]
            g_col_usyd[ind] = row[g_ind_usyd]
            cls_col_usyd.append(row[cls_ind_usyd])
            line_n +=1
            
            
#% set up columns to compare to 
#Right Asenscion
RA_hr_usyd = np.copy(base_col_usyd)
RA_mn_usyd = np.copy(base_col_usyd)
RA_sc_usyd = np.copy(base_col_usyd)
#Declination
DEC_hr_usyd = np.copy(base_col_usyd)
DEC_mn_usyd = np.copy(base_col_usyd)
DEC_sc_usyd = np.copy(base_col_usyd)

#fill in appropriate data values
for ind in np.arange(line_count_usyd-data_start_usyd):
    c = SkyCoord(RA_col_usyd[ind], DEC_col_usyd[ind], frame='icrs', unit='deg')
    RA_hr_usyd[ind],RA_mn_usyd[ind],RA_sc_usyd[ind] = c.ra.hms
    DEC_hr_usyd[ind],DEC_mn_usyd[ind],DEC_sc_usyd[ind] = c.dec.dms
    
(wiyn_ID_6819,usyd_ID_6819) = get_paired_indexes('ngc6819master.sort.sep19')

#        

#%%
get_wiynID = lambda KIC: wiyn_ID_6791[np.where(np.asarray(usyd_ID_6791)==KIC)[0][0]]