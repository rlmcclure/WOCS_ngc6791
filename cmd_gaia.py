#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 13:25:00 2020

@author: rlmcclure
"""

#Plot a color-magnitude diagram for 6791 with all the Kepler outputs
csv_gaia = 'gaia_catalog_1594830066364O-result.csv'

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import scipy.signal as sig
from scipy import stats
import pandas as pd

from astropy import coordinates
from astropy import units as u
from astropy.table import Table

#%%
gaia_6791 = pd.read_csv(csv_gaia,index_col=2)
#%%
p1_6791 = pd.read_csv('Group1_6791.csv',index_col=2)
p2_6791 = pd.read_csv('Group2_6791.csv',index_col=2)
p3_6791 = pd.read_csv('Group3_6791.csv',index_col=2)
p4_6791 = pd.read_csv('Group4_6791.csv',index_col=2)
p5_6791 = pd.read_csv('Group5_6791.csv',index_col=2)
#%%
fig, ax = plt.subplots(figsize=(5,7))
ax.set_ylim(.001,10000)
ax.set_yscale('log')
ax.set_xlim(-1,5)
gaia_6791.plot.scatter('bp_rp','lum_val',c='teff_val',colormap='coolwarm_r',s=1,ax=ax)
plt.show()
#plt.savefig('gaia_cmd.png',dpi=300)
#%%
fig, ax = plt.subplots(figsize=(5,7))
ax.set_ylim(21.5,6)
ax.set_xlim(-1,5)
gaia_6791.plot.scatter('bp_rp','phot_g_mean_mag',c='teff_val',colormap='coolwarm_r',s=1,ax=ax)
gaia_6791.plot.scatter('bp_rp','phot_g_mean_mag',c='grey',s=1,alpha=0.01,ax=ax)
gaia_6791.loc[p1_6791.index.values].plot.scatter('bp_rp','phot_g_mean_mag',c=colorList[0],s=50,marker='+',ax=ax,label='Colman & Good RV')
gaia_6791.loc[p2_6791.index.values].plot.scatter('bp_rp','phot_g_mean_mag',c='g',s=50,marker='+',ax=ax,label='Colman & Maybe RV')
gaia_6791.loc[p3_6791.index.values].plot.scatter('bp_rp','phot_g_mean_mag',c=colorList[9],s=50,marker='+',ax=ax,label='Good RV')
gaia_6791.loc[p4_6791.index.values].plot.scatter('bp_rp','phot_g_mean_mag',c=colorList[11],s=50,marker='+',ax=ax,label='Maybe RV')
gaia_6791.loc[p5_6791.index.values].plot.scatter('bp_rp','phot_g_mean_mag',c=colorList[14],s=50,marker='+',ax=ax,label='Colman Rotators')
plt.legend()
#plt.savefig('gaia_cmd_highpriority_b.png',dpi=300)
plt.show()