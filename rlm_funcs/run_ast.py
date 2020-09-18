#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 15:41:48 2020

@author: rlm
"""
import os

import numpy as np

from astroquery.astrometry_net import AstrometryNet
ast = AstrometryNet()
ast.api_key = 'dnzncgdyolxsukyy'

from astroquery.mast import Catalogs
from astroquery.vizier import Vizier
v = Vizier()

from astropy import coordinates
from astropy import units as u
from astropy.table import Table

import kplr
client = kplr.API()
#%% Astrometry query INTEGRATE THIS INTO THE CODE 
def run_ast(outfile):
    '''
    runs the astrometry on a field NOT COMPLETE
    '''
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