
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

import kplr
client = kplr.API()