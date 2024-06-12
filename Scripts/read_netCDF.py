"""
Script to read
- NetCDF files and transpose to dataframe
- Monthly mean of Temperature and primary production at PAP-SO
"""
import pandas as pd
## IMPORTS
import os
from scipy.interpolate import griddata
from parameters import watercolumn, dz
import xarray as xr
import numpy as np
import cmocean
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
from pyhdf.SD import SD, SDC
import xarray as xr
import rioxarray as rxr
import pprint


path_to_plot = Path('~/simul_dvm2/Plots').expanduser()
path_to_data = Path('~/simul_dvm2/Data').expanduser()
os.chdir(path_to_data)

"""
1) Temperature from PAP station (48.75°N, 16.17°W) : calculation of the monthly temperature along depth
"""

arr = xr.open_dataset('~/simul_dvm2/Data/cmems_mod_glo_phy-thetao_anfc_0.083deg_PT6H-i_1681919076982.nc')
df_temp = arr.to_dataframe()
df_temp = df_temp.reset_index()
df_temp['time'] = df_temp['time'].dt.month #convert time to day of the year(1-365)
df_temp = df_temp.groupby(['time', 'depth'])['thetao'].mean().reset_index(name='daily_T')

xt = np.arange(744, 8760, 1)
yi = watercolumn
xi, yi = np.meshgrid(xt, yi)
x, y, z = df_temp.time * 24 * 30.5, df_temp.depth, df_temp.daily_T

annualT = griddata((x, y), z, (xi, yi), method='linear')
#Remove NaNs above -0.5m
annualT[np.isnan(annualT)] = 0

fig = plt.figure(figsize=(4, 4), layout="constrained")
CS = plt.contourf(xi/24/30.5, -yi, annualT, cmap=cmocean.cm.thermal)
cbar = fig.colorbar(CS, location='top')
plt.xlabel('Months')
cbar.set_label("T [°C]")
plt.ylim(-800, 0)
plt.xlim(1, 12)
plt.xticks([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], ["Jan", "Fev", "Mar", "Apr", "May", "June", "July", "Aug", "Sep", "Oct", "Nov", "Dec"], size=10)
plt.tick_params(axis='x', labelrotation=45)
plot_temp = path_to_plot / "Temp_annual_PAP.pdf"
plt.savefig(plot_temp)
plt.close()

"""
2) PP at PAP station (48-49°N, 16-17°W) : calculation of the monthly PP along depth during 2022
From cmems_mod_ibi_bgc_my_0.083deg-3D_P1M-m in mgC m-3
"""

arr = xr.open_dataset('~/simul_dvm2/Data/cmems_mod_ibi_bgc_my_0.083deg-3D_P1M-m_1716798448048.nc')
df_nppv = arr.to_dataframe()
df_nppv = df_nppv.reset_index()
df_nppv['month'] = df_nppv['time'].dt.month
df_nppv['phyc'] = df_nppv['phyc']*12 #conversion fo mmolC in mgC
df_temp = df_nppv.groupby(['month', 'depth'])['phyc'].mean().reset_index(name='daily_nppv')
x, y, z = df_temp.month * 24 * 30.5, df_temp.depth, df_temp.daily_nppv

annualNPP = griddata((x, y), z, (xi, yi), method='linear')
#Remove NaNs above -0.5m
annualNPP[np.isnan(annualNPP)] = 0
#Integrate the NPP along depth
int_annualNPP = np.max(annualNPP, axis=0)
sum_annualNPP = np.sum(annualNPP, axis=0)*dz

fig = plt.figure(figsize=(4, 4), layout="constrained")
CS = plt.contourf(xi/24/30.5, -yi, annualNPP, cmap=cmocean.cm.algae)
cbar = fig.colorbar(CS, location="top")
plt.xlim(1, 12)
plt.xticks([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], ["Jan", "Fev", "Mar", "Apr", "May", "June", "July", "Aug", "Sep", "Oct", "Nov", "Dec"], size=10)
plt.xlabel('Months')
plt.ylabel('Depth [$m$]')
plt.ylim(-200, 0)
plt.tick_params(axis='x', labelrotation=45)
#plt.xlim(1, 12)
cbar.set_label("$mgC$ $m^{-3}$")
plot_temp = path_to_plot / "NPP_annual_PAP.pdf"
plt.savefig(plot_temp)
plt.close()



