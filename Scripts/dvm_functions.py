"""
Functions used to parametrize the DVM model
"""
import numpy as np
import math
from pathlib import Path
from numba import jit

path_to_plot = Path('~/model_dvm_carbon/Plots').expanduser()


## Irradiance function from J-C Poggiale.
def I_sinusF(t, n=25, a=200):
    T = 48
    omega = (2 * np.pi) / T
    I0 = 1 - (np.exp(-a * np.sin(omega * t)**n)) + 0.5
    return I0

def I_sinus_annualF(t):
    n_max = 60
    n_min = 14
    T = 365
    omega = (2 * np.pi) / T
    phi = T
    Dmig = (n_max-n_min) * (np.cos(omega*(t+phi))+1)/2 + n_min
    return Dmig

@jit(nopython=True)
def beer_lambert(I0, z, k=0.001):
    I_z = I0 * np.exp(-k*z)
    return I_z

@jit(nopython=True)
def I_depthF(Nt12, dt_inf, zz):
    I0 = np.zeros(Nt12)
    for t in range(Nt12):
        I0[t] = I_richards(t * dt_inf + dt_inf)
    I0 = np.hstack((I0, np.flip(I0)))
    # Create a 2d array (temps*depth) to represent the irradiance along depth
    I_depth = np.zeros((zz, Nt12 * 2))
    for z in range(zz):
        I_depth[z] = beer_lambert(I0, z)
    return I_depth

## Dictionary that references the parameters of the functional groups : (0)exponent, (1)multiplicator, (2)factor of conversion to ww in mg, (3)conversion in cw in mgC,
# (4)min length in mm, (5)max length in mm, (6) V/L, (7)gut transit in h-1, (8-12) parameters for respi_func from Ikeda 2016 (ceph and fish) and Ikeda 2013 (crust), (13) conversion efficiency, (14) taxonomic ratio, (15)
#                                   0        1       2          3     4     5    6     7      8       9      10     11    12      13   14  15
micron_dict = {    "bm_ceph":    [2.611,   -3.5,        1e3,  55.44,  15,   50, 1.25,  0.5, 24.461, 0.868, -6.424, -0.261, 0.95, 0.9, 0.2, 0.1], # ww in g to cw (conversion from g to mg)
                  "bm_crust":    [3.23,  -3.261, (100/18.3),  0.419,  15,   50,  1.1, 0.67, 23.079, 0.813, -6.248, -0.136, 0.97, 0.8, 0.15, 0.2], # dw in mg to cw for Euphausiids (R2=0.98)
                    "bm_fish":   [2.902, -1.797,          1,  0.092,  15,   50,  1.8,  0.2, 30.767, 0.870, -8.515, -0.088,  0.9, 0.8, 0.65, 0.16]} # ww in mg to cw for Myctophidae (R2=0.935)

def length_weight_func(length_mm, bm_id):
    # Gives the weight of micronekton in mg and mgC from the length in mm

    # References in Kwong et al. 2020
    a = micron_dict[bm_id][0]
    b = micron_dict[bm_id][1]
    rw = pow(10, a * np.log10(length_mm) + b)
    ww = rw * micron_dict[bm_id][2] #in mg
    cw = rw * micron_dict[bm_id][3] #in mgC
    return ww, cw, rw

### METABOLIC FUNCTIONS ###

#Respiration rate from Ikeda 2014 and 2016
@jit(nopython=True)
def RMR_func(CW, T, D, a0, a1, a2, a3, RQ): #CW in mgC, T in K, D for depth in m
    if D == 0:
        ln_RO = 0
    else:
        ln_RO = a0 + a1 * np.log(CW) + a2 * (1000/T) + a3 * np.log(D) #in µL O2 h-1
    Rmc = np.exp(ln_RO) * RQ * (12/22.4) #in µgC h-1
    R = Rmc/(CW*1000) #in h-1
    return R

#Van't Hoff equation for the growth rate of copepods
def vant_hoff(Tc):
    Qdix = 3
    Tref = 15 #reference temperature in °C
    Rref = 0.008 #reference growth rate in h-1
    r = Qdix ** ((Tc - Tref) / 10) * Rref
    return r

### ALLOMETRIC FUNCTIONS
# Allometric determination of the swimming speed in cm/s
def swimming_speed_func(BL, bm_id):
    a = micron_dict[bm_id][6]
    V = a * BL
    V = (V / 100) * 3600
    return V

def speedF(Nt12, dt_inf, vmax):
    dI_r_max = 0
    dI_r = np.zeros(Nt12)

    It_1 = np.zeros(Nt12)
    It = np.zeros(Nt12)
    for t in range(Nt12):
        It_1[t] = I_sinusF(t * dt_inf + dt_inf)
        It[t] = I_sinusF(t * dt_inf)

    for t in range(Nt12):
        if It[t] <= 1e-3:
            dI_r[t] = 0
        else:
            dI_r[t] = (1 / dt_inf) * (It_1[t] / It[t] - 1)
        if np.abs(dI_r[t]) > dI_r_max:
            dI_r_max = np.abs(dI_r[t])
    dI_r = np.hstack((dI_r, -np.flip(dI_r)))

    # Computation of the speed
    V = np.zeros(2 * Nt12)
    plot_time24 = []
    for t in range(2 * Nt12):
        V[t] = vmax * dI_r[t] / dI_r_max
        plot_time24.append(t)
    return V #m/h


@jit(nopython=True)
def remineratization_rateF(Tc): #in h-1
    r = (0.005*Tc + 0.011)*24
    return r

# Zhang and Dam (1997) adaption of Peterson and Wroblewski (1984) used to estimate natural mortality in day-1
def mortality_allometricF(DW, mort_coef=5.26): #DW in g, B in mmol C m-2 and t number of hours per day that diel migrant stay below the euphotic zone ( about 12 h d-1).
    DM = (mort_coef*1e-3) * DW**-0.25
    M = DM/24
    return M

def natural_mortality_L96(WW): #From Lorenzen 1996 for fish in ocean
    M = 3.69 * (WW**-0.305)
    M = M/365/24
    return M

def natural_mortality_L22(L, L_inf=10.24, K=1): #From Lorenzen 2022
    # median values of myctophidae from fishbase
    a = 0.28
    b = -1.30
    c = 1.08
    lnM = a + b * np.log(L/L_inf) + c * np.log(K)
    Man = np.exp(lnM)
    M = Man / 365 / 24  # mortality in h-1
    return M

def mortality_rate_pauly_func(T): #L_inf in cm, K in year-1 and T in °C
    L_inf = 10.24
    K = 0.532 #median values of myctophidae from fishbase
    lnM = -0.0066 - 0.279*np.log(L_inf) + 0.6543*np.log(K) + 0.4634*np.log(T)
    Man = np.exp(lnM)
    M = Man/365/24 #mortality in h-1
    return M

# Get the current date and time
def sunrise_sunset_calc(date_, lat, lon):
    from datetime import time
    from datetime import datetime
    import ephem
    """input needs to be date and time formats from datetime, output in UTC"""
    time_ = time(12,0)
    obs = ephem.Observer()
    #PyEphem takes and returns only UTC times.
    obs.date = str(date_) + " " + str(time_)
    d_d = datetime.combine(date_, time_)
    #Location of observer
    obs.lon = str(lon) #Note that lon should be in string format
    obs.lat = str(lat)      #Note that lat should be in string format
    #Elevation of observer
    obs.elev = 0
    #To get U.S. Naval Astronomical Almanac values, use these settings
    obs.pressure = 0
    #obs.horizon = '-0:34'
    sunrise = obs.previous_rising(ephem.Sun()) #Sunrise
    noon = obs.next_transit(ephem.Sun(), start=sunrise) #Solar noon
    sunset = obs.next_setting(ephem.Sun()) #Sunset
    return sunrise.datetime().time(), sunset.datetime().time(), noon.datetime().time()

def convert_float_time_format(float_sec):
    hours, seconds = divmod(float_sec, 3600)  # split to hours and seconds
    minutes, seconds = divmod(seconds, 60)
    duree_mig = "{:02.0f}:{:02.0f}".format(hours, minutes)
    return duree_mig

# Function to calculate moving average
def windowed_mean(arr, n):
    import scipy.signal as sig
    dims = len(arr.shape)
    s = sig.convolve(arr, np.ones((2*n+1,)*dims), mode='same')
    d = sig.convolve(np.ones_like(arr), np.ones((2*n+1,)*dims), mode='same')
    return s/d

#Define the Gaussian function
@jit(nopython=True)
def gaussian(x_array, amp, cen, sigma):
    #return amp*(1/(alpha*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x_array-cen)/alpha)**2)))
    return amp * np.exp(-(x_array-cen)**2/(2*sigma**2))

#Function to normalize data inside an array
def normalize(x):
    norm = (x - np.amin(x)) / (np.amax(x) - np.amin(x))
    return norm

# A revoir
def zenith_angleF(lon, lat, J, t):
    SD = -23.44 * np.cos((math.radians(360)/365)*(J+10))
    f = 279.575 + 0.9856 * J
    ET = (-104.7 * math.sin(f) + 596.2*math.sin(2*f) + 4.3 * math.sin(3*f) - 12.7*math.sin(4*f) - 429.3*math.cos(f) - 429.3*math.cos(2*f) + 19.3*math.cos(3*f))/3600
    t0 = 12 - 1/15*lon - ET
    cos_ZA = math.sin(lat) * math.sin(SD) + math.cos(lat) * math.cos(SD) * math.cos(15*(t-t0))
    ZA = math.degrees(np.arccos(cos_ZA))
    return ZA
