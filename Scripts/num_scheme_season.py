"""
Script of the numerical scheme for the model of carbon production by migrant micronekton with annual variation fo the carrying capacity
"""

## IMPORTS
from numba import jit
from parameters import *
from dvm_functions import remineratization_rateF, RMR_func, SDA_func, beer_lambert

#Define the Gaussian function
@jit(nopython=True)
def gaussian(x_array, amp, cen, sigma):
    return amp * np.exp(-(x_array-cen)**2/(2*sigma**2))


# The numerical scheme for the model for dvm_model_seasonal_v1 with variation fo the carrying capacity
@jit(nopython=True)
def AN_RGCD_K2(R0, C0, Wc_ind, V, TkL, d, coef_alpha, I_norm, beta, dt_inf, NPP):
    # Initialisation
    R = np.zeros((tt, zz))
    C = np.zeros((tt, zz))
    G = np.zeros((tt, zz))

    Dm = np.zeros((tt, zz))
    Dg = np.zeros((tt, zz))

    R[0] = R1 = np.copy(R0)
    C[0] = C1 = np.copy(C0)

    G1 = np.zeros(zz)
    D1g = np.zeros(zz)
    D1m = np.zeros(zz)

    t = 0 #time in hour
    nh = 0 #number of hours
    nbhy = 334*24 #number of hours in a year
    rmr_last = RMR_func(Wc_ind, TkL[-1], zz * dz, a0, a1, a2, a3, RQ)

    ## Loop over temps
    while t < tmax - 2 * dt_save:
        R2 = np.zeros(zz)
        C2 = np.zeros(zz)
        G2 = np.zeros(zz)
        w = V[int((t % 24) / dt_inf)]
        if w == 0:
            dt = dt_sup
        else:
            dt = min(dt_sup, dz / np.abs(w) * CFL)
        nu = w * dt / dz

        alphaZ = I_norm[:, int((t % 24) / dt_inf)] * coef_alpha
        Kt = NPP[:, int(t % nbhy)]

        ## Loop over space
        for z in range(zz):
            depth = z * dz
            # Maintenance rate : RMR + AMR + excretion
            rmr = RMR_func(Wc_ind, TkL[z], depth, a0, a1, a2, a3, RQ) # mgC ind-1 h-1

            if bm_id == "bm_crust":
                m = (1 + 0.31 + (np.abs(w) / w_ref))
            else:
                m = (1 + (np.abs(w) / w_ref))

            D1m[z] += dt * rmr * (C1[z]+G1[z]) * m + dt * e * d * G1[z] * Csda
            D1g[z] += dt * d * (1-e) * G1[z]

            if Kt[z] != 0:
                R2[z] = R1[z] + dt * (rho * R1[z] * (1 - R1[z] / Kt[z])) - dt * alphaZ[z] * R1[z] * C1[z] / (1 + beta * R1[z])

        for z in range(1, zz - 1):
            rmr = RMR_func(Wc_ind, TkL[z], z * dz, a0, a1, a2, a3, RQ)
            if w >= 0:
                G2[z] = (1 - nu) * G1[z] + nu * G1[z - 1] + dt * alphaZ[z] * R1[z] * C1[z] / (1 + beta * R1[z]) - dt * (
                            d + mu) * G1[z]
                C2[z] = (1 - nu) * C1[z] + nu * C1[z - 1] + dt * d * e * G1[z] * (1 - Csda) - dt * rmr * m * C1[z] - dt * mu * C1[z]
            else:
                G2[z] = (1 + nu) * G1[z] - nu * G1[z + 1] + dt * alphaZ[z] * R1[z] * C1[z] / (1 + beta * R1[z]) - dt * (
                            d + mu) * G1[z]
                C2[z] = (1 + nu) * C1[z] - nu * C1[z + 1] + dt * d * e * G1[z] * (1 - Csda) - dt * rmr * m * C1[z] - dt * mu * C1[z]

        if w >= 0:
            rmr = 0
            G2[0] = (1 - nu) * G1[0] + dt * alphaZ[z] * R1[0] * C1[0] / (1 + beta * R1[0]) - dt * (d + mu) * G1[0]
            C2[0] = (1 - nu) * C1[0] + dt * e * d * G1[0] * (1 - Csda) - dt * rmr * m * C1[0] - dt * mu * C1[0]

            G2[-1] = G1[-1] + nu * G1[-2] + dt * alphaZ[z] * R1[-1] * C1[-1] / (1 + beta * R1[-1]) - dt * (d + mu) * G1[-1]
            C2[-1] = C1[-1] + nu * C1[-2] + dt * e * d * G1[-1] * (1 - Csda) - dt * rmr_last * m * C1[-1] - dt * mu * C1[-1]
        else:
            rmr = 0
            G2[0] = G1[0] - nu * G1[1] + dt * alphaZ[z] * R1[0] * C1[0] / (1 + beta * R1[0]) - dt * (d + mu) * G1[0]
            C2[0] = C1[0] - nu * C1[1] + dt * e * d * G1[0] * (1 - Csda) - dt * rmr * m * C1[0] - dt * mu * C1[0]

            G2[-1] = (1 + nu) * G1[-1] + dt * alphaZ[z] * R1[-1] * C1[-1] / (1 + beta * R1[-1]) - dt * (d + mu) * G1[-1]
            C2[-1] = (1 + nu) * C1[-1] + dt * e * d * G1[-1] * (1 - Csda) - dt * rmr_last * m * C1[-1] - dt * mu * C1[-1]

        R1 = np.copy(R2)

        # Keep the same gaussian distribution as set at the beginning
        if np.round(t % 24, 1) == 0:
            C1 = np.zeros(zz)
            G1 = np.zeros(zz)
            C1[0:500] = gaussian(watercolumn[0:500], np.max(C2), cen, sigma)
            G1[0:500] = gaussian(watercolumn[0:500], np.max(G2), cen, sigma)
        else:
            C1 = np.copy(C2)
            G1 = np.copy(G2)

        # Data are stored every dt_save timestep
        if t <= temps[nh] < round(t + dt, 7):
            nh = nh + 1
            R[nh] = np.copy(R2)
            G[nh] = np.copy(G2)
            C[nh] = np.copy(C2)

            Dm[nh] = np.copy(D1m)
            Dg[nh] = np.copy(D1g)

            D1m = np.zeros(zz)
            D1g = np.zeros(zz)
        t = round(t + dt, 7)
    return R, G, C, Dm, Dg

# The numerical scheme for the model for dvm_model_seasonal_v2, variation of carrying capacity and temperature over one year
@jit(nopython=True)
def AN_RGCD_T2(R0, C0, V, annual_Tk, d, coef_alpha, beta, I_norm, dt_inf, Kt, rhoT):
    # Initialisation
    R = np.zeros((tt, zz))
    C = np.zeros((tt, zz))
    G = np.zeros((tt, zz))

    Dm = np.zeros((tt, zz))
    Dg = np.zeros((tt, zz))

    R[0] = R1 = np.copy(R0)
    C[0] = C1 = np.copy(C0)

    G1 = np.zeros(zz)
    D1g = np.zeros(zz)
    D1m = np.zeros(zz)

    t = 0
    nt = 0
    nbhy = 334*24

    ## Loop over time: Euler explicit
    while t < tmax - 2 * dt_save:
        R2 = np.zeros(zz)
        C2 = np.zeros(zz)
        G2 = np.zeros(zz)
        w = V[int((t % 24) / dt_inf)]
        if w == 0:
            dt = dt_sup
        else:
            dt = min(dt_sup, dz / np.abs(w) * CFL)
        nu = w * dt / dz

        alphaZ = I_norm[:, int((t % 24) / dt_inf)] * coef_alpha

        # Select T and rho for one timestep along depth
        TkL = annual_Tk[:, int(t % nbhy)]
        rhoL = rhoT[:, int(t % nbhy)]

        rmr_last = RMR_func(Wc_ind, TkL[-1], zz * dz, a0, a1, a2, a3, RQ)
        ## Loop over space: upwind
        for z in range(zz):
            depth = z * dz
            # Maintenance rate : RMR + AMR + excretion
            rmr = RMR_func(Wc_ind, TkL[z], depth, a0, a1, a2, a3, RQ)  # h-1
            # Remineralization rate of FP:
            r = remineratization_rateF(TkL[z] - 273.15)

            if bm_id == "bm_crust":
                m = (1 + 0.31 + (np.abs(w) / w_ref))
            else:
                m = (1 + (np.abs(w) / w_ref))

            sda = 0.175

            D1m[z] += dt * m * rmr * (C1[z]+G1[z]) + dt * e * d * G1[z] * sda
            D1g[z] += dt * (1 - e) * d * G1[z] - r*D1g[z-1]*dt

            if Kt[z] != 0:
                R2[z] = R1[z] + dt * (rhoL[z] * R1[z] * (1 - R1[z] / Kt[z])) - dt * alphaZ[z] * R1[z] * C1[z] / (1 + beta * R1[z])

        for z in range(1, zz - 1):
            rmr = RMR_func(Wc_ind, TkL[z], z * dz, a0, a1, a2, a3, RQ)
            if w >= 0:
                G2[z] = (1 - nu) * G1[z] + nu * G1[z - 1] + dt * alphaZ[z] * R1[z] * C1[z] / (1 + beta * R1[z]) - dt * (
                            d + mu) * G1[z]
                C2[z] = (1 - nu) * C1[z] + nu * C1[z - 1] + dt * d * e * G1[z] * (1 - sda) - dt * rmr * m * C1[z] - dt * mu * C1[z]
            else:
                G2[z] = (1 + nu) * G1[z] - nu * G1[z + 1] + dt * alphaZ[z] * R1[z] * C1[z] / (1 + beta * R1[z]) - dt * (
                            d + mu) * G1[z]
                C2[z] = (1 + nu) * C1[z] - nu * C1[z + 1] + dt * d * e * G1[z] * (1 - sda) - dt * rmr * m * C1[z] - dt * mu * C1[z]

        if w >= 0:
            rmr = 0
            G2[0] = (1 - nu) * G1[0] + dt * alphaZ[z] * R1[0] * C1[0] / (1 + beta * R1[0]) - dt * (d + mu) * G1[0]
            C2[0] = (1 - nu) * C1[0] + dt * e * d * G1[0] * (1 - sda) - dt * rmr * m * C1[0] - dt * mu * C1[0]

            G2[-1] = G1[-1] + nu * G1[-2] + dt * alphaZ[z] * R1[-1] * C1[-1] / (1 + beta * R1[-1]) - dt * (d + mu) * G1[-1]
            C2[-1] = C1[-1] + nu * C1[-2] + dt * e * d * G1[-1] * (1 - sda) - dt * rmr_last * m * C1[-1] - dt * mu * C1[-1]
        else:
            rmr = 0
            G2[0] = G1[0] - nu * G1[1] + dt * alphaZ[z] * R1[0] * C1[0] / (1 + beta * R1[0]) - dt * (d + mu) * G1[0]
            C2[0] = C1[0] - nu * C1[1] + dt * e * d * G1[0] * (1 - sda) - dt * rmr * m * C1[0] - dt * mu * C1[0]

            G2[-1] = (1 + nu) * G1[-1] + dt * alphaZ[z] * R1[-1] * C1[-1] / (1 + beta * R1[-1]) - dt * (d + mu) * G1[-1]
            C2[-1] = (1 + nu) * C1[-1] + dt * e * d * G1[-1] * (1 - sda) - dt * rmr_last * m * C1[-1] - dt * mu * C1[-1]

        R1 = np.copy(R2)
        # Keep the same gaussian distribution as set at the beginning
        if np.round(t % 24, 1) == 0:
            C1 = np.zeros(zz)
            G1 = np.zeros(zz)
            C1[0:500] = gaussian(watercolumn[0:500], np.max(C2), cen, sigma)
            G1[0:500] = gaussian(watercolumn[0:500], np.max(G2), cen, sigma)
        else:
            C1 = np.copy(C2)
            G1 = np.copy(G2)


        # Data are stored every dt_save timestep
        if t <= temps[nt] < round(t + dt, 7):
            nt = nt + 1
            R[nt] = np.copy(R2)
            G[nt] = np.copy(G2)
            C[nt] = np.copy(C2)

            Dm[nt] = np.copy(D1m)
            Dg[nt] = np.copy(D1g)

            D1m = np.zeros(zz)
            D1g = np.zeros(zz)
        t = round(t + dt, 7)

    return R, G, C, Dg, Dm

# The numerical scheme for the model for dvm_model_seasonal_v3 : variation of surface irradiance
@jit(nopython=True)
def AN_RGCD_I2(R0, C0, VL, TkL, Kd_annual, I_norm, d, coef_alpha, beta, dt_inf, KL):
    # Initialisation
    R = np.zeros((tt, zz))
    C = np.zeros((tt, zz))
    G = np.zeros((tt, zz))
    Dm = np.zeros((tt, zz))
    Dg = np.zeros((tt, zz))
    R[0] = R1 = np.copy(R0)
    C[0] = C1 = np.copy(C0)
    G1 = np.zeros(zz)
    D1g = np.zeros(zz)
    D1m = np.zeros(zz)

    t = 0
    nt = 0

    alphaL = []

    ## Loop over temps: Euler explicit
    while t < tmax - 2 * dt_save:

        R2 = np.zeros(zz)
        C2 = np.zeros(zz)
        G2 = np.zeros(zz)

        day = int(t/24)

        # Define speed, surface irradiance and Kd490
        w = VL[day, int((t % 24) / dt_inf)]
        I = I_norm[day, int((t % 24) / dt_inf)]
        #Kd = Kd_annual[day]

        if w == 0:
            dt = dt_sup
        else:
            dt = min(dt_sup, dz / np.abs(w) * CFL)
        nu = w * dt / dz

        rmr_last = RMR_func(Wc_ind, TkL[-1], zz * dz, a0, a1, a2, a3, RQ)
        ## Loop over space: upwind
        alphaZ = np.zeros(zz)
        for z in range(zz):
            depth = z * dz
            alphaZ[z] = beer_lambert(I, depth) * coef_alpha

            if depth == cen:
                if (np.round(t % 24, 1)) == 0:
                    alphaL.append(beer_lambert(I, depth) * coef_alpha)

            rmr = RMR_func(Wc_ind, TkL[z], depth, a0, a1, a2, a3, RQ)  # in [R] h-1
            sda = 0.175

            # Maintenance rate : RMR + AMR + excretion
            if bm_id == "bm_crust":
                m = (1 + 0.31 + (np.abs(w) / w_ref))
            else:
                m = (1 + (np.abs(w) / w_ref))

            D1m[z] += dt * rmr * (C1[z]+G1[z]) * m + dt * e * d * G1[z] * sda
            D1g[z] += dt * d * (1 - e) * G1[z]

            if KL[z] != 0:
                R2[z] = R1[z] + dt * (rho * R1[z] * (1 - R1[z] / KL[z])) - dt * alphaZ[z] * R1[z] * C1[z] / (1 + beta * R1[z])

        for z in range(1, zz - 1):
            rmr = RMR_func(Wc_ind, TkL[z], z * dz, a0, a1, a2, a3, RQ)
            if w >= 0:
                G2[z] = (1 - nu) * G1[z] + nu * G1[z - 1] + dt * alphaZ[z] * R1[z] * C1[z] / (1 + beta * R1[z]) - dt * (
                            d + mu) * G1[z]
                C2[z] = (1 - nu) * C1[z] + nu * C1[z - 1] + dt * d * e * G1[z] * (1 - sda) - dt * rmr * m * C1[z] - dt * mu * C1[z]
            else:
                G2[z] = (1 + nu) * G1[z] - nu * G1[z + 1] + dt * alphaZ[z] * R1[z] * C1[z] / (1 + beta * R1[z]) - dt * (
                            d + mu) * G1[z]
                C2[z] = (1 + nu) * C1[z] - nu * C1[z + 1] + dt * d * e * G1[z] * (1 - sda) - dt * rmr * m * C1[z] - dt * mu * C1[z]

        if w >= 0:
            rmr = 0
            G2[0] = (1 - nu) * G1[0] + dt * alphaZ[z] * R1[0] * C1[0] / (1 + beta * R1[0]) - dt * (d + mu) * G1[0]
            C2[0] = (1 - nu) * C1[0] + dt * e * d * G1[0] * (1 - sda) - dt * rmr * m * C1[0] - dt * mu * C1[0]

            G2[-1] = G1[-1] + nu * G1[-2] + dt * alphaZ[z] * R1[-1] * C1[-1] / (1 + beta * R1[-1]) - dt * (d + mu) * G1[-1]
            C2[-1] = C1[-1] + nu * C1[-2] + dt * e * d * G1[-1] * (1 - sda) - dt * rmr_last * m * C1[-1] - dt * mu * C1[-1]
        else:
            rmr = 0
            G2[0] = G1[0] - nu * G1[1] + dt * alphaZ[z] * R1[0] * C1[0] / (1 + beta * R1[0]) - dt * (d + mu) * G1[0]
            C2[0] = C1[0] - nu * C1[1] + dt * e * d * G1[0] * (1 - sda) - dt * rmr * m * C1[0] - dt * mu * C1[0]

            G2[-1] = (1 + nu) * G1[-1] + dt * alphaZ[z] * R1[-1] * C1[-1] / (1 + beta * R1[-1]) - dt * (d + mu) * G1[-1]
            C2[-1] = (1 + nu) * C1[-1] + dt * e * d * G1[-1] * (1 - sda) - dt * rmr_last * m * C1[-1] - dt * mu * C1[-1]

        R1 = np.copy(R2)
        # Keep the same gaussian distribution as set at the beginning
        if np.round(t % 24, 1) == 0:
            C1 = np.zeros(zz)
            G1 = np.zeros(zz)
            C1[0:500] = gaussian(watercolumn[0:500], np.max(C2), cen, sigma)
            G1[0:500] = gaussian(watercolumn[0:500], np.max(G2), cen, sigma)
        else:
            C1 = np.copy(C2)
            G1 = np.copy(G2)

        # Data are stored every dt_save timestep
        if t <= temps[nt] < round(t + dt, 10):
            nt = nt + 1
            R[nt] = np.copy(R2)
            G[nt] = np.copy(G2)
            C[nt] = np.copy(C2)

            Dm[nt] = np.copy(D1m)
            Dg[nt] = np.copy(D1g)

            D1m = np.zeros(zz)
            D1g = np.zeros(zz)
        t = round(t + dt, 10)
    return R, G, C, Dg, Dm, alphaL

# The numerical scheme for the model for dvm_model_seasonal_v : variation of the 3 environmental parameters
@jit(nopython=True)
def AN_RGCD_KTI2(R0, C0, VL, annual_Tk, Kd_annual, I_norm, d, coef_alpha, beta, dt_inf, rhoT, NPP):
    # Initialisation
    R = np.zeros((tt, zz))
    C = np.zeros((tt, zz))
    G = np.zeros((tt, zz))

    Dm = np.zeros((tt, zz))
    Dg = np.zeros((tt, zz))

    R[0] = R1 = np.copy(R0)
    C[0] = C1 = np.copy(C0)

    G1 = np.zeros(zz)
    D1g = np.zeros(zz)
    D1m = np.zeros(zz)

    alphaL = []

    nbhy = 334 * 24
    t = 0
    nt = 0

    ## Loop over time: Euler explicit
    while t < tmax - 2 * dt_save:
        R2 = np.zeros(zz)
        C2 = np.zeros(zz)
        G2 = np.zeros(zz)

        day = int(t / 24)

        # Define speed, surface irradiance and Kd490
        w = VL[day, int((t % 24) / dt_inf)]
        I = I_norm[day, int((t % 24) / dt_inf)]
        Kd = Kd_annual[day]

        if w == 0:
            dt = dt_sup
        else:
            dt = min(dt_sup, dz / np.abs(w) * CFL)
        nu = w * dt / dz

        # Select T for one timestep along depth
        TkL = annual_Tk[:, int(t % nbhy)]
        rhoL = rhoT[:, int(t % nbhy)]
        Kt = NPP[:, int(t % nbhy)]

        rmr_last = RMR_func(Wc_ind, TkL[-1], zz * dz, a0, a1, a2, a3, RQ)
        ## Loop over space: upwind
        alphaZ = np.zeros(zz)
        for z in range(zz):
            depth = z * dz
            alphaZ[z] = beer_lambert(I, depth, k=Kd) * coef_alpha

            if depth == cen:
                if (np.round(t % 24, 1)) == 0:
                    alphaL.append(beer_lambert(I, depth, k=Kd) * coef_alpha)

            rmr = RMR_func(Wc_ind, TkL[z], depth, a0, a1, a2, a3, RQ)  # in [R] h-1
            sda = 0.175

            # Remineralization rate of FP:
            r = remineratization_rateF(TkL[z] - 273.15)

            # Maintenance rate : RMR + AMR + excretion
            m = (1 + (np.abs(w) / w_ref))

            D1m[z] += dt * rmr * C1[z] * m + dt * e * d * G1[z] * sda
            D1g[z] += dt * (1 - e) * d * G1[z] - r * D1g[z-1] * dt

            if Kt[z] != 0:
                R2[z] = R1[z] + dt * (rhoL[z] * R1[z] * (1 - R1[z] / Kt[z])) - dt * alphaZ[z] * R1[z] * C1[z] / (1 + beta * R1[z])

        for z in range(1, zz - 1):
            rmr = RMR_func(Wc_ind, TkL[z], z * dz, a0, a1, a2, a3, RQ)
            if w >= 0:
                G2[z] = (1 - nu) * G1[z] + nu * G1[z - 1] + dt * alphaZ[z] * R1[z] * C1[z] / (1 + beta * R1[z]) - dt * (
                            d + mu) * G1[z]
                C2[z] = (1 - nu) * C1[z] + nu * C1[z - 1] + dt * d * e * G1[z] * (1 - sda) - dt * rmr * m * C1[z] - dt * mu * C1[z]
            else:
                G2[z] = (1 + nu) * G1[z] - nu * G1[z + 1] + dt * alphaZ[z] * R1[z] * C1[z] / (1 + beta * R1[z]) - dt * (
                            d + mu) * G1[z]
                C2[z] = (1 + nu) * C1[z] - nu * C1[z + 1] + dt * d * e * G1[z] * (1 - sda) - dt * rmr * m * C1[z] - dt * mu * C1[z]

        if w >= 0:
            rmr = 0
            G2[0] = (1 - nu) * G1[0] + dt * alphaZ[z] * R1[0] * C1[0] / (1 + beta * R1[0]) - dt * (d + mu) * G1[0]
            C2[0] = (1 - nu) * C1[0] + dt * e * d * G1[0] * (1 - sda) - dt * rmr * m * C1[0] - dt * mu * C1[0]

            G2[-1] = G1[-1] + nu * G1[-2] + dt * alphaZ[z] * R1[-1] * C1[-1] / (1 + beta * R1[-1]) - dt * (d + mu) * G1[-1]
            C2[-1] = C1[-1] + nu * C1[-2] + dt * e * d * G1[-1] * (1 - sda) - dt * rmr_last * m * C1[-1] - dt * mu * C1[-1]
        else:
            rmr = 0
            G2[0] = G1[0] - nu * G1[1] + dt * alphaZ[z] * R1[0] * C1[0] / (1 + beta * R1[0]) - dt * (d + mu) * G1[0]
            C2[0] = C1[0] - nu * C1[1] + dt * e * d * G1[0] * (1 - sda) - dt * rmr * m * C1[0] - dt * mu * C1[0]

            G2[-1] = (1 + nu) * G1[-1] + dt * alphaZ[z] * R1[-1] * C1[-1] / (1 + beta * R1[-1]) - dt * (d + mu) * G1[-1]
            C2[-1] = (1 + nu) * C1[-1] + dt * e * d * G1[-1] * (1 - sda) - dt * rmr_last * m * C1[-1] - dt * mu * C1[-1]

        R1 = np.copy(R2)
        # Keep the same gaussian distribution as set at the beginning
        if np.round(t % 24, 1) == 0:
            C1 = np.zeros(zz)
            G1 = np.zeros(zz)
            C1[0:500] = gaussian(watercolumn[0:500], np.max(C2), cen, sigma)
            G1[0:500] = gaussian(watercolumn[0:500], np.max(G2), cen, sigma)
        else:
            C1 = np.copy(C2)
            G1 = np.copy(G2)

        # Data are stored every dt_save timestep
        if t <= temps[nt] < round(t + dt, 10):
            nt = nt + 1
            R[nt] = np.copy(R2)
            G[nt] = np.copy(G2)
            C[nt] = np.copy(C2)

            Dm[nt] = np.copy(D1m)
            Dg[nt] = np.copy(D1g)

            D1m = np.zeros(zz)
            D1g = np.zeros(zz)
        t = round(t + dt, 10)
    return R, G, C, Dg, Dm, alphaL
