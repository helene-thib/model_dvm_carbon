"""
Script of the numerical scheme for the model of carbon production by migrant micronekton with visual predation
"""

## IMPORTS
from numba import jit
from parameters import *
from dvm_functions import RMR_func, SDA_func, gaussian

from pathlib import Path
path_to_plot = Path('~/simul_dvm2/Plots').expanduser()


@jit(nopython=True)
def remineratization_rateF(Tc):
    r_day = 0.005*Tc + 0.011
    r = r_day/24
    return r


## The numerical scheme for the model
@jit(nopython=True)
def AN_RGC_VP2(R0, C0, Wc_ind, V, TkL, d, e, mu, w_ref, sda, beta, dt_inf, KL, I_norm, coef_alpha, RQ, a0, a1, a2, a3):
    # Initialisation
    R = np.zeros((tt, zz))
    C = np.zeros((tt, zz))
    G = np.zeros((tt, zz))

    R[0] = R1 = np.copy(R0)
    C[0] = C1 = np.copy(C0)

    RtotL = np.zeros(tt)
    CtotL = np.zeros(tt)

    G1 = np.zeros(zz)

    t = 0
    nt = 0
    rmr_last = RMR_func(Wc_ind, TkL[-1], Zmax, a0, a1, a2, a3, RQ)

    ### Loop over time
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

        ### Loop over depth, temperature, carrying capacity
        for z in range(zz):
            # Maintenance rate : RMR + AMR + excretion
            if bm_id == "bm_crust":
                m = (1 + 0.31 + (np.abs(w) / w_ref))
            else:
                m = (1 + (np.abs(w) / w_ref))

            if KL[z] != 0:
                R2[z] = R1[z] + dt * (rho * R1[z] * (1 - R1[z] / KL[z])) - dt * alphaZ[z] * R1[z] * C1[z] / (
                            1 + beta * R1[z])

        for z in range(1, zz - 1):
            rmr = RMR_func(Wc_ind, TkL[z], z*dz, a0, a1, a2, a3, RQ)
            if w >= 0:
                G2[z] = (1 - nu) * G1[z] + nu * G1[z - 1] + dt * alphaZ[z] * R1[z] * C1[z] / (1 + beta * R1[z]) - dt * (
                        d + mu) * G1[z]
                C2[z] = (1 - nu) * C1[z] + nu * C1[z - 1] + dt * d * e * G1[z] * (1 - sda) - dt * rmr * m * C1[
                    z] - dt * mu * C1[z]
            else:
                G2[z] = (1 + nu) * G1[z] - nu * G1[z + 1] + dt * alphaZ[z] * R1[z] * C1[z] / (1 + beta * R1[z]) - dt * (
                        d + mu) * G1[z]
                C2[z] = (1 + nu) * C1[z] - nu * C1[z + 1] + dt * d * e * G1[z] * (1 - sda) - dt * rmr * m * C1[
                    z] - dt * mu * C1[z]

        if w >= 0:
            rmr = 0
            G2[0] = (1 - nu) * G1[0] + dt * alphaZ[z] * R1[0] * C1[0] / (1 + beta * R1[0]) - dt * (d + mu) * G1[0]
            C2[0] = (1 - nu) * C1[0] + dt * e * d * G1[0] * (1 - sda) - dt * rmr * m * C1[0] - dt * mu * C1[0]

            G2[-1] = G1[-1] + nu * G1[-2] + dt * alphaZ[z] * R1[-1] * C1[-1] / (1 + beta * R1[-1]) - dt * (d + mu) * G1[
                -1]
            C2[-1] = C1[-1] + nu * C1[-2] + dt * e * d * G1[-1] * (1 - sda) - dt * rmr_last * m * C1[-1] - dt * mu * C1[
                -1]
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
            RtotL[nt] = np.sum(R2)
            CtotL[nt] = np.sum(C2)

            R[nt] = np.copy(R2)
            G[nt] = np.copy(G2)
            C[nt] = np.copy(C2)

            if nt > (24 * 61)/dt_save:
                last_index = nt - int((24*60) / dt_save)
                arrR = RtotL[last_index:nt]
                arrC = CtotL[last_index:nt]

                ratioR = np.std(arrR) / np.mean(arrR) * 100
                ratioC = np.std(arrC) / np.mean(arrC) * 100

                if ratioR < 1.5 and ratioC < 1.3 or CtotL[nt] < 0.001:
                    print("a l'equilibre", "ratioC=", ratioC, " et ratioR=", ratioR)
                    #last_day = nt - int(24 / dt_save)
                    return CtotL[nt], RtotL[nt], np.max(C2), np.max(R2)
                    #return R, G, C
                    #return R[nt], C[nt], np.max(R2), np.max(C2)
                    #return R[nt], C[nt], np.sum(R[nt]), np.sum(C[nt]) #for global descriptors

        t = round(t + dt, 7)
    print("pas a l'equilibre", "ratioC=", ratioC, " et ratioR=", ratioR)
    #last_day = nt - int(24 / dt_save)
    return CtotL[nt], RtotL[nt], np.max(C2), np.max(R2)
    #return R, G, C
    #return R[nt], C[nt]
    #return R[nt], C[nt], np.max(R2), np.max(C2)
    #return R[nt], C[nt], np.sum(R[nt]), np.sum(C[nt]) #for global descriptors


@jit(nopython=True)
## The numerical scheme for the model with the computation of the detritus
def AN_RGCD_VP(R0, C0, V, TkL, d, e, mu, w_ref, beta, dt_inf, KL, I_norm, coef_alpha, RQ, a0, a1, a2, a3, Wc_ind, sda):

    # initialisation
    R = np.zeros((tt, zz))
    C = np.zeros((tt, zz))
    G = np.zeros((tt, zz))

    Drmr = np.zeros((tt, zz))
    Damr = np.zeros((tt, zz))
    Dsda = np.zeros((tt, zz))
    Dg = np.zeros((tt, zz))

    R[0] = R1 = np.copy(R0)
    C[0] = C1 = np.copy(C0)

    G1 = np.zeros(zz)
    D1g = np.zeros(zz)
    D1rmr = np.zeros(zz)
    D1amr = np.zeros(zz)
    D1sda = np.zeros(zz)

    t = 0
    nt = 0
    rmr_last = RMR_func(Wc_ind, TkL[-1], zz * dz, a0, a1, a2, a3, RQ)

    # Loop over time
    while t < tmax - 2 * dt_save:
        R2 = np.zeros(zz)
        C2 = np.zeros(zz)
        G2 = np.zeros(zz)
        w = V[int((t % 24) / dt_inf)]
        if w == 0:
            dt = dt_sup
        else:
            dt = min(dt_sup, dz / np.abs(w) * CFL)
        # print('dt='+str(dt))
        nu = w * dt / dz

        alphaZ = I_norm[:, int((t % 24) / dt_inf)] * coef_alpha
        # Loop over space
        for z in range(zz):
            depth = z * dz
            # Maintenance rate : RMR + AMR + excretion
            rmr = RMR_func(Wc_ind, TkL[z], depth, a0, a1, a2, a3, RQ) # mgC ind-1 h-1

            if bm_id == "bm_crust":
                m = (1 + 0.31 + (np.abs(w) / w_ref))
            else:
                m = (1 + (np.abs(w) / w_ref))

            D1rmr[z] += dt * rmr * (C1[z] + G1[z])
            D1sda[z] += dt * d * e * G1[z] * sda
            D1amr[z] += dt * rmr * (C1[z] + G1[z]) * (np.abs(w) / w_ref)
            D1g[z] += dt * (1-e) * d * G1[z]

            if KL[z] != 0:
                R2[z] = R1[z] + dt * (rho * R1[z] * (1 - R1[z] / KL[z])) - dt * alphaZ[z] * R1[z] * C1[z] / (1 + beta * R1[z])

        for z in range(1, zz - 1):
            rmr = RMR_func(Wc_ind, TkL[z], z*dz, a0, a1, a2, a3, RQ)
            if w >= 0:
                G2[z] = (1 - nu) * G1[z] + nu * G1[z - 1] + dt * alphaZ[z] * R1[z] * C1[z] / (1 + beta * R1[z]) - dt * (
                        d + mu) * G1[z]
                C2[z] = (1 - nu) * C1[z] + nu * C1[z - 1] + dt * d * e * G1[z] * (1 - sda) - dt * rmr * m * C1[z] - dt*mu*C1[z]
            else:
                G2[z] = (1 + nu) * G1[z] - nu * G1[z + 1] + dt * alphaZ[z] * R1[z] * C1[z] / (1 + beta * R1[z]) - dt * (
                        d + mu) * G1[z]
                C2[z] = (1 + nu) * C1[z] - nu * C1[z + 1] + dt * d * e * G1[z] * (1 - sda) - dt * rmr * m * C1[z] - dt*mu*C1[z]

        if w >= 0:
            rmr = 0
            G2[0] = (1 - nu) * G1[0] + dt * alphaZ[z] * R1[0] * C1[0] / (1 + beta * R1[0]) - dt * (d + mu) * G1[0]
            C2[0] = (1 - nu) * C1[0] + dt * e * d * G1[0] * (1 - sda) - dt * rmr * m * C1[0] - dt*mu*C1[0]

            G2[-1] = G1[-1] + nu * G1[-2] + dt * alphaZ[z] * R1[-1] * C1[-1] / (1 + beta * R1[-1]) - dt * (d + mu) * G1[-1]
            C2[-1] = C1[-1] + nu * C1[-2] + dt * e * d * G1[-1] * (1 - sda) - dt * rmr_last * m * C1[-1] - dt*mu*C1[-1]
        else:
            rmr = 0
            G2[0] = G1[0] - nu * G1[1] + dt * alphaZ[z] * R1[0] * C1[0] / (1 + beta * R1[0]) - dt * (d + mu) * G1[0]
            C2[0] = C1[0] - nu * C1[1] + dt * e * d * G1[0] * (1 - sda) - dt * rmr * m * C1[0] - dt*mu*C1[0]

            G2[-1] = (1 + nu) * G1[-1] + dt * alphaZ[z] * R1[-1] * C1[-1] / (1 + beta * R1[-1]) - dt * (d + mu) * G1[-1]
            C2[-1] = (1 + nu) * C1[-1] + dt * e * d * G1[-1] * (1 - sda) - dt * rmr_last * m * C1[-1] - dt*mu*C1[-1]

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
            nt += 1
            R[nt] = np.copy(R2)
            G[nt] = np.copy(G2)
            C[nt] = np.copy(C2)

            Drmr[nt] = np.copy(D1rmr)
            Damr[nt] = np.copy(D1amr)
            Dsda[nt] = np.copy(D1sda)
            Dg[nt] = np.copy(D1g)

            D1rmr = np.zeros(zz)
            D1amr = np.zeros(zz)
            D1sda = np.zeros(zz)
            D1g = np.zeros(zz)

        t = round(t + dt, 7)
    #last_day = nt - int(24/ dt_save)
    #return np.sum(R[last_day:nt], axis=0), np.sum(G[last_day:nt], axis=0), np.sum(C[last_day:nt], axis=0), np.sum(Drmr[last_day:nt], axis=0), np.sum(Dg[last_day:nt], axis=0), np.sum(Dsda[last_day:nt], axis=0), np.sum(Damr[last_day:nt], axis=0)
    return R, G, C, Dg, Drmr, Dsda, Damr