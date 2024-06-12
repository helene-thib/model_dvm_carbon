"""
Plot respiration
"""
import pandas as pd
from matplotlib import pyplot as plt
from sklearn import linear_model
import matplotlib as mpl
import seaborn as sns
from parameters import *
from dvm_functions import RMR_func, length_weight_func
import cmocean
path_to_plot = Path('~/simul_dvm2/Plots').expanduser()
from plotting_funcs import nonlinear_colormap, contour_levels_func, get_continuous_cmap
cmap = cmocean.cm.thermal

RMR_L = []
RMR_size = []
RMR_T = []
RMR_D = []

for size in sizeL:
    # Weight
    Wg, Wc_ind, Wd = length_weight_func(size, bm_id)  # Wg adn Wd in mg and Wc in mgC
    for Tk in [0, 5, 10, 15, 20, 25]:
        for depth in np.arange(5, 1000, 200):
            rmr = RMR_func(Wc_ind, Tk+273.15, depth, a0, a1, a2, a3, RQ)
            RMR_L.append(rmr)
            RMR_size.append(size)
            RMR_T.append(Tk)
            RMR_D.append(depth)

zipped = list(zip(RMR_L, RMR_size, RMR_T, RMR_D))
df = pd.DataFrame(zipped, columns=['RMR', 'Size', 'T [$°C$]', 'Depth [m]'])

f, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(10, 5), layout='constrained')
sns.lineplot(data=df, x="Size", y="RMR", hue="T [$°C$]", ax=ax1)
sns.lineplot(data=df, x="Size", y="RMR", hue="Depth [m]", ax=ax2)
ax1.set_xlabel('Size [$mm$]')
ax2.set_xlabel('Size [$mm$]')
ax1.set_ylabel('RMR [$mgC ind^{-1} h^{-1}$]')
#cbar = ax1.figure.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=df['T [$°C$]'].min(),vmax=df['T [$°C$]'].max(),clip=False), cmap=cmap),  label='T [$°C$]')

name_dvmR = 'Respi_size_temp.pdf'
plot_dvmR = path_to_plot / name_dvmR
plt.savefig(plot_dvmR)
plt.close()

RMR_taxo = []
RMR_weight = []
RMR_L = []

#Respiration rate from Ikeda 2014 and 2016
def RMR_02_func(CW, T, D, a0, a1, a2, a3, RQ): #CW in mgC, T in K, D for depth in m
    if D == 0:
        ln_RO = 0
    else:
        ln_RO = a0 + a1 * np.log(CW) + a2 * (1000/T) + a3 * np.log(D) #in µL O2 h-1
    R = np.exp(ln_RO) #in µgC h-1
    return R

for size in [20, 30, 40, 50, 60, 70, 80, 90, 100]:
    # Weight
    Wg, Wc_ind, Wd = length_weight_func(size, bm_id)  # Wg adn Wd in mg and Wc in mgC
    for bm_id in bm_idL:
        # Respiration coefficients for respi_func
        a0 = micron_dict[bm_id][8]
        a1 = micron_dict[bm_id][9]
        a2 = micron_dict[bm_id][10]
        a3 = micron_dict[bm_id][11]
        RQ = micron_dict[bm_id][12]

        rmr = RMR_02_func(Wc_ind, 295, 50, a0, a1, a2, a3, RQ)
        RMR_L.append(rmr)
        RMR_weight.append(Wg)
        RMR_taxo.append(bm_id)

zipped = list(zip(RMR_L, RMR_weight, RMR_taxo))
df2 = pd.DataFrame(zipped, columns=['RMR', 'Weight', 'Taxo'])
f, ax = plt.subplots(figsize=(10, 5), layout='constrained')
sns.lineplot(data=df2, x=np.log(df2["Weight"]), y=np.log(df2["RMR"]), hue="Taxo")
plt.legend()
ax.set_xlabel('Size [$mm$]')
ax.set_ylabel('RMR [$mgC ind^{-1} h^{-1}$]')
plt.show()
quit()

## Generating coef alpha coefficient for fish
lim_size = np.array([2, 8]).reshape(-1, 1)
lim_alpha_fish = [80, 130]
lm = linear_model.LinearRegression()
lm.fit(lim_size, lim_alpha_fish)
alphaL = np.arange(2, 10, 1).reshape(-1, 1)
coef_alpha_fish = lm.predict(alphaL)
a_fish = lm.intercept_
b_fish = lm.coef_
equation_fish = 'Fish :y='+str(np.round(a_fish, 1))+"x"+"+"+str(np.round(b_fish[0], 1))
print(coef_alpha_fish)
print(equation_fish)

## Generating coef alpha coefficient for crust
lim_size_crust = np.array([2, 5, 8]).reshape(-1, 1)
lim_alpha_crust = [32, 130, 180]
lm = linear_model.LinearRegression()
lm.fit(lim_size, lim_alpha_crust)
coef_alpha_crust = lm.predict(alphaL)
a_crust = lm.intercept_
b_crust = lm.coef_
equation_crust = 'Crustaceans :y='+str(np.round(a_crust, 1))+"x"+"+"+str(np.round(b_crust[0], 1))
print(coef_alpha_crust)

## Generating coef alpha coefficient for ceph
lim_alpha_ceph = [9, 40]
lm = linear_model.LinearRegression()
lm.fit(lim_size, lim_alpha_ceph)
coef_alpha_ceph = lm.predict(alphaL)
a_ceph = lm.intercept_
b_ceph = lm.coef_
equation_ceph = 'Squids :y='+str(np.round(a_ceph, 1))+"x"+"+"+str(np.round(b_ceph[0], 1))
print(coef_alpha_ceph)

plt.plot(alphaL, coef_alpha_fish, label=equation_fish)
plt.plot(alphaL, coef_alpha_crust, label=equation_crust)
plt.plot(alphaL, coef_alpha_ceph, label=equation_ceph)
plt.xlabel('Size [cm]')
plt.ylabel('$C_α$')
plt.legend()
name_fit = 'fit_alpha_size.pdf'
plot_fit = path_to_plot / name_fit
plt.savefig(plot_fit)



