"""
Plot the results of the model from dvm_model_vp.py
"""
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from pylab import meshgrid, imshow, colorbar
from parameters import *
from plotting_funcs import nonlinear_colormap, contour_levels_func, get_continuous_cmap


import cmocean
cm = nonlinear_colormap()

path_to_results = Path('~/simul_dvm2/Results').expanduser()
path_to_plot = Path('~/simul_dvm2/Plots').expanduser()
path_to_data = Path('~/simul_dvm2/Data').expanduser()

# Load results
#path_zoo = path_to_data / 'zoo_annual.txt'
path_R = path_to_results / 'results_R.txt'
path_G = path_to_results / 'results_G.txt'
path_C = path_to_results / 'results_C.txt'
path_RMR = path_to_results / 'results_M.txt'
path_Dg = path_to_results / 'results_Dg.txt'
path_SDA = path_to_results / 'results_SDA.txt'
path_AMR = path_to_results / 'results_AMR.txt'

#zoo = np.loadtxt(path_zoo)
R = np.loadtxt(path_R)
Req = R[-1, :]
C = np.loadtxt(path_C)
Ceq = C[-1, :]
G = np.loadtxt(path_G)
print("Req=", np.max(Req))
print("Ceq=", np.max(Ceq))
print("C/R=", np.max(Ceq)/np.max(Req))

if bm_id == "bm_fish":
    path_Ceq = path_to_results / 'Ceq_fish.txt'
    path_Req = path_to_results / 'Req_fish.txt'
elif bm_id == "bm_crust":
    path_Ceq = path_to_results / 'Ceq_crust.txt'
    path_Req = path_to_results / 'Req_crust.txt'
else:
    path_Ceq = path_to_results / 'Ceq_ceph.txt'
    path_Req = path_to_results / 'Req_ceph.txt'
np.savetxt(path_Req, Req)
np.savetxt(path_Ceq, Ceq)

# Integration of the state variables and detritus along depth
Ctot = C.sum(axis=1)*dz
Rtot = R.sum(axis=1)*dz
Gtot = G.sum(axis=1)*dz

Ctot_end = Ctot[-1]
Rtot_end = Rtot[-1]

# Last 24h values of R and C of the simulation
plt.ion()
fig = plt.figure()
for i in range(-int((24/dt_save)), -1, 2):
    plt.title("R and C at equilibrium")
    plt.plot(R[i, :], -watercolumn, label='Resource')
    plt.plot(C[i, :], -watercolumn, label='Consumers')
    plt.legend()

    # to flush the GUI events
    plt.pause(0.001)
    plot_name_RCeq = path_to_plot / "RCeq.pdf"
    plt.savefig(plot_name_RCeq)
    plt.clf()
plt.close('all')

# Plot of the integrated biomass over time of the 3 state variables
month = temps[:len(Ctot)]/24/30.5
fig, ax = plt.subplots()
#ax.plot(month, Rtot, label='Resource')
ax.plot(month, Gtot, label='Gut')
#ax.plot(month, Ctot, label='Consumer')
#ax.plot(month, zoo, label='Zooplankton')
ax.set_xlabel('Month ($h$)')
ax.set_ylabel('Carbon mass ($mgC$ $m^{-2}$)')
textstr = '\n'.join((
    r'$\mathrm{Vmax (m/h)}=%.1f$' % (int(vmax), ),
    r'$\mathrm{µ (h^{-1})}=%.4f$' % (mu, ),
    r'$\mathrm{\alpha (h^{-1})}=%.1f$' % (coef_alpha, ),
    r'$\mathrm{\beta (h R^{-1})}=%.2f$' % (beta, ),
    r'$\mathrm{e}=%.2f$' % (e, )))
# these are matplotlib.patch.Patch properties
props = dict(boxstyle='round', edgecolor='none', facecolor='blue', alpha=0.05)
# place a text box in upper left in axes coords
ax.text(0.05, 0.55, textstr, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', bbox=props)
plt.legend()
name_RGCtot = 'RGCtot.pdf'
plot_RGCtot = path_to_plot / name_RGCtot
plt.savefig(plot_RGCtot)

# Pattern of consumers migration for each day
colors = ["#1f0954", "#1e488f", "#7bc8f6", "#6fc276", "#fcb001", "#fff917"]
cmap_dvm = get_continuous_cmap(colors)
dt = dt_save

fig, ax = plt.subplots(1, 1, figsize=(5, 3), layout="constrained")
props = dict(boxstyle='round', edgecolor='black', facecolor='white', alpha=0.9)
ims = []
tmax = len(Ctot)*dt
for i in range(int(tmax/24)):
    im = ax.imshow(np.transpose(C[int(i*24/dt):int((i*24/dt)+24/dt)]), aspect='auto', cmap=cmap_dvm, extent=[-12, 12, Zmax, 0], animated=True)
    t = ax.text(-10, 150, "Day " + str(i), bbox=props)
    plt.xticks([-12, -6, 0, 6, 12])
    plt.xlabel('Time relative to local noon [$h$]')
    plt.ylabel('Depth [$m$]')
    """plt.pause(0.01)
    plot_dvmC = path_to_plot / 'colormapDVMC.pdf'
    plt.savefig(plot_dvmC)
    plt.clf()"""
    if i == 0:
        #cbar = fig.colorbar(im, boundaries=np.linspace(0,0.05,6))
        cbar = fig.colorbar(im)
        cbar.set_label("$mgC$ $m^{-2}$")
    ims.append([im, t])

ani = animation.ArtistAnimation(fig, ims, interval=100, blit=True, repeat=False)
plt.show()
plot_dvmC = path_to_plot / 'colormapDVMC.gif'
ani.save(plot_dvmC, writer="pillow")

plt.close('all')

# Last day resource concentration
fig = plt.figure(figsize=(5, 3), layout="constrained")
y = np.copy(temps)
x = np.copy(watercolumn)
X, Y = meshgrid(x, y)  # grid of point
im = imshow(np.transpose(C[-int(24/dt):-1]), aspect='auto', extent=[-12, 12, Zmax, 0], cmap=cmap_dvm)  # drawing the function # cmap=cm.RdBu
plt.xticks([-12, -6, 0, 6, 12])
plt.xlabel('Time relative to local noon [$h$]')
cbar = fig.colorbar(im)
cbar.set_label("$mgC$ $m^{-2}$")  # adding the colobar on the right
plt.ylabel('Depth [$m$]')
name_dvmR = 'colormapDVMC.pdf'
plot_dvmR = path_to_plot / name_dvmR
plt.savefig(plot_dvmR)
plt.close('all')

# Last day resource concentration
fig = plt.figure(figsize=(6, 5))
y = np.copy(temps)
x = np.copy(watercolumn)
X, Y = meshgrid(x, y)  # grid of point
im = imshow(np.transpose(R[-int(24/dt):-1]), aspect='auto', extent=[0, 24, Zmax, 0])  # drawing the function # cmap=cm.RdBu,
colorbar(im)  # adding the colobar on the right
plt.xlabel('Time [$h$]')
plt.ylabel('Depth [$m$]')
name_dvmR = 'colormapDVMR.pdf'
plot_dvmR = path_to_plot / name_dvmR
plt.savefig(plot_dvmR)
plt.close('all')

### DETRITUS ###
Dg = np.loadtxt(path_Dg)
Drmr = np.loadtxt(path_RMR)
Dsda = np.loadtxt(path_SDA)
Damr = np.loadtxt(path_AMR)
Dresp = Dsda + Drmr + Damr

#Dead bodies of consumers and their gut
Du = (C+G) * mu*dt

# Number of days
nj = tmax/24

det_names = {"Dresp": Dresp, "Dg": Dg, "Du": Du}
plt.figure(figsize=(7, 5))

tt = len(Ctot)
temps = np.arange(0, tt*dt_save, dt_save)
# Production of detritus over time integrated along the depth
for Dname, nd in det_names.items():
    D = nd
    # Production of detritus over time integrated along the depth
    Dtot = np.zeros(tt)
    for t in range(tt):
        Dtot[t] = np.sum(D[t])
    Dtot *= dz

    fig = plt.figure()
    plt.title('Detritus production')
    plt.plot(temps, Dtot)
    plt.xlabel('Time (h)')
    plt.ylabel('Carbon mass (arbitrary unit)')

    name_det = Dname + '_tot.pdf'
    plot_name_det = path_to_plot / name_det
    plt.savefig(plot_name_det)
    plt.close()

# Integration along depth
Dg_tot = Dg.sum(axis=0)*dz
Maint_tot = Dresp.sum(axis=0)*dz
Mort_tot = Du.sum(axis=0)*dz
Drmr_tot = Drmr.sum(axis=0)*dz
Dsda_tot = Dsda.sum(axis=0)*dz
Damr_tot = Damr.sum(axis=0)*dz

# Integration over time
#FP
Dg_1d = np.reshape(Dg, -1)
arr3D = np.reshape(Dg_1d, (int(nj), int((24/dt_save)), zz))
Dg_time = np.sum(arr3D, axis=0).sum(1)
Dg_time *= dt
#Dead bodies
Du_1d = np.reshape(Du, -1)
arr3D = np.reshape(Du_1d, (int(nj), int((24/dt_save)), zz))
Du_time = np.sum(arr3D, axis=0).sum(1)
Du_time *= dt
#RMR
Drmr_1d = np.reshape(Drmr, -1)
arr3D = np.reshape(Drmr_1d, (int(nj), int((24/dt_save)), zz))
Drmr_time = np.sum(arr3D, axis=0).sum(1)
Drmr_time *= dt
#SDA
Dsda_1d = np.reshape(Dsda, -1)
arr3D = np.reshape(Dsda_1d, (int(nj), int((24/dt_save)), zz))
Dsda_time = np.sum(arr3D, axis=0).sum(1)
Dsda_time *= dt
#AMR
Damr_1d = np.reshape(Damr, -1)
arr3D = np.reshape(Damr_1d, (int(nj), int((24/dt_save)), zz))
Damr_time = np.sum(arr3D, axis=0).sum(1)
Damr_time *= dt

#Total Respiration
Maint_time = Drmr_time + Damr_time + Dsda_time

### Plot of the integrated biomass along depth of the detritus
det_dict = {'Respiration': Maint_tot.ravel().tolist(),
            'Fecal pellets': Dg_tot.ravel().tolist(),
            #'Excretion': Exc_tot.ravel().tolist(),
            'Dead bodies': Mort_tot.ravel().tolist()}
# Plot of the integrated biomass over time of the detritus
det_dict_time = {'Respiration': Maint_time.ravel().tolist(),
            'Fecal pellets': Dg_time.ravel().tolist(),
            #'Excretion': Exc_tot.ravel().tolist(),
            'Dead bodies': Du_time.ravel().tolist()}

fig, ax = plt.subplots()
color_map = ["papayawhip", "rosybrown", "dimgrey"]
ax.stackplot(watercolumn, det_dict.values(),
             labels=det_dict.keys(), alpha=0.8, colors=color_map)
plt.ylabel('Detritus concentrations [$mgC$ $m^{-2}$ $h^{-1}$]')
plt.xlabel('Depth [$m$]')
plt.legend()
fig.tight_layout()
name_det_tot = 'Det_tot.pdf'
plot_det_tot = path_to_plot / name_det_tot
plt.savefig(plot_det_tot)
plt.close()

#All detritus integrated over time
xtemps = np.arange(0, 24, dt_save)
fig, ax = plt.subplots()

ax.plot(xtemps, det_dict_time['Fecal pellets'], label='Fecal pellets')
ax.plot(xtemps, det_dict_time['Dead bodies'], label='Dead bodies')
ax.plot(xtemps, det_dict_time['Respiration'], label='Respiration')
plt.ylabel('Detritus concentrations [$mgC$ $m^{-2}$ $h^{-1}$]')
plt.xlabel('Time [$h$]')
plt.legend()
fig.tight_layout()
name_det_tot = 'Det_tot_time.pdf'
plot_det_tot = path_to_plot / name_det_tot
plt.savefig(plot_det_tot)

# Plot of the respiration integrated over time (respD_dict) and along depth (respT_dict)
respD_dict = pd.DataFrame({'RMR': Drmr_tot.ravel().tolist(),
            'AMR': Damr_tot.ravel().tolist(),
            'SDA': Dsda_tot.ravel().tolist(), })

time_last24 = temps[-int(24/dt):-1]
respT_dict = pd.DataFrame({
                           'RMR': Drmr_time.ravel().tolist(),
                           'AMR': Damr_time.ravel().tolist(),
                           'SDA': Dsda_time.ravel().tolist(), })

respD_perc = respD_dict.divide(respD_dict.sum(axis=1), axis=0)*100
respT_perc = respT_dict.divide(respT_dict.sum(axis=1), axis=0)*100

## Proportions of respiration pathways
fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, layout="constrained", figsize=(6,3))
color_map = ["papayawhip", "lightskyblue", "indianred"]
ax2.stackplot(watercolumn,  respD_perc["RMR"],  respD_perc["AMR"],  respD_perc["SDA"], labels=['RMR', "AMR", "SDA"], alpha=0.8, colors=color_map)
ax1.stackplot(xtemps,  respT_perc["RMR"],  respT_perc["AMR"],  respT_perc["SDA"], labels=['RMR', "AMR", "SDA"], alpha=0.8, colors=color_map)
#ax.stackplot(watercolumn,  resp_dict.values(), labels=resp_dict.keys(), alpha=0.8, colors=color_map)
ax1.set_ylabel('Respiration [%]')
ax2.set_xlabel('Depth [$m$]')
ax1.set_xlabel('Time [$h$]')
plt.legend(loc='lower right')
ax1.set_xlim([0, 24])
ax2.set_xlim([30, Zmax])
ax1.set_ylim([0, 100])
ax2.set_ylim([0, 100])
ax1.set_xticks(np.arange(0, 25, 3), np.arange(0, 25, 3))
fig.set_size_inches(10, 4)
name_resp_tot = 'Resp_perc.pdf'
plot_resp_tot = path_to_plot / name_resp_tot
plt.savefig(plot_resp_tot)

### Plots of integrated respiration pathways over time and along depth
"""fig, (ax1, ax2) = plt.subplots(1, 2, layout="constrained")
color_map = ["moccasin", "lightskyblue", "indianred"]
ax2.stackplot(watercolumn,  respD_dict["RMR"],  respD_dict["AMR"],  respD_dict["SDA"], labels=['RMR', "AMR", "SDA"], alpha=0.8, colors=color_map)
ax1.stackplot(xtemps,  respT_dict["RMR"],  respT_dict["AMR"],  respT_dict["SDA"], labels=['RMR', "AMR", "SDA"], alpha=0.8, colors=color_map)
ax1.set_ylabel('Respiration [$mgC$ $m^{-2}$ $h^{-1}$]')
ax2.set_xlabel('Depth [$m$]')
ax1.set_xlabel('Time [$h$]')
plt.legend()
ax1.set_xticks(np.arange(0, 25, 3), np.arange(0, 25, 3))
fig.set_size_inches(10, 4)
name_resp_tot = 'Resp_tot.pdf'
plot_resp_tot = path_to_plot / name_resp_tot
plt.savefig(plot_resp_tot)"""