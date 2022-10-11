import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from math import log10

##### read data from file

# arrays
with open('arrays_info.txt', 'r') as file:
    file.readline()
    tarr_info = file.readline().strip('\n').split()
    rarr_info = file.readline().strip('\n').split()
    aarr_info = file.readline().strip('\n').split()

# tarr 
tst, tend, Nt = float(tarr_info[2]), float(tarr_info[3]), int(tarr_info[4])
if tarr_info[1] == 'linear':
    tarr = np.linspace(tst, tend, Nt)
elif tarr_info[1] == 'logarithmic':
    tarr = np.logspace(log10(tst), log10(tend), Nt)
# print(tarr)

# rarr
rst, rend, Nr = float(rarr_info[2]), float(rarr_info[3]), int(rarr_info[4])
if rarr_info[1] == 'linear':
    rarr = np.linspace(rst, rend, Nr)
elif rarr_info[1] == 'logarithmic':
    rarr = np.logspace(log10(rst), log10(rend), Nr)
rarr = np.linspace(rst, rend, Nr)
# print(rarr)

# aarr
ast, aend, Na = float(aarr_info[2]), float(aarr_info[3]), int(aarr_info[4])
if aarr_info[1] == 'linear':
    aarr = np.linspace(ast, aend, Na)
elif aarr_info[1] == 'logarithmic':
    aarr = np.logspace(log10(ast), log10(aend), Na)
# print(aarr)

# dust temperature
Tarr_shape = (Nt, Nr, Na)
Tarr_plt = np.loadtxt('dust_temperature.txt', skiprows=1)
Tarr_plt = np.reshape(Tarr_plt, Tarr_shape)


##### settings for the plot

# find the maximum temperature to set the ylimit
Tmax = np.max(Tarr_plt)

# the time slice to be plotted for each distance
pltt_st, pltt_end = 1, Nt
pltt = np.arange(pltt_st, pltt_end, int((pltt_end - pltt_st)/3))

# find the range of r to be plotted
# distances at which all dust has sublimated will not be in the range

# pltr_st, pltr_end: the indices of the smallest and biggest distances in plot
# pltr_nonsub: the index of the smallest distance where no dust is sublimated
pltr_st = np.argwhere(Tarr_plt[pltt[0], :, -1] == 0).flatten()[-1]
pltr_nonsub = np.argwhere(Tarr_plt[pltt[2], :, 0] != 0).flatten()[0]
pltr_end = Nr

# use different interval for distances with and without sublimation
pltr = np.append(np.arange(pltr_st, pltr_nonsub, 1), np.arange(pltr_nonsub, pltr_end, 5))

# the number of lines to be plotted
NTplt = len(pltr)

# set the sublimated dust temperature to nan to show a vertical line on plot
Tarr_plt[Tarr_plt==0] = np.nan

# set up the colormap
colors = cm.viridis(np.linspace(0, 1, NTplt))[::-1]

##### plot
fig, ax = plt.subplots(figsize=(13, 8))
for i_pltr, rcolor in zip(pltr, colors):
    ax.semilogx(aarr, Tarr_plt[pltt[0], i_pltr, :], '-', color=rcolor, alpha=0.3, linewidth=2)
    vline1 = np.argwhere(np.isnan(Tarr_plt[pltt[0], i_pltr, :]))
    if (len(vline1) > 0) and (len(vline1) < Na):
        vline1 = max(vline1)+1
        ax.vlines(x=aarr[vline1], ymin=Tarr_plt[pltt[0], i_pltr, vline1], ymax=Tmax,
                   ls='-', color='k', alpha=0.3, lw=1.5)
    ax.semilogx(aarr, Tarr_plt[pltt[1], i_pltr, :], '--', color=rcolor, alpha=0.6, linewidth=2.5)
    vline2 = np.argwhere(np.isnan(Tarr_plt[pltt[1], i_pltr, :]))
    if (len(vline2) > 0) and (len(vline2) < Na):
        vline2 = max(vline2)+1
        ax.vlines(x=aarr[vline2], ymin=Tarr_plt[pltt[1], i_pltr, vline2], ymax=Tmax,
                   ls='--', color='k', alpha=0.6, lw=1.5)
    ax.semilogx(aarr, Tarr_plt[pltt[2], i_pltr, :], ':', color=rcolor, alpha=1, linewidth=3)
    vline3 = np.argwhere(np.isnan(Tarr_plt[pltt[2], i_pltr, :]))
    if (len(vline3) > 0) and (len(vline3) < Na):
        vline3 = max(vline3)+1
        ax.vlines(x=aarr[vline3], ymin=Tarr_plt[pltt[2], i_pltr, vline3], ymax=Tmax, 
                   ls=':', color='k', alpha=1, lw=1.5)
    ax.text(x=aarr[-1], y=Tarr_plt[0, i_pltr, -1], s='{:<3.2f} pc'.format(rarr[i_pltr]), size=10)
ax.set_xlim(0.7*aarr[0], 1.8*aarr[-1])
ax.set_xlabel('a [$\\mu m$]', size=15)
ax.set_ylim(ymax=Tmax)
ax.set_ylabel('T [K]', size=15)
ax.grid(True, 'both')
label_handle1 = ax.vlines(x=[], ymin=[], ymax=[], ls='-', color='k', alpha=0.3, lw=1.5,\
     label='t = {:.0f}s'.format(tarr[pltt[0]]))
label_handle2 = ax.vlines(x=[], ymin=[], ymax=[], ls='--', color='k', alpha=0.6, lw=1.5,\
     label='t = {:.0f}s'.format(tarr[pltt[1]]))
label_handle3 = ax.vlines(x=[], ymin=[], ymax=[], ls=':', color='k', alpha=1, lw=1.5,\
     label='t = {:.0f}s'.format(tarr[pltt[2]]))
ax.legend(handles=[label_handle1, label_handle2, label_handle3], fontsize=12, loc='lower left')
ax.set_title('Dust Temperature')

plt.show()