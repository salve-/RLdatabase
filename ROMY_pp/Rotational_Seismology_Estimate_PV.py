#get_ipython().magic(u'matplotlib inline')
from __future__ import print_function
import matplotlib.pylab as plt
from obspy.signal.cross_correlation import xcorr
import numpy as np
from obspy import read, read_events
from obspy.geodetics.base import gps2dist_azimuth
from obspy.taup import TauPyModel
import pylab as pl
import matplotlib.patches as mpatches
plt.style.use('ggplot')
#plt.rcParams['figure.figsize'] = 12, 8
plt.rcParams['axes.linewidth'] = 1.5




AC = read('acc_Tohoku_preproc.mseed')
RLAS = read('rot_Tohoku_preproc.mseed')
cat = read_events('xml_Tohoku.xml')
event = cat[0]

# event location from event info
source_latitude = event.origins[0].latitude
source_longitude = event.origins[0].longitude

# station location (Wettzell)
station_latitude = 49.144001
station_longitude = 12.8782

# theoretical backazimuth and distance
baz = gps2dist_azimuth(source_latitude, source_longitude, station_latitude, station_longitude)


sampling_rate = int(RLAS[0].stats.sampling_rate)
sec = 120  # window length for correlation (teleseismic event)

# calculate correlation coefficients
corrcoefs = []
for ic in xrange(0, len(RLAS[0]) // (int(sampling_rate * sec))):
        coeffs = xcorr(RLAS[0].data[sampling_rate * sec * ic : sampling_rate * sec * (ic + 1)],
                       AC[0].data[sampling_rate * sec * ic : sampling_rate * sec * (ic + 1)], 0)
        corrcoefs.append(coeffs[1])






TauPy_model = TauPyModel('ak135')
arrivals_p = TauPy_model.get_travel_times(distance_in_degree=0.001 * baz[0] / 111.11, 
                                        source_depth_in_km=event.origins[0].depth*0.001,
                                       phase_list=["P","p","Pdiff","PP","PKiKP","PKIKP","Pn","Pg"])
arrivals_s = TauPy_model.get_travel_times(distance_in_degree=0.001 * baz[0] / 111.11, 
                                        source_depth_in_km=event.origins[0].depth*0.001,
                                       phase_list=["S","s","Sdiff","SS","SKiKS","SKIKS","Sn","Sg"])
tiemp = []
tiems = []
for i in range(0,len(arrivals_p)): tiemp.append(arrivals_p[i].time)
for ii in range(0,len(arrivals_s)): tiems.append(arrivals_s[ii].time)

# first arrivals
arriv_p = min(tiemp)
arriv_s = min(tiems)
print("P-wave arrival: ", arriv_p, "sec")
print("S-wave arrival: ", arriv_s, "sec")


# **Estimate phase velocities only for time windows after S-Waves**
sampling_rate = int(RLAS[0].stats.sampling_rate)

# calculate Love wave phase velocities [km/s] for time windows featuring correlation coefficients > 0.75
surf = int(arriv_s/120.)+2
phasv = []
for iph in xrange(surf, len(corrcoefs)):
    if corrcoefs[iph] >= 0.75:
        phas_v = 0.001 * 0.5 * max(AC[0].data[sampling_rate * sec * iph : sampling_rate * sec * (iph + 1)]) / max(RLAS[0].data[sampling_rate * sec * iph : sampling_rate * sec * (iph + 1)])
    else:
        phas_v = np.NaN
    phasv.append(phas_v)


# **Plot phase velocities**

fig, (ax,ax3) = plt.subplots(2, sharex=True, figsize=(20/2.54,8/2.54))
#ax = plt.subplot(211)
#ax3 = plt.subplot(212)

#plt.figure(figsize=(15,3))
#ax=plt.subplot(111)
cmap= 'RdYlBu_r'
rot, = ax.plot(RLAS[0].times(), -0.2+RLAS[0].data/np.max(np.abs(RLAS[0].data)), 'r', label='vertical rotation rate')
acc, = ax.plot(AC[0].times(), 0.2+AC[0].data/np.max(np.abs(AC[0].data)), 'k', label='transverse acceleration')

ax.set_ylabel('norm. amp.')
ax.set_ylim(-1.0,1.2)
ax.axhline(y=.75, linewidth=1, c='k',ls='dashed')
ax.legend(loc=3, prop={"size":10}, ncol=2, bbox_to_anchor=(-0.01, -.06),framealpha=.8)
for xi in range(sec,sec * len(corrcoefs)-1, sec):
    ax.axvline(x=xi, color='.7')
    ax3.axvline(x=xi, color='.7')

ax2=ax.twinx()
sc=ax2.scatter(np.arange(60,sec * len(corrcoefs),sec),corrcoefs,c=corrcoefs,cmap=cmap,s=50,edgecolors='w',linewidth=1.5, vmin=-1, vmax=1)
ax2.set_xlim(0, RLAS[0].times()[-1])
ax2.set_ylim(-1,1.2)
ax2.set_ylabel('\nx-corr coef.')
ax2.annotate('0.75 threshold', xy=(50,.8),xycoords='data')
ax2.grid(visible=False)
ax2.spines['bottom'].set_color('k')
ax2.spines['top'].set_color('k')
ax2.spines['left'].set_color('k')
ax2.spines['right'].set_color('k')

# ax1.plot(RLAS[0].times(), RLAS[0].data)
# ax1.set_ylabel('vert. rot. rate \n[nrad/s]')

# ax2.plot(AC[0].times(), AC[0].data, 'k')
# ax2.set_ylabel('transv. acc. \n[nm/s]')
# ax2.yaxis.major.formatter.set_powerlimits((-1,2))

# ax1.axvline(arriv_p);ax1.annotate('P-arrival', xy=(arriv_p+20,np.max(RLAS[0].data)),xycoords='data');
# ax1.axvline(arriv_s);ax1.annotate('S-arrival', xy=(arriv_s+20,np.max(RLAS[0].data)),xycoords='data');
# ax2.axvline(arriv_p)
# ax2.axvline(arriv_s)    

ax3.scatter(np.arange(surf*120+60, sec * (len(phasv)+surf)+60, sec), phasv,c=corrcoefs[surf:], vmin=-1, vmax=1, s=35,cmap=cmap, marker='D',edgecolors='w',linewidth=1) ## +60 to locate in middle of window
ax3.set_xlim(0, RLAS[0].stats.delta * len(RLAS[0].data))
ax3.set_ylabel('phase vel. \n[km/s]')
ax3.set_xlabel('time [s]')

#rot = mpatches.Patch(color='red', label='vertical rotation rate')
#acc = mpatches.Patch(color='k', label='transverse acceleration')
#ax3.legend(handles=[rot,acc], loc=2, prop={'size':10})


# colorbar
cax = fig.add_axes([0.89, 0.2, 0.035, 0.32])
cbar=fig.colorbar(sc, cax=cax, orientation='vertical', ticks=[-1,0,1], label='x-corr coef.')
cbar.ax.tick_params(labelsize=9)



# add P- and S-wave arrivals 
ax3.axvline(arriv_p,linewidth=1,c='k');ax.axvline(arriv_p,linewidth=1,c='k')
ax3.annotate('P-arrival', xy=(arriv_p+20,4.5),xycoords='data');
ax3.axvline(arriv_s,linewidth=1,c='k');ax.axvline(arriv_s,linewidth=1,c='k')
ax3.annotate('S-arrival', xy=(arriv_s+20,4.5),xycoords='data');
ax3.spines['bottom'].set_color('k')
ax3.spines['top'].set_color('k')
ax3.spines['left'].set_color('k')
ax3.spines['right'].set_color('k')

pl.subplots_adjust(bottom=0.2,left=0.1,right=0.88, wspace=0.3)
plt.show()

