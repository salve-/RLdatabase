 # -*- coding: utf-8 -*-
from __future__ import print_function
import matplotlib.pylab as plt
from obspy.clients.fdsn import Client as fdsnClient
from obspy.clients.arclink.client import Client as arclinkClient
from obspy.core import UTCDateTime
from obspy.core.stream import Stream
from obspy.geodetics.base import gps2dist_azimuth
import numpy as np
from obspy.taup import TauPyModel
from obspy.signal.cross_correlation import xcorr
from obspy.signal.rotate import rotate_ne_rt
import matplotlib as mpl
import matplotlib.patches as patches
import matplotlib.lines as mlines
from pylab import Rectangle



#get_ipython().magic(u'matplotlib inline')
plt.style.use('ggplot')
plt.rcParams['figure.figsize'] = 15, 8


# **Download data streams**
# + Download event information via FDSN client
# + Download streams via Arclink client
# + Only need a one-liner to obtain each component 
# + Wettzell Ring Laser data is available on different channels 

# In[2]:


c_fdsn = fdsnClient('IRIS')
cat = c_fdsn.get_events(minmagnitude=7.0, starttime=UTCDateTime(2011,10,23), endtime=UTCDateTime(2011,10,24), catalog='GCMT')
event = cat[0]
print(cat)
print(event.event_descriptions[0]['type'], ': ',event.event_descriptions[0]['text'])

c = arclinkClient(user='test@obspy.org')

start = event.origins[0].time
print('Origin time: ', start)
end = start + 3600

RLAS = c.get_waveforms(network='BW', station='RLAS', location='', channel='BJZ', starttime=start, endtime=end)

BHE = c.get_waveforms(network='GR', station='WET', location='', channel='BHE', starttime=start, endtime=end)
BHN = c.get_waveforms(network='GR', station='WET', location='', channel='BHN', starttime=start, endtime=end)
BHZ = c.get_waveforms(network='GR', station='WET', location='', channel='BHZ', starttime=start, endtime=end)

AC = Stream(traces=[BHE[0],BHN[0],BHZ[0]])
ac = AC.copy()


# **Remove the instrument responses of the instruments from the recordings + convert units**
# - convert Ring Laser recordings to nrad/s units using a conversion factor
# - remove the seismometer response using poles and zeros + convert from velocity to acceleration [nm/s^2] in one step
# - trim the traces to make sure start- and endtimes match for both instruments

# In[3]:

RLAS.detrend(type='linear')
RLAS[0].data = RLAS[0].data * 1/6.3191 * 1e-3

AC.detrend(type='linear')
AC.taper(max_percentage=0.05)


paz_sts2 = {'poles': [(-0.0367429 + 0.036754j), (-0.0367429 - 0.036754j)],
            'sensitivity': 0.944019640,
            'zeros': [0j],
            'gain': 1.0}

AC.simulate(paz_remove=paz_sts2, remove_sensitivity=True)

startaim = max([tr.stats.starttime for tr in (AC + RLAS)])
endtaim = min([tr.stats.endtime for tr in (AC + RLAS)])

AC.trim(startaim, endtaim, nearest_sample=True)
RLAS.trim(startaim, endtaim, nearest_sample=True)


# **Resample, Filter and Rotate**

# In[4]:



RLAS.decimate(factor=4)
AC.decimate(factor=4)
f_cutoff = 1.0

RLAS.filter('lowpass', freq=f_cutoff, corners=2, zerophase=True)
AC.filter('lowpass', freq=f_cutoff, corners=2, zerophase=True)

# event location from event info
source_latitude = event.origins[0].latitude
source_longitude = event.origins[0].longitude

# station location (Wettzell)
station_latitude = 49.144001
station_longitude = 12.8782

# theoretical backazimuth and distance
baz = gps2dist_azimuth(source_latitude, source_longitude, station_latitude, station_longitude)

print('Epicentral distance [m]: ',baz[0])
print('Theoretical azimuth [deg]: ', baz[1])
print('Theoretical backazimuth [deg]: ', baz[2])

# rotate E-N component seismometer recordings to radial[1]-transverse[0] components using the theoretical BAz
AC_original = AC.copy()
#normalize 
AC_original.normalize()
RLAS.normalize()
AC.normalize()
AC.rotate(method='NE->RT',back_azimuth=baz[2])

sampling_rate = int(RLAS[0].stats.sampling_rate)
time = np.linspace(0, len(AC[0].data)/sampling_rate,len(AC[0].data))


# **Estimate correlation coefficients for different time windows and estimate BAZ**

# In[5]:



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

sampling_rate = int(RLAS[0].stats.sampling_rate)
sec = 60  # window length for correlation

# calculate correlation coefficients
corrcoefs = []
for ic in xrange(0, len(RLAS[0]) // (int(sampling_rate * sec))):
        coeffs = xcorr(RLAS[0].data[sampling_rate * sec * ic : sampling_rate * sec * (ic + 1)],
                       AC[0].data[sampling_rate * sec * ic : sampling_rate * sec * (ic + 1)], 0)
        corrcoefs.append(coeffs[1])

# estimate the Backazimuth for each time window 
step = 10
backas = np.linspace(0, 360 - step, 360 / step)
corrbaz = []
ind=None
for i6 in xrange(0, len(backas)):
    for i7 in xrange(0, len(corrcoefs)):
        corrbazz = xcorr(RLAS[0][sampling_rate * sec * i7 : sampling_rate * sec * (i7 + 1)],
                             rotate_ne_rt(AC_original.select(component='N')[0].data, 
                                          AC_original.select(component='E')[0].data, backas[i6])
                             [1][sampling_rate * sec * i7 : sampling_rate * sec * (i7 + 1)],0)
        corrbaz.append(corrbazz[1])
corrbaz = np.asarray(corrbaz)
corrbaz = corrbaz.reshape(len(backas), len(corrcoefs))

maxcorr = []
for l1 in xrange(0, len(corrcoefs)):
    maxcor_r = backas[corrbaz[:, l1].argmax()]
    maxcorr.append(maxcor_r)
maxcorr = np.asarray(maxcorr)
X, Y = np.meshgrid(np.arange(0, sec * len(corrcoefs), sec), backas)




# **Plot backazimuth estimation**

# In[6]:

# vertical rotation rate
#plt.subplot2grid((3, 30), (0, 0), colspan=29)
#plt.plot(time[700:1000], RLAS[0].data[700:1000],label='vertical rotation rate')
#plt.xlim(0, time[-1])
#plt.ylabel('vert. rot. rate \n[nrad/s]')
#plt.legend()

# add P- and S-wave arrivals 
#plt.axvline(arriv_p);plt.annotate('P-arrival', xy=(arriv_p+20,np.max(RLAS[0].data)),xycoords='data');
#plt.axvline(arriv_s);plt.annotate('S-arrival', xy=(arriv_s+20,np.max(RLAS[0].data)),xycoords='data');


def get_colors(inp, colormap, vmin=None, vmax=None):
    norm = plt.Normalize(vmin, vmax)
    return colormap(norm(inp))

colors = get_colors(np.linspace(0,36,37), plt.cm.RdYlBu) #rainbow #RdYlBu #bwr


fig=plt.figure(figsize=(20/2.54,8/2.54))
ax=plt.subplot2grid((4, 30), (0, 0), colspan=7, rowspan=2, projection='polar', frameon=False)
#ax=plt.subplot(111,  projection='polar')
theta=np.linspace(0,np.pi-np.pi/37.,37)
radii=np.ones(37)
width=np.pi/37
#r = np.arange(0, 3.0, 0.01)
#theta = 2 * np.pi * r
ax.set_theta_zero_location("N")
bars=ax.bar(theta,radii,width=width,bottom=0.8)
ax.bar(1.7453,1,width=width,bottom=0.8, color='k')
ax.bar(1.7453,1,width=width,bottom=0.9, color='r')
for iii in range(0,37):
  bars[iii].set_facecolor(colors[iii])
#ax.plot(theta, r, color='r', linewidth=3)
ax.set_yticks([])
ax.set_rmax(1.0)
ax.set_theta_direction(-1)
ax.set_xlim(0,190)
ax.set_xticks([0,np.pi/4.,np.pi/2.,np.pi*.75,np.pi])

ax.text(4.5,.6,'Backazimuth\n Grid Search',fontsize=14,fontweight='bold')

#legend
k_line = mlines.Line2D([], [], color='k',label='transv. acc.',linewidth=2)
r_line = mlines.Line2D([], [], color='r', label='vert. rot. rate',linewidth=2)
#ax.legend(handles=[k_line,r_line], loc=6)
ax.grid(True)


# transverse acceleration
tra=plt.subplot2grid((4, 30), (0, 8), colspan=21, rowspan=2)
#plt.plot(time[870*sampling_rate:1100*sampling_rate], AC[0].data[870*sampling_rate:1100*sampling_rate], 'k',label='transverse acceleration', linewidth=3)
for ii in range(3,39): # center the baz values around the estimated best value of 104 deg
  actf = AC_original.copy()
  act = actf.rotate(method='NE->RT',back_azimuth=ii*5)[0].data
  if np.abs(ii*5-baz[2])<2:
    tra.plot(time[500*sampling_rate:1100*sampling_rate], -0.1*ii+act[500*sampling_rate:1100*sampling_rate], c='k', linewidth=2.5)
  elif ii%2==0:
    tra.plot(time[500*sampling_rate:1100*sampling_rate], -0.1*ii+act[500*sampling_rate:1100*sampling_rate], c=colors[ii-3])

tra.plot(time[500*sampling_rate:1100*sampling_rate], -5+RLAS[0][500*sampling_rate:1100*sampling_rate], c='r', linewidth=1)
tra.set_xlim(900,1100)
tra.set_yticks([])

## Frame
autoAxis = tra.axis()
rec = Rectangle((autoAxis[0],autoAxis[2]),(autoAxis[1]-autoAxis[0]),(autoAxis[3]-autoAxis[2]),fill=False,lw=4,color='k')
rec = tra.add_patch(rec)
rec.set_clip_on(False)

#plt.xticks()
#plt.ylabel('norm. transv. acc. \n')
#plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
tra.legend()
tra.axvline(930,c='k',linewidth=2)
tra.axvline(1015,c='k',linewidth=2)
#plt.plot([884,889],[-4,.05], linewidth=2, c='k') 
#plt.axvline(arriv_p)
#plt.axvline(arriv_s)    
#print len(time), len(RLAS), len(AC[0])

wav=plt.subplot2grid((4, 30), (2, 0), colspan=29)
wav.plot(time, RLAS[0],'r',linewidth=.75)
wav.plot(time, 0.5+AC[0].data,'k',linewidth=.75)
#plt.yticks([0,60,120,180,240,300,360])
wav.set_xlim(0, time[-1])
wav.legend(handles=[k_line,r_line], loc=7)
#wav.set_xlabel('time [s]')
wav.add_patch(
    patches.Rectangle(
        (900, -.9),
        200,
        2.3,
        edgecolor='k',
        fill=False,   # remove background
        linewidth=2,
    )
)
# wav.set_xticklabels(())


transFigure = fig.transFigure.inverted()
coord1 = transFigure.transform(tra.transData.transform([899.5,-5.87]))
coord2 = transFigure.transform(wav.transData.transform([900,1.26]))
coord3 = transFigure.transform(tra.transData.transform([1100,-5.87]))
coord4 = transFigure.transform(wav.transData.transform([1100,1.26]))
line = mpl.lines.Line2D((coord1[0],coord2[0]),(coord1[1],coord2[1]),
                               transform=fig.transFigure,lw=2,c='k')
line2 = mpl.lines.Line2D((coord3[0],coord4[0]),(coord3[1],coord4[1]),
                               transform=fig.transFigure,lw=2,c='k')
fig.lines = line,line2,



# backazimuth estimation plot
plt.subplot2grid((4, 30), (3, 0), colspan=29)
im = plt.pcolor(X, Y, corrbaz, cmap=plt.cm.RdYlGn_r)
xx=np.arange(0, sec * len(corrcoefs) + 1, sec)
eba=np.ones(len(xx))*104
plt.plot(xx, eba, '--',c='.5', lw=2)
plt.plot(np.arange(0, sec * len(corrcoefs), sec), maxcorr, '.k',linewidth=6)
#plt.yticks([0,60,120,180,240,300,360])
plt.xlim(0, time[-1])
plt.ylim(0, 360)
plt.ylabel(u'estimated \nbackazimuth [deg]')
plt.xlabel('time [s]')
plt.text(20, 104, u'EBAz=104°',
                 bbox={'facecolor': 'k', 'alpha': 0.8}, color='w')


# add colorbar
fig = plt.subplot2grid((4, 30), (3, 29))
norm = mpl.colors.Normalize(vmin=-1, vmax=1)
cb1 = mpl.colorbar.ColorbarBase(fig, cmap=plt.cm.RdYlGn_r, norm=norm, orientation='vertical')
cb1.set_label(ur'X-corr. coeff.', fontweight='bold')
plt.subplots_adjust(hspace = .5)

plt.show()