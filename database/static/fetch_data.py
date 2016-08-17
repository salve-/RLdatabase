#!/usr/bin/env python
# -*- coding: utf-8 -*-

from obspy.core import UTCDateTime
from obspy.core.stream import Stream
from obspy.clients.fdsn import Client as fdsnClient
import numpy as np
import matplotlib.pyplot as plt

############################### FUNCTIONS  ##################################

# function that fetches the data streams
def fetch_data(net,sta,chan,start,end,source):
	c = fdsnClient(source)
	st = c.get_waveforms(network=net, station=sta, location='', channel=chan,
                         starttime=start, endtime=end)
	return st

# function that removes instrument response and converts units
def remove_instr_resp(rt, ac):
    """
    This function removes the instrument response from the original signal
    and checks if starttime and endtime match for both instruments.
    """
    rt.detrend(type='linear')

    rt[0].data = rt[0].data * 1. / 6.3191 * 1e-3  # Rotation rate in nrad/s

    ac.detrend(type='linear')

    # TAPER
    taper_percentage = 0.05
    for comp in ac:
	    taper = np.blackman(np.int(len(comp.data) * taper_percentage))
	    taper_left, taper_right = np.array_split(taper, 2)
	    taper = np.concatenate([taper_left,
	                            np.ones(len(comp.data)-len(taper)),
	                            taper_right])
	    comp.data = comp.data * taper

    # acceleration in nm/s^2
    # note: single zero to go from velocity to acceleration
    paz_sts2 = {'poles': [(-0.0367429 + 0.036754j),
                          (-0.0367429 - 0.036754j)],
                'sensitivity': 0.944019640, 'zeros': [0j], 'gain': 1.0}
    ac.simulate(paz_remove=paz_sts2, remove_sensitivity=True)  # nm/s^2


    # make sure start and endtimes match for both instruments
    startaim = max([tr.stats.starttime for tr in (ac + rt)])
    endtaim = min([tr.stats.endtime for tr in (ac + rt)])

    ac.trim(startaim, endtaim, nearest_sample=True)
    rt.trim(startaim, endtaim, nearest_sample=True)

    return rt, ac


######################## Main program #############################

### Fetch seismic data ###


# start time from which data is fetched
starttime = UTCDateTime("2011-03-11 05:46:00")
endtime = starttime + 3600

# fetch Wettzell Ring laser data
RLAS = fetch_data('BW','RLAS', 'BJZ', starttime, endtime, 'LMU')  # Vertical rotation rate

# fetch Wettzell Broadband seismometer data (3 components, ground velocity)
BHE = fetch_data("GR","WET", "BHE", starttime, endtime, 'BGR')  # East
BHN = fetch_data("GR","WET", "BHN", starttime, endtime, 'BGR')  # North
BHZ = fetch_data("GR","WET", "BHZ", starttime, endtime, 'BGR')  # Vertical

# make a stream for simplification
AC = Stream(traces=[BHE[0],BHN[0],BHZ[0]])

# remove instrument response, convert units, broad band seismometer: velocity->acceleration
RLAS, AC = remove_instr_resp(RLAS,AC)

### Plotting part ###
fig, (ax1,ax2,ax3,ax4) = plt.subplots(4, sharex=True)

ax1.set_title(str(starttime.date) + ':   ' + str(starttime.time) + " - " + str(endtime.time))
ax1.plot(RLAS[0].times(), RLAS[0],'r',label='vertical rotation rate')
ax1.set_ylabel('rot. rate [nrad/s]')

ax2.plot(AC[0].times(),AC[0],'k',label='horizontal acc. E-component')
ax2.set_ylabel('acc. [nm/s^2]')

ax3.plot(AC[1].times(),AC[1],'k',label='horizontal acc. N-component')
ax3.set_ylabel('acc. [nm/s^2]')

ax4.plot(AC[2].times(),AC[2],'k',label='vertical acc.')
ax4.set_ylabel('acc. [nm/s^2]')
ax4.set_xlabel('time [s]')

for ax in [ax1,ax2,ax3,ax4]:
	ax.legend(loc=2, prop={"size":12})

fig.tight_layout()

plt.show()
