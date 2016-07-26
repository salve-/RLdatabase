#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
ROTATIONAL SEISMOLOGY ROUTINES. Following the theory and assuming a
transversely polarized plane wave, it compares the transversal acceleration
and the vertical rotation rate of an event through:

1) direct waveform comparison in different time windows (P-coda, S-waves and
surface waves) using broadband data from WET (Wettzell) station for the
acceleration and from Wettzell ringlaser (RLAS station) for the vertical
rotation rate.

2) zero-lag cross-correlation analysis of the waveforms using the theoretical
backazimuth of the event. Aside from the correlation-coefficients, it returns
the local horizontal phase velocity for the time windows with a correlation
factor larger than 0.75. It also does the correlation analysis for different
azimuths and estimates the backazimuth as the one that has the highest
correlation coefficients.

Additionally, the routine generates a .xml file that stores data for each
event, like peak values (acceleration, rotation rate, zero-lag correlation
coefficient), signal-to-noise ratio, backazimuths.


INFORMATION: This script can be used to generate the figures for the website.
Concerning dispersion curves there are different updated versions, that include
more frequency bands and improved filtering. The problem is that automated velocity
picking does not work very well. Therefore:

+ Choose the built in velocity picker to pick phase velocities for single events 
by hand. Set: --mode velocity_picker

+ It allows to filter the to be processed events by different parameters,
including a distance-to-magnitude ratio. --> set: --filter True

+ YOU CAN ALSO READ IN QUAKEML FILES DIRECTLY, E.G. TO GET EVENTS FROM OTHER 
CATALOGES (ISC, ...) by setting: --mode qmlfile

+ YOU CAN CHOOSE EVENTS FROM THE IRIS FDSN CATALOG, which is usually faster
and more abundant especially for older events. Set: --mode fdsn

+ The events are bandstoppped for the secondary microseism (5-12s) if they are 
non-local.

+ maximum correlation coefficients for the estimated BAz are added 
in the subplot 3 on page 3.

+ only 0.5 hour recordings shown for local and close events!

+ P-coda windows (sec_p) now shorter for local and close events (2s).

+ saves peak displacement!
"""

import matplotlib as mpl
mpl.use('Agg')
import matplotlib
import matplotlib.pylab as plt
import argparse
import datetime
import json
from obspy import read_events
from obspy.clients.arclink.client import Client as arclinkClient
from obspy.clients.fdsn import Client as fdsnClient
from obspy.signal.rotate import rotate_ne_rt
from obspy.signal.cross_correlation import xcorr
from obspy.core.util.attribdict import AttribDict
from obspy.geodetics.base import gps2dist_azimuth, locations2degrees
from obspy.taup import TauPyModel
from obspy.core import read
from obspy.core.stream import Stream
# from obspy.core.event import readEvents
from obspy.imaging.mopad_wrapper import Beach
from obspy.core.utcdatetime import UTCDateTime
from mpl_toolkits.basemap import Basemap
import os
import numpy as np
import shutil

import urllib2
from xml.dom.minidom import parseString
from collections import OrderedDict


# if matplotlib.__version__ < '1.0':  # Matplotlib 1.0 or newer is necessary
#     raise ValueError('I need Matplotlib version 1.0 or newer.')


class RotationalProcessingException(Exception):
    pass


def data_from_file(net, sta, loc, chan, starttime, endtime):
    """
    Fetches seismogram data from mseed-file in a directory given by a path
    (DFM_dir).

    :type net: str
    :param net: Network code, e.g. ``'BW'``.
    :type sta: str)
    :param sta: Station code, e.g. ``'PFORL'``.
    :type loc: str
    :param loc: Location code, e.g. ``'01'``. Location code may
        contain wild cards.
    :type chan: str
    :param chan: Channel code, e.g. ``'BJZ.D'``. Channel code may
        contain wild cards.
    :type starttime: :class: `~obspy.core.utcdatetime.UTCDateTime`
    :param starttime: Starttime of the fetched seismogram.
    :type endtime: :class: `~obspy.core.utcdatetime.UTCDateTime`
    :param endtime: Endtime of the fetched seismogram.
    :return: Stream object :class: `~obspy.core.stream.Stream`
    """

    DFM_dir = "/home/johannes/waveformCompare/PFORL"  # directory of PFO data
    st = Stream()

    date = starttime
    while date <= endtime:
        print 'Loading day %s ...' % date.strftime("%d.%m")
        print 'PFO digital demodulated.'
        mseed_file = "%s/%s.%s..%s.%s.%s" % \
            (DFM_dir, net, sta, chan, date.strftime("%Y"), date.strftime("%j"))

        if os.path.exists(mseed_file):
            try:
                st += read(mseed_file)
                print 'Reading day %s successful!' % date.strftime("%Y.%m.%d")
            except IOError:
                print 'File not available for %s.' % date.strftime("%Y.%m.%d")
        else:
            print 'You chose a wrong path or the file does not exist!'
        date = date + 24 * 3600.

    return st


def download_data(origin_time, net, sta, loc, chan, source):
    """
    It downloads the data from seismic stations for the desired event(s).
    Inputs are the origin time (UTC), network, station, location and channel
    of the event. Returns a stream object fetched from Arclink. If Arclink
    does not work data is alternatively fetched from Seishub.

    :type origin_time: :class: `~obspy.core.utcdatetime.UTCDateTime`
    :param origin_time: origin time of the event.
    :type net: str
    :param net: Network code, e.g. ``'BW'``.
    :type sta: str
    :param sta: Station code, e.g. ``'WET'``.
    :type loc: str
    :param loc: Location code, e.g. ``'01'``. Location code may
        contain wild cards.
    :type chan: str
    :param chan: Channel code, e.g. ``'EHE'``. Channel code may
        contain wild cards.
    :return: Stream object :class: `~obspy.core.stream.Stream`
    """
    # try:
    #     c = arclinkClient(user='test@obspy.org')
    #     st = c.get_waveforms(network=net, station=sta, location='', channel=chan,
    #                          starttime=origin_time-190,
    #                          endtime=origin_time+3*3600+10)
    try:
        print "trying to use fdsn client service..."
        c = fdsnClient(source)
        st = c.get_waveforms(network=net, station=sta, location='', channel=chan,
                             starttime=origin_time-190,
                             endtime=origin_time+3*3600+10)

    except:
        print "trying to fetch data from file..."        
        dataDir_get = '/import/netapp-m-02-bay200/mseed_online/archive/'
        fileName = ".".join((net, sta, "." + chan + ".D",
                             origin_time.strftime("%Y.%j")))
        filePath = os.path.join(dataDir_get, origin_time.strftime("%Y"),
                                net, sta, chan + '.D', fileName)

        if os.path.isfile(filePath):
            st = read(filePath, starttime = origin_time - 180,
                      endtime = origin_time + 3 * 3600)
        else:
            print "++++ cannot find the following file: \n %s \n++++" % filePath

        if not st:
            raise RotationalProcessingException('Data not available for this'
                                                ' event...')
    st.trim(starttime=origin_time-180, endtime=origin_time+3*3600)

    print 'Download of', st[0].stats.station, st[0].stats.channel, \
        'data successful!'
    return st


def event_info_data(event, station, mode):

    """
    It extracts information from the event and generates variables containing
    the event latitude, longitude, depth, and origin time.
    Ringlaser (RLAS) and broadband signals (WET) are received from the
    download_data function.
    The great circle distance (in m and °) between event location and station
    in Wetzell, as well as the theoretical backazimuth are computed.

    :type event: :class: `~obspy.core.event.Event`
    :param event: Contains the event information.
    :type station: str
    :param station: Station from which data are fetched ('WET' or 'PFO').
    :type mode: str
    :param mode: Defines if WET data are fetched from Neries ('neries')
        or an IRIS link ('link').
    :rtype latter: float
    :return latter: Latitude of the event in degrees.
    :rtype lonter: float
    :return lonter: Longitude of the event in degrees.
    :rtype depth: float
    :return depth: Hypocenter depth in km
    :type startev: :class: `~obspy.core.utcdatetime.UTCDateTime`
    :return startev: Origin time of the event.
    :rtype rt: :class: `~obspy.core.stream.Stream`
    :return rt: Rotational signal from ringlaser.
    :rtype ac: :class: `~obspy.core.stream.Stream`
    :return ac: Three component broadband station signal.
    :rtype baz: tuple
    :return baz: [0] great circle distance in m, [1] theoretical azimuth,
        [2] theoretical backazimuth.
    :rtype gcdist: float
    :return gcdist: Great circle distance in degrees.
    """
    origin = event.preferred_origin() or event.origins[0]
    latter = origin.latitude
    lonter = origin.longitude
    startev = origin.time
    #from IPython.core.debugger import Tracer; Tracer(colors="Linux")()
    if station == 'WET' and mode == 'link':
        depth = origin.depth * 0.001  # Depth in km
    else:
        depth = origin.depth * 0.001  # Depth in km
    if station == 'WET':
        source = 'http://erde.geophysik.uni-muenchen.de' # if erde doesn't work, try 'BGR'
        net_r = 'BW'
        net_s = 'GR' #GR'
        sta_r = 'RLAS'
        sta_s = 'WET'#'WETR'
        loc_r = ''
        loc_s = ''
        if origin.time < UTCDateTime(2010,04,16):
            chan1 = 'BAZ'
        else: 
            chan1 = 'BJZ'
        chan2 = 'BHE'
        chan3 = 'BHN'
        chan4 = 'BHZ'
        # ringlaser signal
        rt = download_data(startev, net_r, sta_r, loc_r, chan1, source)
        #rt[0].data = rt[0].data * (-1) ## apply only for periods of flipped data due to instrument errors
        # broadband station signal
        acE = download_data(startev, net_s, sta_s, loc_s, chan2, source)
        acN = download_data(startev,  net_s, sta_s, loc_s, chan3, source)
        acZ = download_data(startev,  net_s, sta_s, loc_s, chan4, source)
        ac = Stream(traces=[acE[0], acN[0], acZ[0]])
        for ca in [ac[0], ac[1], ac[2], rt[0]]:
            ca.stats.coordinates = AttribDict()
            ca.stats.coordinates['longitude'] = 12.8782
            ca.stats.coordinates['latitude'] = 49.144001
            ca.stats['starttime'] = startev - 180
            ca.stats['sampling_rate'] = 20.

    else:
        net_r = 'BW'
        net_s = 'I*'
        sta_r = 'PFORL'
        sta_s = 'PFO'
        chan1 = 'BJZ.D'
        chan2 = 'BHE'
        chan3 = 'BHN'
        chan4 = 'BHZ'
        source = 'http://erde.geophysik.uni-muenchen.de'
        rt = data_from_file(net_r, sta_r, '', chan1, startev - 180,
                            startev + 3 * 3600)
        rt[0].stats.coordinates = AttribDict()
        rt[0].stats.coordinates['longitude'] = -116.453611
        rt[0].stats.coordinates['latitude'] = 33.606389
        rt[0].stats['starttime'] = startev - 180
        rt[0].stats['sampling_rate'] = 20.
        c = arclinkClient(user='test@obspy.org')
        ac = c.get_waveforms(network=net_s, station=sta_s, location='00',
                             channel='BH*', starttime=startev - 180,
                             endtime=startev + 3 * 3600, attach_response=True)

    # theoretical event backazimuth and distance
    baz = gps2dist_azimuth(latter, lonter, rt[0].stats.coordinates.latitude,
                          rt[0].stats.coordinates.longitude)
    # Great circle distance

    gcdist = locations2degrees(latter, lonter,
                               rt[0].stats.coordinates.latitude,
                               rt[0].stats.coordinates.longitude)
    return latter, lonter, depth, startev, rt, ac, baz, gcdist, net_r, net_s,\
        chan1, chan2, chan3, chan4, sta_r, sta_s, loc_r, loc_s, source


def station_components(station):

    """
    The East and North components have different labels in the WET and PFO
    data. They are standardized by this function.

    :type station: str
    :param station: Station from which data are fetched ('WET' or 'PFO').
    :rtype compE: str
    :return compE: Label for East component depending on station.
    :rtype compN: str
    return compN: Label for North component depending on station.
    """

    if station == 'WET':
        compE = 'E'
        compN = 'N'
    else:
        compE = '1'
        compN = '2'
    return compE, compN


def is_local(baz):

    """
    Checks whether the event is close (< 333.33 km), local (< 1111.1 km) or
    non-local.

    :type baz: tuple
    :param baz: Great circle distance in m, azimuth A->B in degrees,
        azimuth B->A in degrees.
    :rtype: str
    :return: Self-explaining string for event distance.
    """
    if 0.001 * baz[0] / 111.11 < 10.0:
        if 0.001 * baz[0] / 111.11 < 3.0:
            is_local = 'close'
        else:
            is_local = 'local'
    else:
        is_local = 'non-local'
    return is_local


def Get_MomentTensor_Magnitude(link):

    """
    Extracts the moment tensor and magnitude for an event if an IRIS-xml file
    is given by a url (link).

    :type link: str
    :param link: Link to the IRIS-xml, where event and moment tensor data are
        fetched.
    :rtype MomentTensor: list of floats
    :return MomentTensor: List of the six independent components of the moment
        tensor.
    :rtype Magnitude: float
    :return Magnitude: Moment magnitude (Mw) of the event.
    """
    file = urllib2.urlopen(link)
    data = file.read()
    file.close()
    dom = parseString(data)
    xmlTag = dom.getElementsByTagName('value')
    xmlTag2 = dom.getElementsByTagName('text')
    Magnitude = float(xmlTag[19].firstChild.nodeValue)
    # 19th value in xml-file
    Region = str(xmlTag2[0].firstChild.nodeValue)
    MomentTensor = []
    for el in range(1, 7):
        value = float(xmlTag[el].firstChild.nodeValue)
        MomentTensor.append(value)
    return MomentTensor, Magnitude, Region


def resample(is_local, baz, rt, ac):

    """
    Resamples signal accordingly with sampling rates and cut-off frequencies
    dependent on the location of the event (5 sec and 2Hz for local events,
    60 sec and 1 Hz for non-local events).

    :type is_local: str
    :param is_local: Self-explaining string for event distance.
    :type baz: tuple
    :param baz: Great circle distance in m, azimuth A->B in degrees,
        azimuth B->A in degrees.
    :type rt: :class: `~obspy.core.stream.Stream`
    :param rt: Rotational signal from ringlaser.
    :type ac: :class: `~obspy.core.stream.Stream`
    :param ac: Three component broadband station signal.
    :rtype rt: :class: `~obspy.core.stream.Stream`
    :return rt: Decimated rotational signal from ringlaser.
    :rtype ac: :class: `~obspy.core.stream.Stream`
    :return ac: Decimated three component broadband station signal.
    :rtype rt_pcoda: :class: `~obspy.core.stream.Stream`
    :return rt_pcoda: (Decimated) copy of rotational signal from ringlaser
        for p-coda calculation.
    :rtype ac_pcoda: :class: `~obspy.core.stream.Stream`
    :return ac_pcoda: (Decimated) copy of the three component broadband
        station signal for p-coda calculation.
    :rtype sec: int
    :return sec: Sampling rate.
    :rtype cutoff: float
    :return cutoff: Cut-off frequency for the lowpass filter.
    :rtype cutoff_pc: float
    :return cutoff_pc: Cut-off frequency for the highpass filter in P-coda.
    """

    cutoff_pc = 0.5  # Cut-off frequency for the highpass filter in the P-coda
    if is_local == 'local':
        for trr in (rt + ac):
            trr.data = trr.data[0: 1800 * rt[0].stats.sampling_rate]
        rt_pcoda = rt.copy()
        ac_pcoda = ac.copy()
        rt.decimate(factor=2)
        ac.decimate(factor=2)
        sec = 5
        cutoff = 2.0  # Cut-off freq for the lowpass filter for local events
    elif is_local == 'non-local':
        rt_pcoda = rt.copy()
        rt_pcoda.decimate(factor=2)
        ac_pcoda = ac.copy()
        ac_pcoda.decimate(factor=2)
        rt.decimate(factor=4)
        ac.decimate(factor=4)
        sec = 120
        cutoff = 1.0  # Cut-off freq for the lowpass filter for non-loc events
    else:
        for trr in (rt + ac):
            trr.data = trr.data[0: 1800 * rt[0].stats.sampling_rate]
        rt_pcoda = rt.copy()
        ac_pcoda = ac.copy()
        rt.decimate(factor=2)
        ac.decimate(factor=2)
        sec = 3
        cutoff = 4.0  # Cut-off freq for the lowpass filter for local events
    return rt, ac, rt_pcoda, ac_pcoda, sec, cutoff, cutoff_pc


def remove_instr_resp(rt, ac, rt_pcoda, ac_pcoda, station, startev):

    """
    This function removes the instrument response from the original signal
    and checks if starttime and endtime match for both instruments.

    :type rt: :class: `~obspy.core.stream.Stream`
    :param rt: Rotational signal from ringlaser.
    :type ac: :class: `~obspy.core.stream.Stream`
    :param ac: Three component broadband station signal.
    :type rt_pcoda: :class: `~obspy.core.stream.Stream`
    :param rt_pcoda: Copy of rotational signal from ringlaser
        for p-coda calculation.
    :type ac_pcoda: :class: `~obspy.core.stream.Stream`
    :param ac_pcoda: Copy of the three component broadband
        station signal for p-coda calculation.
    :type station: str
    :param station: Station from which data are fetched ('WET' or 'PFO').
    :type startev: :class: `~obspy.core.utcdatetime.UTCDateTime`
    :param startev: Origin time of the event.
    :rtype rt: :class: `~obspy.core.stream.Stream`
    :return rt: Detrended and trimmed rotational signal from ringlaser.
    :rtype ac: :class: `~obspy.core.stream.Stream`
    :return ac: Detrended and trimmed three component broadband station signal.
    :rtype rt_pcoda: :class: `~obspy.core.stream.Stream`
    :return rt_pcoda: Detrended and trimmed copy of rotational signal from
        ringlaser for p-coda calculation.
    :rtype ac_pcoda: :class: `~obspy.core.stream.Stream`
    :return ac_pcoda: Detrended and trimmed copy of the three component
        broadband station signal for p-coda calculation.
    """

    sta_s = 'WET' #'WETR'
    
    rt.detrend(type='linear')
    rt_pcoda.detrend(type='linear')
    if station == 'WET':
        rt[0].data = rt[0].data * 1. / 6.3191 * 1e-3  # Rotation rate in nrad/s
        rt_pcoda[0].data = rt_pcoda[0].data * 1. / 6.3191 * 1e-3

        ac.detrend(type='linear')
        ac_pcoda.detrend(type='linear')
        # TAPER
        ac.taper(max_percentage=0.05)
        rt.taper(max_percentage=0.05)
        # acceleration in nm/s^2
        # note: single zero to go from velocity to acceleration

        displacement = ac.copy() # copy 'raw' velocity measurements to use as displacement later

        paz_sts2 = {'poles': [(-0.0367429 + 0.036754j),
                              (-0.0367429 - 0.036754j)],
                    'sensitivity': 0.944019640, 'zeros': [0j], 'gain': 1.0}
        paz_sts2_displacement = {'poles': [(-0.0367429 + 0.036754j),
                              (-0.0367429 - 0.036754j)],
                    'sensitivity': 0.944019640, 'zeros': [0j,0j,0j], 'gain': 1.0}
        paz_lennartz = {'poles': [(-0.22 + 0.235j),
                              (-0.22 - 0.235j), (-0.23 + 0.0j)],
                    'sensitivity': 1, 'zeros': [(0+0j),(0+0j)], 'gain': 1.0}
        paz_lennartz_displacement = {'poles': [(-0.22 + 0.235j),
                              (-0.22 - 0.235j), (-0.23 + 0.0j)],
                    'sensitivity': 1, 'zeros': [(0+0j),(0+0j),(0+0j),(0+0j)], 'gain': 1.0}
        if sta_s == 'WETR':
            ac.simulate(paz_remove=paz_lennartz, remove_sensitivity=True)  # nm/s^2
            ac_pcoda.simulate(paz_remove=paz_lennartz, remove_sensitivity=True)
            displacement.simulate(paz_remove=paz_lennartz_displacement, remove_sensitivity=True)
            ac.filter('highpass', freq=0.04, zerophase=True, corners=3)
            ac_pcoda.filter('highpass', freq=0.04, zerophase=True, corners=3)
            displacement.filter('highpass', freq=0.04, zerophase=True, corners=3)
        else:
            ac.simulate(paz_remove=paz_sts2, remove_sensitivity=True)  # nm/s^2
            ac_pcoda.simulate(paz_remove=paz_sts2, remove_sensitivity=True)
            displacement.simulate(paz_remove=paz_sts2_displacement, remove_sensitivity=True)

    else:
        rt[0].data = rt[0].data * 1. / 2.5284 * 1e-3  # Rotation rate in nrad/s
        rt_pcoda[0].data = rt_pcoda[0].data * 1. / 2.5284 * 1e-3  # Rotation
        # rate in nrad/s
        ac.detrend(type='linear')
        ac_pcoda.detrend(type='linear')

        displacement = ac.copy()
        # TAPER
        ac.taper(max_percentage=0.05)
        rt.taper(max_percentage=0.05)

        ac.remove_response(output='ACC', pre_filt=(0.005, 0.006, 30., 35.))
        displacement.remove_response(output='DISP', pre_filt=(0.005, 0.006, 30., 35.))
        ac_pcoda.remove_response(output='VEL',
                                 pre_filt=(0.005, 0.006, 30., 35.))

        # to nm/s^2
        for traza in (ac + ac_pcoda):
            traza.data = 1e9 * traza.data

    # make sure start and endtimes match for both instruments
    startaim = max([tr.stats.starttime for tr in (ac + rt)])
    endtaim = min([tr.stats.endtime for tr in (ac + rt)])

    ac.trim(startaim, endtaim, nearest_sample=True)
    displacement.trim(startaim, endtaim, nearest_sample=True)
    rt.trim(startaim, endtaim, nearest_sample=True)
    ac_pcoda.trim(startaim, endtaim, nearest_sample=True)
    rt_pcoda.trim(startaim, endtaim, nearest_sample=True)

    return rt, ac, rt_pcoda, ac_pcoda, displacement

def gaussianfilter(sigarray, delta, bandwidth, freq0):
    """
    sigarray = signal array (much faster if the length of this is a power of 2)
    delta    = time sampling interval (seconds)
    bandwidth    = filter df (>0.)
    freq0    = center frequency (Hz)
    """
    #1 : prepare the frequency domain filter
    n = len(sigarray)  #number of samples
    freq = fftfreq(n, delta) #exact (!) frequency array
    # we construct our gaussian according the constQ criterion of Archambeau et al.
    beta = np.log(2.)/2.
    g = np.sqrt(beta/np.pi)*np.exp(-beta * (np.abs(freq - freq0) / bandwidth) ** 2.) #do not forget negative frequencies


    #2 : convolve your signal by the filter in frequency domain 
    sigarray_fourier = fft(sigarray) 
    sigarray_fourier_filtered = sigarray_fourier * g

    #3 : back to time domain
    sigarray_filtered = np.real(ifft(sigarray_fourier_filtered))
    sigarray_filtered = highpass(sigarray_filtered, freq=0.0033, df=5, corners=3, zerophase=True)
    sigarray_filtered = detrend(sigarray_filtered)
    return sigarray_filtered

def filter_and_rotate(ac, rt, baz, rt_pcoda, ac_pcoda, cutoff, cutoff_pc,
                      station, is_local):
    """
    Filters trace data using the cut-off frequencies and
    lowpass/ highpass and zerophase filters. Rotates the horizontal components
    of the signal to theoretical backazimuth.

    :type rt: :class: `~obspy.core.stream.Stream`
    :param rt: Rotational signal from ringlaser.
    :type ac: :class: `~obspy.core.stream.Stream`
    :param ac: Three component broadband station signal.
    :type baz: tuple
    :param baz: Great circle distance in m, azimuth A->B in degrees,
        azimuth B->A in degrees.
    :type rt_pcoda: :class: `~obspy.core.stream.Stream`
    :param rt_pcoda: Copy of rotational signal from ringlaser
        for p-coda calculation.
    :type ac_pcoda: :class: `~obspy.core.stream.Stream`
    :param ac_pcoda: Copy of the three component broadband
        station signal p-coda calculation.
    :type cutoff: float
    :param cutoff: Cut-off frequency for the lowpass filter.
    :type cutoff_pc: float
    :param cutoff_pc: Cut-off frequency for the highpass filter in P-coda.
    :type station: str
    :param station: Station from which data are fetched ('WET' or 'PFO').
    :rtype rotate: :class: `~obspy.core.stream.Stream`
    :return rotate: Stream object of the broadband station signal with
        rotated horizontal components.
    :rtype pcod_rotate: :class: `~obspy.core.stream.Stream`
    :return pcod_rotate: Stream object of the broadband station signal with
        rotated horizontal components for P-coda calculations.
    :rtype pcoda_rotate: :class: `~obspy.core.stream.Stream`
    :return pcoda_rotate: Stream object of the broadband station signal with
        rotated horizontal components for P-coda calculations.
    :rtype frtp: :class: `~obspy.core.stream.Stream`
    :return frtp: Filtered rt_pcoda.
    :rtype facp: :class: `~obspy.core.stream.Stream`
    :return facp: Filtered ac_pcoda.
    :rtype frotate: :class: `~obspy.core.stream.Stream`
    :return frotate: Filtered and rotated ac_pcoda.
    :rtype cop_rt: :class: `~obspy.core.stream.Stream`
    :return cop_rt: Highpass filtered rotational trace (rt).
    """
    compE, compN = station_components(station)

    # new ac_band for rotate_band

    ac_band1 = ac.copy()
    ac_band2 = ac.copy()
    ac_band3 = ac.copy()
    ac_band4 = ac.copy()
    ac_band5 = ac.copy()
    ac_band6 = ac.copy()
    ac_band7 = ac.copy()
    ac_band8 = ac.copy()
    rt1 = rt.copy()
    rt2 = rt.copy()
    rt3 = rt.copy()
    rt4 = rt.copy()
    rt5 = rt.copy()
    rt6 = rt.copy()
    rt7 = rt.copy()
    rt8 = rt.copy()

    cop_ac = ac.copy()
    cop_rt = rt.copy()
    cop_ac.filter('highpass', freq=cutoff_pc, corners=2, zerophase=True)
    cop_rt.filter('highpass', freq=cutoff_pc, corners=2, zerophase=True)
    ac.filter('lowpass', freq=cutoff, corners=2, zerophase=True)
    rt.filter('lowpass', freq=cutoff, corners=2, zerophase=True)
    ac.filter('highpass', freq=0.005, corners=2, zerophase=True)
    rt.filter('highpass', freq=0.005, corners=2, zerophase=True)

    # filter out secondary microseism only for non-local events
    if is_local == "non-local":    
        ac.filter('bandstop', freqmin=0.083, freqmax=0.2, corners=4, zerophase=True)
        rt.filter('bandstop', freqmin=0.083, freqmax=0.2, corners=4, zerophase=True)
        ac.taper(max_percentage=0.05)
        rt.taper(max_percentage=0.05)


    # rotate translational signal to theoretical event backazimuth
    rotate = rotate_ne_rt(ac.select(component=compN)[
                          0].data, ac.select(component=compE)[0].data, baz[2])
    pcod_rotate = rotate_ne_rt(cop_ac.select(component=compN)[
        0].data, cop_ac.select(component=compE)[0].data, baz[2])
    pcoda_rotate = rotate_ne_rt(ac_pcoda.select(component=compN)[0].data,
                                ac_pcoda.select(component=compE)[0].data,
                                baz[2])
    # for each frequency band we need a separate rt and rotate:
    rt_band1 = rt1.filter('bandpass', freqmin=0.01, freqmax=0.02,
                          corners=3, zerophase=True)
    ac_band1 = ac_band1.filter('bandpass', freqmin=0.01, freqmax=0.02,
                               corners=3, zerophase=True)
    rotate_band1 = rotate_ne_rt(ac_band1.select(component=compN)[
        0].data, ac_band1.select(component=compE)[0].data, baz[2])
    rt_band2 = rt2.filter('bandpass', freqmin=0.02, freqmax=0.04,
                          corners=3, zerophase=True)
    ac_band2 = ac_band2.filter('bandpass', freqmin=0.02, freqmax=0.04,
                               corners=3, zerophase=True)
    rotate_band2 = rotate_ne_rt(ac_band2.select(component=compN)[
        0].data, ac_band2.select(component=compE)[0].data, baz[2])
    rt_band3 = rt3.filter('bandpass', freqmin=0.04, freqmax=0.1,
                          corners=3, zerophase=True)
    ac_band3 = ac_band3.filter('bandpass', freqmin=0.04, freqmax=0.1,
                               corners=3, zerophase=True)
    rotate_band3 = rotate_ne_rt(ac_band3.select(component=compN)[
        0].data, ac_band3.select(component=compE)[0].data, baz[2])
    rt_band4 = rt4.filter('bandpass', freqmin=0.1, freqmax=0.2,
                          corners=3, zerophase=True)
    ac_band4 = ac_band4.filter('bandpass', freqmin=0.1, freqmax=0.2,
                               corners=3, zerophase=True)
    rotate_band4 = rotate_ne_rt(ac_band4.select(component=compN)[
        0].data, ac_band4.select(component=compE)[0].data, baz[2])
    rt_band5 = rt5.filter('bandpass', freqmin=0.2, freqmax=0.3,
                          corners=3, zerophase=True)
    ac_band5 = ac_band5.filter('bandpass', freqmin=0.2, freqmax=0.3,
                               corners=3, zerophase=True)
    rotate_band5 = rotate_ne_rt(ac_band5.select(component=compN)[
        0].data, ac_band5.select(component=compE)[0].data, baz[2])
    rt_band6 = rt6.filter('bandpass', freqmin=0.3, freqmax=0.4,
                          corners=3, zerophase=True)
    ac_band6 = ac_band6.filter('bandpass', freqmin=0.3, freqmax=0.4,
                               corners=3, zerophase=True)
    rotate_band6 = rotate_ne_rt(ac_band6.select(component=compN)[
        0].data, ac_band6.select(component=compE)[0].data, baz[2])
    rt_band7 = rt7.filter('bandpass', freqmin=0.4, freqmax=0.6,
                          corners=3, zerophase=True)
    ac_band7 = ac_band7.filter('bandpass', freqmin=0.4, freqmax=0.6,
                               corners=3, zerophase=True)
    rotate_band7 = rotate_ne_rt(ac_band7.select(component=compN)[
        0].data, ac_band7.select(component=compE)[0].data, baz[2])
    rt_band8 = rt8.filter('bandpass', freqmin=0.6, freqmax=1.0,
                          corners=3, zerophase=True)
    ac_band8 = ac_band8.filter('bandpass', freqmin=0.6, freqmax=1.0,
                               corners=3, zerophase=True)
    rotate_band8 = rotate_ne_rt(ac_band8.select(component=compN)[
        0].data, ac_band8.select(component=compE)[0].data, baz[2])
    ###

    frtp = rt_pcoda.copy()
    frtp.filter('highpass', freq=cutoff_pc, corners=2, zerophase=True)
    facp = ac_pcoda.copy()
    facp.filter('highpass', freq=cutoff_pc, corners=2, zerophase=True)
    frotate = rotate_ne_rt(facp.select(component=compN)[0].data,
                           facp.select(component=compE)[0].data, baz[2])
    print 'Done'
    return rotate, pcod_rotate, pcoda_rotate, frtp, facp, frotate, cop_rt,\
        rt_band1, rt_band2, rt_band3, rt_band4, rt_band5, rt_band6, rt_band7,\
        rt_band8, rotate_band1, rotate_band2, rotate_band3, rotate_band4,\
        rotate_band5, rotate_band6, rotate_band7, rotate_band8


def ps_arrival_times(distance, depth, init_sec):

    """
    Obtains the arrival times (in seconds after the start time of the fetched
    data) of the first P an S waves of the event. The inputs are the
    epicentral distance in degrees, the depth in km and the initial time in
    seconds (starttime_of_the_event - data_starttime)

    :type distance: float
    :param distance: Great circle distance between earthquake source and
        receiver station.
    :type depth: float
    :param depth: Hypocenter depth in km.
    :type init_sec: float
    :param init_sec: Initial time of the event in sec in the fetched data.
    :rtype arriv_p: float
    :return arriv_p: Arrival time of the first P-wave.
    :rtype arriv_s: float
    :return arriv_s: Arrival time of the first S-wave.
    """
    # use taup to get the theoretical arrival times for P & S
    TauPy_model = TauPyModel('iasp91')
    tt = TauPy_model.get_travel_times(
        distance_in_degree=0.001 * distance / 111.11, source_depth_in_km=depth)
    tiemp = []
    tiems = []
    # from all possible P arrivals select the earliest one
    for i2 in xrange(0, len(tt)):
        if tt.__getitem__(i2).__dict__['name'] == 'P' or tt.__getitem__(i2).__dict__[
                'name'] == 'p' or tt.__getitem__(i2).__dict__['name'] ==\
                'Pdiff' or tt.__getitem__(i2).__dict__['name'] == 'PKiKP' or\
                tt.__getitem__(i2).__dict__['name'] == 'PKIKP' or tt.__getitem__(i2).__dict__[
                'name'] == 'PP' or tt.__getitem__(i2).__dict__['name'] ==\
                'Pb' or tt.__getitem__(i2).__dict__['name'] == 'Pn' or\
                tt.__getitem__(i2).__dict__['name'] == 'Pg':
                    tiem_p = tt.__getitem__(i2).__dict__['time']
                    tiemp.append(tiem_p)
    arriv_p = np.floor(init_sec + min(tiemp))

    # from all possible S arrivals select the earliest one
    for i3 in xrange(0, len(tt)):
        if tt.__getitem__(i3).__dict__['name'] == 'S' or tt.__getitem__(i3).__dict__[
                'name'] == 's' or tt.__getitem__(i3).__dict__['name'] ==\
                'Sdiff' or tt.__getitem__(i3).__dict__['name'] == 'SKiKS' or\
                tt.__getitem__(i3).__dict__['name'] == 'SKIKS' or tt.__getitem__(i3).__dict__[
                'name'] == 'SS' or tt.__getitem__(i3).__dict__['name'] ==\
                'Sb' or tt.__getitem__(i3).__dict__['name'] == 'Sn' or\
                tt.__getitem__(i3).__dict__['name'] == 'Sg':
                    tiem_s = tt.__getitem__(i3).__dict__['time']
                    tiems.append(tiem_s)
    arriv_s = np.floor(init_sec + min(tiems))
    return arriv_p, arriv_s


def time_windows(baz, arriv_p, arriv_s, init_sec, is_local):
    """
    Determines time windows for arrivals and subplots for P-waves,
    S-waves, initial and latter surface waves.

    :type baz: tuple
    :param baz: Great circle distance in m, azimuth A->B in degrees,
        azimuth B->A in degrees.
    :type arriv_p: float
    :param arriv_p: Arrival time of the first P-wave.
    :type arriv_s: float
    :param arriv_s: Arrival time of the first S-wave.
    :type init_sec: float
    :param init_sec: Initial time of the event in sec in the fetched data.
    :type is_local: str
    :param is_local: Self-explaining string for event distance.
    :rtype min_pw: float
    :return min_pw: Starttime for P-waves window.
    :rtype max_pw: Endtime for P-waves window.
    :return min_sw: Starttime for S-waves window.
    :rtype max_sw: Endtime for S-waves window.
    :return min_lwi: Starttime for initial surface-waves window.
    :rtype max_lwi: Endtime for initial surface-waves window.
    :return min_lwf: Starttime for latter surface-waves window.
    :rtype max_lwf: Endtime for latter surface-waves window.
    """

    # TIME WINDOWS (for arrivals and subplots)
    # Window lengths dependent on event distance
    if is_local == 'non-local':
        min_pw = arriv_p
        max_pw = min_pw + (arriv_s - arriv_p) // 4
        min_sw = arriv_s - 0.001 * (arriv_s - arriv_p)
        max_sw = arriv_s + 150
        min_lwi = surf_tts(baz[0], init_sec) - 20
        t1 = (baz[0]/1000000) * 50
        # window length grows 50 sec per 1000 km.
        max_lwi = min_lwi + t1
        min_lwf = max_lwi
        t2 = (baz[0]/1000000) * 60
        # window length grows 60 sec per 1000 km.
        max_lwf = min_lwf + t2
    elif is_local == 'local':
        min_pw = arriv_p
        max_pw = min_pw + 20
        min_sw = arriv_s - 5
        max_sw = min_sw + 20
        min_lwi = surf_tts(baz[0], init_sec) + 20
        max_lwi = min_lwi + 50
        min_lwf = max_lwi
        max_lwf = min_lwf + 80
    else:
        min_pw = arriv_p
        max_pw = min_pw + 7
        min_sw = arriv_s
        max_sw = min_sw + 7
        min_lwi = surf_tts(baz[0], init_sec) + 5
        max_lwi = min_lwi + 12
        min_lwf = max_lwi
        max_lwf = min_lwf + 80
    print(min_pw, min_lwi)
    return min_pw, max_pw, min_sw, max_sw, min_lwi, max_lwi, min_lwf, max_lwf


def surf_tts(distance, start_time):

    """
    Uses arrival times for different epicentral distances based on the IASP91
    travel times model to estimate a curve of travel times for surface waves
    and get the arrival time of the surface waves of the event. Inputs are the
    epicentral distance in degrees and the time in seconds at which the event
    starts in the fetched data.

    :type distance: float
    :param distance: Epicentral distance in degrees between earthquake source
        and receiver station.
    :type start_time: float
    :param start_time: Starttime of the event in the fetched seismogram.
    :rtype arrival: float
    :return arrival: Arrival time of the surface waves of the event.
    """
    deltas = np.arange(0., 140., 5.)
    tts = 60. * np.array(
        [0., 2., 4., 6.2, 8.4, 11., 13., 15.2, 17.8, 19.4, 22., 24.1, 26.6,
         28.6, 30.8, 33., 35.6, 37.4, 39.8, 42., 44.2, 46.4, 48.8, 50.9, 53.6,
         55.2, 57.8, 60.])
    (mval, nval) = np.polyfit(deltas, tts, 1)
    # calculate surface wave travel times for degrees 1 to 180 ?
    surftts = mval * np.arange(0., 180.1, 0.01)
    difer = []
    for i4 in xrange(0, len(surftts)):
        dife_r = abs(0.001 * distance / 111.11 - np.arange(0., 180.1, 0.01)
                     [i4])
        difer.append(dife_r)
    # love wave arrival: event time + surftts for closest degree??
    # (smallest difference between distance for surftts and actual distance of
    #  event)
    arriv_lov = np.floor(start_time + surftts[np.asarray(difer).argmin()])
    diferans = []
    for i1 in xrange(len(deltas)):
        dif2 = abs(np.arange(0., 180.1, 0.01)[np.asarray(difer).argmin()] -
                   deltas[i1])
        diferans.append(dif2)
    # arrival = love wave arrival - p arrival?
    peq = surftts[np.asarray(difer).argmin()] - \
        tts[np.asarray(diferans).argmin()]
    arrival = arriv_lov + peq
    return arrival


def Get_corrcoefs(rotra, rodat, acstr, rotate_array, sec, station):

    """
    Calculates the zero-lag correlation coefficients between the ringlaser
    data and the broadband station data.

    :type rotra: :class: `~obspy.core.trace.Trace`
    :param rotra: Trace of the rotational data from ringlaser.
    :type rodat: numpy.ndarray
    :param rodat: Rotational data ...
    :type acstr: :class: `~obspy.core.stream.Stream`
    :param acstr: Three component broadband station signal.
    :type rotate_array: numpy.ndarray
    :param rotate_array:
    :type sec: int
    :param sec: Sampling rate.
    :type station: str
    :param station: Station from which data are fetched ('WET' or 'PFO').
    :rtype corrcoefs: numpy.ndarray
    :return corrcoefs: Correlation coefficients.
    :rtype thres: numpy.ndarray
    :return thres: Array for plotting dashed line of 75'%' correlation.
    """

    compE, compN = station_components(station)
    corrcoefs = []
    for i5 in xrange(0, len(rodat) // (int(rotra.stats.sampling_rate) * sec)):
        coeffs = xcorr(rodat[rotra.stats.sampling_rate * sec *
                             i5:rotra.stats.sampling_rate * sec * (i5 + 1)],
                       rotate_array[acstr.select(component=compN)[0].stats.
                                    sampling_rate * sec * i5:acstr.
                                    select(component=compN)[0].stats.
                                    sampling_rate * sec * (i5 + 1)], 0)
        corrcoefs.append(coeffs[1])
    corrcoefs = np.asarray(corrcoefs)
    thres = 0.75 * np.ones(len(corrcoefs) + 1)
    return corrcoefs, thres


def backas_analysis(rotra, rodat, acstr, sec, corrcoefs, ind, station):

    """
    Backazimuth analysis: Computes the correlation coefficients for
    the backazimuth and backazimuth values.

    :type rotra: :class: `~obspy.core.trace.Trace`
    :param rotra: Trace of the rotational data from ringlaser.
    :type rodat: numpy.ndarray
    :param rodat: Rotational data ...
    :type acstr: :class: `~obspy.core.stream.Stream`
    :param acstr: Three component broadband station signal.
    :type sec: int
    :param sec: Sampling rate.
    :type corrcoefs: numpy.ndarray
    :param corrcoefs: Correlation coefficients.
    :type ind: int
    :param ind: Index for stream data selection.
    :type station: str
    :param station: Station from which data are fetched ('WET' or 'PFO').
    :rtype corrbaz: numpy.ndarray
    :return corrbaz: Correlation coefficients for each backazimuth.
    :rtype maxcorr: numpy.ndarray
    :return maxcorr: Backazimuth values for maximum correlation for each time
        window.
    :rtype backas: numpy.ndarray
    :return backas: Vector containing backazimuths (step: 10°).
    """
    compE, compN = station_components(station)
    step = 10
    backas = np.linspace(0, 360 - step, 360 / step)
    corrbaz = []
    for i6 in xrange(0, len(backas)):
        for i7 in xrange(0, len(corrcoefs)):
            corrbazz = xcorr(rodat[rotra.stats.sampling_rate * sec *
                                   i7:rotra.stats.sampling_rate * sec *
                                   (i7 + 1)],
                             rotate_ne_rt(acstr.select(component=compN)[0].
                                          data[0:ind], acstr.
                                          select(component=compE)[0].
                                          data[0:ind], backas[i6])
                             [1][rotra.stats.sampling_rate * sec * i7:rotra.
                                 stats.sampling_rate * sec * (i7 + 1)],
                             0)
            corrbaz.append(corrbazz[1])
    corrbaz = np.asarray(corrbaz)
    corrbaz = corrbaz.reshape(len(backas), len(corrcoefs))

    maxcorr = []
    for l1 in xrange(0, len(corrcoefs)):
        maxcor_r = backas[corrbaz[:, l1].argmax()]
        maxcorr.append(maxcor_r)
    maxcorr = np.asarray(maxcorr)

    coefs = [] # array containing max. x-corr coef for each window
    for m1 in range(0,len(corrcoefs)):
        coefs.append(np.max(corrbaz[:,m1]))

    return corrbaz, maxcorr, backas, coefs


def backas_est(rt, ac, min_sw, max_lwf, station):
    """
    Calculates the sum of all correlation coefficients above a certain
    threshold (0.9) within S-waves and surface waves for each backazimuth.

    :type rt: :class: `~obspy.core.stream.Stream`
    :param rt: Rotational signal from ringlaser.
    :type ac: :class: `~obspy.core.stream.Stream`
    :param ac: Three component broadband station signal.
    :type min_sw: float
    :param min_sw: Starttime for S-waves window.
    :type max_lwf: float
    :param max_lwf: Endtime for latter surface waves window.
    :type station: str
    :param station: Station from which data are fetched ('WET' or 'PFO').
    :rtype corrsum: list of floats
    :return corrsum: Sum of all correlation coefficients above a certain
        threshold (0.9) within S-waves and surface waves for each backazimuth.
    :rtype backas2: numpy.ndarray
    :return backas2: Vector containing backazimuths (step: 1°).
    """
    compE, compN = station_components(station)
    rt2 = rt[0].data[min_sw * rt[0].stats.sampling_rate:
                     max_lwf * rt[0].stats.sampling_rate]
    acn2 = ac.select(component=compN)[0].data[min_sw *
                                              rt[0].stats.sampling_rate:max_lwf
                                              * rt[0].stats.
                                              sampling_rate]
    ace2 = ac.select(component=compE)[0].data[min_sw *
                                              rt[0].stats.sampling_rate:max_lwf
                                              * rt[0].stats.sampling_rate]
    sec2 = 30
    step2 = 1
    backas2 = np.linspace(0, 360 - step2, 360 / step2) # BAz array
    corrbaz2 = []

    for j3 in xrange(len(backas2)):
        for j4 in xrange(len(rt2) // (int(rt[0].stats.sampling_rate) * sec2)):
            corrbazz2 = xcorr(rt2[int(rt[0].stats.sampling_rate) * sec2 *
                                  j4:int(rt[0].stats.sampling_rate) * sec2 *
                                  (j4 + 1)],
                              rotate_ne_rt(acn2, ace2, backas2[j3])[1]
                              [int(rt[0].stats.sampling_rate) * sec2 * j4:
                               int(rt[0].stats.sampling_rate) * sec2 *
                               (j4 + 1)],
                              0)
            corrbaz2.append(corrbazz2[1]) # correlation coefficients for BAz-array

    corrbaz2 = np.asarray(corrbaz2)
    corrbaz2 = corrbaz2.reshape(len(backas2), len(corrbaz2) / len(backas2))
    corrsum = []
    
    for j1 in xrange(len(corrbaz2[:, 0])):
        bazsum = []
        for j2 in xrange(len(corrbaz2[0, :])):
            if corrbaz2[j1, j2] >= 0.9:
                corradj = corrbaz2[j1, j2]
            else:
                corradj = 0.0
            bazsum.append(corradj)
        bazsum = np.asarray(bazsum)
        bazsum = sum(bazsum)
        corrsum.append(bazsum)

    best_ebaz = backas2[np.asarray(corrsum).argmax()] ## = EBA!
    max_ebaz_xcoef = np.max(corrbaz2[best_ebaz]) ## maximum correlation coefficient for EBA

    return corrsum, backas2, max_ebaz_xcoef, best_ebaz


def phase_vel(rt, sec, corrcoefs, rotate, corrsum, backas2, ind_band,
              ind_surf):

    """
    Calculates the phase velocities and the estimated backazimuth.

    :type rt: :class: `~obspy.core.stream.Stream`
    :param rt: Rotational signal from ringlaser.
    :type sec: int
    :param sec: Sampling rate.
    :type corrcoefs: numpy.ndarray
    :param corrcoefs: Correlation coefficients.
    :type rotate: :class: `~obspy.core.stream.Stream`
    :param rotate: Stream object of the broadband station signal with
        rotated horizontal components.
    :type corrsum: list of floats
    :param corrsum: Sum of all correlation coefficients above a certain
        threshold (0.9) within S-waves and surface waves for each backazimuth.
    :type backas2: numpy.ndarray
    :param backas2: Vector containing backazimuths (step: 1°).
    :rtype phasv: numpy.ndarray
    :return phasv: Phase velocities of the seismic signal.
    :rtype EBA: float
    :return EBA: Estimated backazimuth.
    """

    phasv = []
    if not ind_band:  # not dealing with frequency bands
        for i8 in xrange(0, len(corrcoefs)):
            if corrcoefs[i8] >= 0.75:
                # Velocity in km/s
                phas_v = .001 * 0.5 * max(rotate[1][rt[0].stats.sampling_rate *
                                          sec * i8:rt[0].stats.sampling_rate *
                                          sec * (i8 + 1)]) /\
                    max(rt[0].data[rt[0].stats.sampling_rate * sec *
                        i8:rt[0].stats.sampling_rate * sec * (i8 + 1)])
            else:
                phas_v = np.NaN
            phasv.append(phas_v)
        phasv = np.asarray(phasv)  # Velocity in km/s

    if ind_band:  # dealing with frequency bands
        for i8 in xrange(ind_surf, len(corrcoefs)):
            if corrcoefs[i8] >= 0.75:
                # Velocity in km/s
                phas_v = .001 * 0.5 * max(rotate[1][rt[0].stats.sampling_rate *
                                          sec * i8:rt[0].stats.sampling_rate *
                                          sec * (i8 + 1)]) /\
                    max(rt[0].data[rt[0].stats.sampling_rate * sec *
                        i8:rt[0].stats.sampling_rate * sec * (i8 + 1)])
            else:
                phas_v = np.NaN
            phasv.append(phas_v)
        phasv = np.asarray(phasv)  # Velocity in km/s

    if max(np.asarray(corrsum)) == 0.0:
        EBA = np.nan
    else:
        # Estimated backazimuth [°]
        EBA = backas2[np.asarray(corrsum).argmax()]
    return phasv, EBA


def sn_ratio(full_signal, p_arrival, sam_rate):

    """
    Characterizes the signal-to-noise ratio of the event(s) as the ratio of
    the peak amplitude of the whole wave train and the mean amplitude in a
    noise window before the first theoretical arrival, assuming that the noise
    has the same behavior in all the data. The inputs are the data, the
    theoretical time of the first P-arrival (as seconds after the first sample
    of the fetched data) and the sampling rate.

    :type full_signal: numpy.ndarray
    :param full_signal: Amplitude data of the full signal.
    :type p_arrival: float
    :param p_arrival: Arrival time of the first P-wave.
    :type sam_rate: float
    :param sam_rate: Sampling rate.
    :rtype SNR: float
    :return SNR: Signal-to-noise ratio of the seismogram.
    """

    tr_sign = max(full_signal)
    tr_noise = abs(np.mean(full_signal[sam_rate * (p_arrival - 180): sam_rate
                   * (p_arrival - 100)]))
    SNR = tr_sign/tr_noise
    return SNR


def store_info_json(rotate, ac, rt, corrcoefs, baz, arriv_p, EBA, tag_name,
                    station, phasv_means, phasv_stds, startev, event, net_r,
                    net_s, chan1, chan2, chan3, chan4, sta_r, sta_s, source,
                    loc_r, loc_s, event_source, depth, displ, magnitude, 
                    distance, max_ebaz_xcoef):

    """
    Generates a human readable .json file that stores data for each event,
    like peak values (acceleration, rotation rate,
    zero-lag correlation coefficient), signal-to-noise ratio, backazimuths.

    :type rotate: :class: `~obspy.core.stream.Stream`
    :param rotate: Stream object of the broadband station signal with
        rotated horizontal components.
    :type ac: :class: `~obspy.core.stream.Stream`
    :param ac: Three component broadband station signal.
    :type rt: :class: `~obspy.core.stream.Stream`
    :param rt: Rotational signal from ringlaser.
    :type corrcoefs: numpy.ndarray
    :param corrcoefs: Correlation coefficients.
    :type baz: tuple
    :param baz: Great circle distance in m, azimuth A->B in degrees,
        azimuth B->A in degrees.
    :type arriv_p: float
    :param arriv_p: Arrival time of the first P-wave.
    :type EBA: float
    :param EBA: Estimated backazimuth.
    :type tag_name: string
    :param tag_name: Name of the .json file.
    :type station: str
    :param station: Station from which data are fetched ('WET' or 'PFO').
    """
    compE, compN = station_components(station)
    mxpt = np.argmax(rt[0].data, axis=0)
    listn = []
    for val in range(mxpt, mxpt + 100):
        if (rt[0].data[val] > 0):
            listn.append(val)
        else:
            break
    for val2 in range(mxpt + len(listn), mxpt + len(listn) + 100):
        if (rt[0].data[val2] < 0):
            listn.append(val2)
        else:
            break
    mnpt = np.argmin(rt[0].data[mxpt:mxpt + len(listn)])
    sampl_rate = rt[0].stats.sampling_rate
    F_peak = sampl_rate / (2. * mnpt)  # Rotation frequency at peak rotation
    # measured from point of max rotation rate to the next minimum (1/2 period)
    # that's why there is a factor 2
    PDIS = max(displ[1])
    PAT = max(rotate[1])  # Peak transverse acceleration [nm/s]
    PAZ = max(ac.select(component='Z')
              [0].data)  # Peak vertical acceleration [nm/s]
    PRZ = max(rt[0].data)  # Peak vertical rotation rate [nrad/s]
    PCC = max(corrcoefs)  # Peak correlation coefficient
    TBA = baz[2]  # Theoretical backazimuth [°]
    SNT = sn_ratio(rotate[1], arriv_p, ac.select(component=compN)[
                   0].stats.sampling_rate)  # SNR for transverse acceleration
    SNZ = sn_ratio(ac.select(component='Z')[0].data, arriv_p, ac.select(
        component='Z')[0].stats.sampling_rate)  # SNR for vertical acc.
    SNR = sn_ratio(rt[0].data, arriv_p, rt[0].stats.sampling_rate)
    # Signal to noise ratio for vertical rotation rate
#    from IPython.core.debugger import Tracer; Tracer(colors="Linux")()
    dic = OrderedDict([
            ('data', OrderedDict([
                ('rotational', OrderedDict([
                    ('network', net_r),
                    ('station', sta_r),
                    ('loc', loc_r),
                    ('channel', chan1),
                    ('arclink_data_select_source', source)])),
                ('translational', OrderedDict([
                    ('network', net_s),
                    ('station', sta_s),
                    ('loc', loc_s),
                    ('channel_N', chan3),
                    ('channel_E', chan2),
                    ('channel_Z', chan4),
                    ('arclink_data_select_source', source)]))
                ])),
            ('event_id', event.resource_id.id),
            ('event_source', event_source),
            ('starttime', str(startev-180)),
            ('endtime', str(startev+3*3600)),
            ('station_latitude', str(rt[0].stats.coordinates['latitude'])),
            ('station_longitude', str(rt[0].stats.coordinates['longitude'])),
            ('event_latitude', event.preferred_origin().latitude),
            ('event_longitude', event.preferred_origin().longitude),
            ('magnitude', magnitude),
            ('depth', depth),
            ('depth_unit', 'km'),
            ('epicentral_distance', distance),
            ('epicentral_distance_unit', 'km'),
            ('peak_transverse_acceleration', PAT),
            ('peak_transverse_acceleration_unit', 'nm/s^2'),
            ('peak_transverse_displacement', PDIS),
            ('peak_transverse_displacement_unit', 'nm'),
            ('theoretical_backazimuth', TBA),
            ('theoretical_backazimuth_unit', 'degree'),
            ('estimated_backazimuth', EBA),
            ('estimated_backazimuth_unit', 'degree'),
            ('max_xcoef_for_estimated_backazimuth', max_ebaz_xcoef),
            ('frequency_at_peak_vertical_rotation_rate', F_peak),
            ('frequency_at_peak_vertical_rotation_rate_unit', 'Hz'),
            ('peak_vertical_acceleration', PAZ),
            ('peak_vertical_acceleration_unit', 'nm/s^2'),
            ('peak_vertical_rotation_rate', PRZ),
            ('peak_vertical_rotation_rate_unit', 'nrad/s'),
            ('peak_correlation_coefficient', PCC),
            ('transverse_acceleration_SNR', SNT),
            ('vertical_acceleration_SNR', SNZ),
            ('vertical_rotation_rate_SNR', SNR),
            ('phase_velocity', OrderedDict([
                ('band_1', OrderedDict([
                    ('freqmin', 0.01),
                    ('freqmax', 0.02),
                    ('freq_unit', 'Hz'),
                    ('mean_phase_vel', phasv_means[0]),
                    ('vel_std', phasv_stds[0]),
                    ('vel_unit', 'km/s')])),
                ('band_2', OrderedDict([
                    ('freqmin', 0.02),
                    ('freqmax', 0.04),
                    ('freq_unit', 'Hz'),
                    ('mean_phase_vel', phasv_means[1]),
                    ('vel_std', phasv_stds[1]),
                    ('vel_unit', 'km/s')])),
                ('band_3', OrderedDict([
                    ('freqmin', 0.04),
                    ('freqmax', 0.10),
                    ('freq_unit', 'Hz'),
                    ('mean_phase_vel', phasv_means[2]),
                    ('vel_std', phasv_stds[2]),
                    ('vel_unit', 'km/s')])),
                ('band_4', OrderedDict([
                    ('freqmin', 0.10),
                    ('freqmax', 0.20),
                    ('freq_unit', 'Hz'),
                    ('mean_phase_vel', phasv_means[3]),
                    ('vel_std', phasv_stds[3]),
                    ('vel_unit', 'km/s')])),
                ('band_5', OrderedDict([
                    ('freqmin', 0.20),
                    ('freqmax', 0.30),
                    ('freq_unit', 'Hz'),
                    ('mean_phase_vel', phasv_means[4]),
                    ('vel_std', phasv_stds[4]),
                    ('vel_unit', 'km/s')])),
                ('band_6', OrderedDict([
                    ('freqmin', 0.30),
                    ('freqmax', 0.40),
                    ('freq_unit', 'Hz'),
                    ('mean_phase_vel', phasv_means[5]),
                    ('vel_std', phasv_stds[5]),
                    ('vel_unit', 'km/s')])),
                ('band_7', OrderedDict([
                    ('freqmin', 0.40),
                    ('freqmax', 0.60),
                    ('freq_unit', 'Hz'),
                    ('mean_phase_vel', phasv_means[6]),
                    ('vel_std', phasv_stds[6]),
                    ('vel_unit', 'km/s')])),
                ('band_8', OrderedDict([
                    ('freqmin', 0.60),
                    ('freqmax', 1.0),
                    ('freq_unit', 'Hz'),
                    ('mean_phase_vel', phasv_means[7]),
                    ('vel_std', phasv_stds[7]),
                    ('vel_unit', 'km/s')]))
                ]))
            ])
    outfile = open(tag_name + '/' + tag_name[23:] + '.json', 'wt')
    json.dump(dic, outfile, indent=4)

    outfile.close()


def plotWaveformComp(event, station, link, mode, filter_e, event_source):

    """
    Compares vertical rotation rate and transversal acceleration through
    direct waveform comparison in different time windows and through cross-
    correlation analysis. It also stores some values obtained through the
    routine, like peak values (signal amplitudes, correlation coefficients)
    and signal-to-noise ratios

    :type event: :class: `~obspy.core.event.Event`
    :param event: Contains the event information.
    :type station: str
    :param station: Station from which data are fetched ('WET' or 'PFO').
    :type link: string
    :param link: Link to the Iris-xml, where event and moment tensor data are
        fetched.
    :type mode: str
    :param mode: Defines if WET data are fetched from Neries ('neries')
        or an IRIS link ('link').
    """

    tag_name = os.path.join('database/static/OUTPUT', catalog + '_' + str(event.origins[0]['time'].date) +
     'T' + str(event.origins[0]['time'].hour) + ':' + str(event.origins[0]['time'].minute) +
     '_' + str(event.magnitudes[0]['mag']) + '_' + 
     str(event.event_descriptions[0]['text'].splitlines()[0].replace(' ', '_').replace(',', '_')))
    tag_name2 = os.path.join('database/static/OUTPUT/Phase_velocities', str(event)
                             .splitlines()[0].replace('\t', '')
                             .replace(' ', '_').replace('|', '_')
                             .replace(':', '_'))
    # event information:
    latter, lonter, depth, startev, rt, ac, baz, gcdist, net_r, net_s, chan1,\
        chan2, chan3, chan4, sta_r, sta_s, loc_r, loc_s, source =\
        event_info_data(event, station, mode)
    try:
        region = event.event_descriptions[0]['text']
        reg=True
    except:
        reg=False

    Mw = float(str(event).split('\n')[0][57:60])

    if filter_e == 'True':
        a = -0.8583925  # a,b: parameters of logarithmic curve
        b = 0.76685464
        Mag = b * np.log(0.001 * baz[0]) + a  # required magnitude for process
        if Mag > Mw:
            print 'The event does probably not allow for velocity estimates as'\
                'distance is too large or magnitude too small'
            os.rmdir(tag_name)
            return
    
    MomentTensor, Magnitude, Region = Get_MomentTensor_Magnitude(link)
    # first page: map with event location & information
    print 'Plotting map with station, event and great circle...'
    if is_local(baz) == 'local':
        plt.figure(figsize=(18, 9))

        plt.subplot2grid((4, 9), (0, 4), colspan=5, rowspan=4)

        # conic map plot
        map = Basemap(projection='lcc', lat_0=(rt[0].stats.coordinates.
                      latitude + latter) / 2, lon_0=(rt[0].stats.coordinates.
                      longitude + lonter) / 2, resolution='i', width=3000000,
                      height=2000000)
        map.drawparallels(np.arange(0., 90, 5.), labels=[1, 0, 0, 1])
        map.drawmeridians(np.arange(0., 360., 5.), labels=[1, 0, 0, 1])
        map.drawstates(linewidth=0.25)
    elif is_local(baz) == 'non-local' and baz[0] <= 13000000:
        plt.figure(figsize=(18, 9))
        plt.subplot2grid((4, 9), (0, 4), colspan=5, rowspan=4)
        # globe plot
        map = Basemap(projection='ortho', lat_0=(rt[0].stats.coordinates.
                      latitude + latter) / 2, lon_0=(rt[0].stats.coordinates.
                      longitude + lonter) / 2, resolution='l')
        map.drawmeridians(np.arange(0, 360, 30))
        map.drawparallels(np.arange(-90, 90, 30))
    elif is_local(baz) == 'non-local' and baz[0] > 13000000:
        plt.figure(figsize=(18, 9))
        plt.subplot2grid((4, 9), (0, 4), colspan=5, rowspan=4)
        # If the great circle between the station and event is crossing the
        # 180° meridian in the pacific and the stations are far apart the map
        # has to be re-centered, otherwise we would see the wrong side of the
        # globe.
        if abs(rt[0].stats.coordinates.longitude - lonter) > 180:
            lon0 = 180 + (rt[0].stats.coordinates.longitude + lonter) / 2
        else:
            lon0 = (rt[0].stats.coordinates.longitude + lonter) / 2
        map = Basemap(projection='moll', lon_0=lon0, resolution='l')
        map.drawmeridians(np.arange(0, 360, 30))
        map.drawparallels(np.arange(-90, 90, 30))
    else:
        plt.figure(figsize=(26, 13))
        plt.title('Event: %s %s\n \n '
                  % (startev.date, startev.time), fontsize=24, fontweight='bold')
        plt.subplot2grid((4, 9), (0, 4), colspan=5, rowspan=4)
        # conic map plot
        map = Basemap(projection='lcc', lat_0=(rt[0].stats.coordinates.
                      latitude + latter) / 2, lon_0=(rt[0].stats.coordinates.
                      longitude + lonter) / 2, resolution='i', width=600000,
                      height=400000)
        map.drawparallels(np.arange(0., 90, 2.), labels=[1, 0, 0, 1])
        map.drawmeridians(np.arange(0., 360., 2.), labels=[1, 0, 0, 1])
        map.drawrivers(linewidth=0.25, color='b')
        map.drawstates(linewidth=0.25)

    map.drawcoastlines(linewidth=0.25)
    map.drawcountries(linewidth=0.25)
    map.fillcontinents(color='coral', lake_color='lightblue')
    map.drawmapboundary(fill_color='lightblue')
    if is_local(baz) == 'local' or is_local(baz) == 'close':
        map.drawlsmask(land_color='coral', ocean_color='lightblue', lakes=True)
        map.drawcountries(linewidth=0.6)
    map.drawgreatcircle(lonter, latter, rt[0].stats.coordinates.longitude,
                        rt[0].stats.coordinates.latitude, linewidth=3,
                        color='yellow')

    # Add beachballs for the event and station triangle
    if station == 'WET':
        x, y = map(lonter, latter)
        statlon, statlat = map(rt[0].stats.coordinates.longitude, rt[0].stats.
                               coordinates.latitude)

        if is_local(baz) == 'non-local':
            map.scatter(statlon, statlat, 200, color="b", marker="v",
                        edgecolor="k", zorder=100)
            plt.text(statlon + 200000, statlat, 'WET', va="top",
                     family="monospace", weight="bold", zorder=101,
                     color='k', backgroundcolor='white')
            map.scatter(x, y, 200, color="b", marker="*",
                        edgecolor="k", zorder=100)

            if mode == 'link':
                plt.subplot2grid((4, 9), (1, 0), colspan=2)
                plt.title(u'Event: %s \n %s \n \n' % (startev, Region),
                          fontsize=20, weight='bold')
                ax = plt.gca()
                ax.axis('equal')
                ax.axis('off')
                b = Beach(MomentTensor, xy=(0.5, 0.5), facecolor='blue',
                          width=0.5, linewidth=1, alpha=1.0)
                b.set_zorder(200)
                ax.add_collection(b)

                plt.subplot2grid((4, 9), (2, 0), colspan=2)
                plt.title(u'\nMagnitude: %s \n\nDistance: %.2f [km] - %.2f [°]'
                          '\n\nDepth: %.2f [km]'
                          % (str(event).split('\n')[0][57:64], 0.001 * baz[0],
                             0.001 * baz[0] / 111.11, depth), fontsize=18,
                          fontweight='bold')
                ax = plt.gca()
                ax.axis('off')

                plt.subplot2grid((4, 9), (3, 0), colspan=2)
                plt.title(u'Event Information: \n Global Centroid-Moment-Tensor '
                          'Catalog (GCMT)'
                          '\n\n Processing Date:\n'+str(UTCDateTime().date),
                          fontsize=14)
                ax = plt.gca()
                ax.axis('off')

            else:
                map.scatter(x, y, 200, color="b", marker="*", edgecolor="k",
                            zorder=100)

                plt.subplot2grid((4, 9), (1, 0), colspan=2)
                plt.title(u'Event: %s %s\n' % (startev.date, startev.time), fontsize=20,
                          weight='bold')
                ax = plt.gca()
                ax.axis('equal')
                ax.axis('off')

                plt.subplot2grid((4, 9), (2, 0), colspan=2)
                if reg==True:
                    plt.title(u'\n\nRegion: %s \n\nMagnitude: %s \n\nDistance: %.2f [km] - %.2f [°]\n\nDepth: %.2f [km]'
                              % (region, str(event).split('\n')[0][57:64], 0.001 * baz[0],
                                 0.001 * baz[0] / 111.11, depth),
                              fontsize=18, fontweight='bold')
                else:
                    plt.title(u'\n\nMagnitude: %s \n\nDistance: %.2f [km] - %.2f [°]\n\nDepth: %.2f [km]'
                              % (str(event).split('\n')[0][57:64], 0.001 * baz[0],
                                 0.001 * baz[0] / 111.11, depth),
                              fontsize=18, fontweight='bold')
                ax = plt.gca()
                ax.axis('off')

                plt.subplot2grid((4, 9), (3, 0), colspan=2)
                plt.title(u'Event Information: \n Global Centroid-Moment-Tensor '
                          'Catalog (GCMT)'
                          '\n\n Processing Date:\n'+str(UTCDateTime().date),
                          fontsize=14)
                ax = plt.gca()
                ax.axis('off')

        else:
            map.scatter(statlon, statlat, 200, color="b", marker="v",
                        edgecolor="k", zorder=100)
            plt.text(statlon + 35000, statlat, 'WET', fontsize=12, va="top",
                     family="monospace", weight="bold", zorder=101,
                     color='k', backgroundcolor='white')
            map.scatter(x, y, 300, color="b", marker="*",
                        edgecolor="k", zorder=100)

            if mode == 'link':
                plt.subplot2grid((4, 9), (1, 0), colspan=2)
                plt.title(u'Event: %s \n %s \n \n' % (startev, Region),
                          fontsize=20, weight='bold')
                ax = plt.gca()
                ax.axis('equal')
                ax.axis('off')
                b = Beach(MomentTensor, xy=(0.5, 0.5), facecolor='blue',
                          width=0.5, linewidth=1, alpha=1.0)
                b.set_zorder(200)
                ax.add_collection(b)

                plt.subplot2grid((4, 9), (2, 0), colspan=2)

                plt.title(u'\nMagnitude: %s \n\nDistance: %.2f [km] - %.2f [°]'
                          '\n\nDepth: %.2f [km]'
                          % (str(event).split('\n')[0][57:64], 0.001 * baz[0],
                             0.001 * baz[0] / 111.11, depth), fontsize=18,
                          fontweight='bold')
                ax = plt.gca()
                ax.axis('off')
                plt.subplot2grid((4, 9), (3, 0), colspan=2)
                plt.title(u'Event Information: \n Global Centroid-Moment-Tensor '
                          'Catalog (GCMT)'
                          '\n\n Processing Date:\n'+str(UTCDateTime().date),
                          fontsize=14)
                ax = plt.gca()
                ax.axis('off')

            else:
                plt.subplot2grid((4, 9), (1, 0), colspan=2)
                plt.title(u'Event: %s %s\n' % (startev.date, startev.time), fontsize=20,
                          weight='bold')
                ax = plt.gca()
                ax.axis('equal')
                ax.axis('off')

                plt.subplot2grid((4, 9), (2, 0), colspan=2)
                if reg==True:
                    plt.title(u'\n\nRegion: %s \n\nMagnitude: %s \n\nDistance: %.2f [km] - %.2f [°]\n\nDepth: %.2f [km]'
                              % (region, str(event).split('\n')[0][57:64], 0.001 * baz[0],
                                 0.001 * baz[0] / 111.11, depth),
                              fontsize=18, fontweight='bold')
                else:
                    plt.title(u'\nMagnitude: %s \n\nDistance: %.2f [km] - %.2f [°]'
                              '\n\nDepth: %.2f [km]'
                              % (str(event).split('\n')[0][57:64], 0.001 * baz[0],
                                 0.001 * baz[0] / 111.11, depth), fontsize=18,
                              fontweight='bold')
                ax = plt.gca()
                ax.axis('off')
                plt.subplot2grid((4, 9), (3, 0), colspan=2)
                plt.title(u'Event Information: \n Global Centroid-Moment-Tensor '
                          'Catalog (GCMT)'
                          '\n\n Processing Date:\n'+str(UTCDateTime().date),
                          fontsize=14)
                ax = plt.gca()
                ax.axis('off')

    else:
        x, y = map(lonter, latter)
        statlon, statlat = map(rt[0].stats.coordinates.longitude, rt[0].stats.
                               coordinates.latitude)

        if is_local(baz) == 'local':
            map.scatter(x, y, 600, color="b", marker="*", edgecolor="k",
                        zorder=100)
            map.scatter(statlon, statlat, 700, color="b", marker="v",
                        edgecolor="k", zorder=100)
            plt.text(statlon + 27000, statlat, 'PFO', fontsize=18, va="top",
                     family="monospace", weight="bold", zorder=101,
                     color='k', backgroundcolor='white')
        elif is_local(baz) == 'non-local':
            map.scatter(x, y, 200, color="b", marker="*", edgecolor="k",
                        zorder=100)
            map.scatter(statlon, statlat, 300, color="b", marker="v",
                        edgecolor="k", zorder=100)
            plt.text(statlon + 200000, statlat, 'PFO', va="top",
                     family="monospace", weight="bold", zorder=101,
                     color='k', backgroundcolor='white')
        else:
            map.scatter(x, y, 250, color="b", marker="*", edgecolor="k",
                        zorder=100)
            map.scatter(statlon, statlat, 400, color="b", marker="v",
                        edgecolor="k", zorder=100)
            plt.text(statlon + 6000, statlat + 1000, 'PFO', va="top",
                     family="monospace", weight="bold", zorder=101,
                     color='k', backgroundcolor='white', fontsize=18)

        plt.subplot2grid((4, 9), (1, 0), colspan=2)
        plt.title(u'Event: %s %s\n' % (startev.date, startev.time), fontsize=20, weight='bold')
        ax = plt.gca()
        ax.axis('equal')
        ax.axis('off')

        plt.subplot2grid((4, 9), (3, 0), colspan=2)
        if reg==True:
            plt.title(u'\n\nRegion: %s \n\nMagnitude: %s \n\nDistance: %.2f [km] - %.2f [°]\n\nDepth: %.2f [km]'
                      % (region, str(event).split('\n')[0][57:64], 0.001 * baz[0],
                         0.001 * baz[0] / 111.11, depth),
                      fontsize=18, fontweight='bold')
        else:
            plt.title(u'\nMagnitude: %s \n\nDistance: %.2f [km] - %.2f [°]'
                      '\n\nDepth: %.2f [km]'
                      % (str(event).split('\n')[0][57:63], 0.001 * baz[0],
                         0.001 * baz[0] / 111.11, depth), fontsize=18,
                      fontweight='bold')
        ax = plt.gca()
        ax.axis('off')

    print 'Done.'
    plt.savefig(tag_name + '/' + tag_name[23:] + '_page_1.png')
    plt.close()

    # Preprocesing of rotational and translational signal

    # check if event is local
    # + resample signal accordingly
    # + Different cut-off frequency for the lowpass filter,
    #   depending on the epicentral distance
    rt, ac, rt_pcoda, ac_pcoda, sec, cutoff, cutoff_pc =\
        resample(is_local(baz), baz, rt, ac)

    print 'Removing instrument response...'
    rt, ac, rt_pcoda, ac_pcoda, displacement = remove_instr_resp(rt, ac, rt_pcoda, ac_pcoda,
                                                   station, startev)

    print 'Filtering and rotating...'
    rotate, pcod_rotate, pcoda_rotate, frtp, facp, frotate, cop_rt, rt_band1,\
        rt_band2, rt_band3, rt_band4, rt_band5, rt_band6, rt_band7, rt_band8,\
        rotate_band1, rotate_band2, rotate_band3, rotate_band4, rotate_band5,\
        rotate_band6, rotate_band7, rotate_band8 =\
        filter_and_rotate(ac, rt, baz, rt_pcoda, ac_pcoda, cutoff, cutoff_pc,
                          station, is_local(baz))

    # process displacement the same way as acceleration
    compE, compN = station_components(station)
    displacement.filter('lowpass', freq=cutoff, corners=2, zerophase=True)
    displacement.filter('highpass', freq=0.005, corners=2, zerophase=True)
    displacement.taper(max_percentage=0.05)
    displacement_rot = rotate_ne_rt(displacement.select(component=compN)[
                          0].data, displacement.select(component=compE)[0].data, baz[2])
    # plt.plot(displacement_rot[1]/np.abs(max(displacement_rot[1])),'k')
    # plt.plot(rotate[1]/np.abs(max(rotate[1])),'g')
    # plt.show()

    print 'Getting arrival times...'
    init_sec = startev - ac[0].stats.starttime
    # When the event starts in the fetched data

    arriv_p, arriv_s = ps_arrival_times(baz[0], depth, init_sec)

    min_pw, max_pw, min_sw, max_sw, min_lwi, max_lwi, min_lwf, max_lwf =\
        time_windows(baz, arriv_p, arriv_s, init_sec, is_local(baz))


    # test if ring laser signal is flipped, which happened a few times:
    st_rt = rt[0].data[int(rt[0].stats.sampling_rate*min_lwi):int(rt[0].stats.sampling_rate*max_lwi)]
    st_rt = st_rt/np.max(np.abs(st_rt))
    st_ac = rotate[1][int(rt[0].stats.sampling_rate*min_lwi):int(rt[0].stats.sampling_rate*max_lwi)]
    st_ac = st_ac/np.max(np.abs(st_ac))
    testCC = xcorr(st_rt, st_ac, shift_len=100, full_xcorr=True) 

    #if testCC[1] < -0.8:  # if the CC is highly negative, the signal is flipped
     #   rt[0].data = rt[0].data*(-1)  # if flipped, re-flip it!

    # Waveform comparison plot
    print 'Waveform comparison plot...'
    time = rt[0].stats.delta * np.arange(0, len(rt[0].data))  # Time in seconds
    # phase velocity determination
    # TODO minimize misfit or waterlevel method
    # TAPER
    rt.taper(max_percentage=0.05)
    fact1 = 2 * max(rt[0].data)  # Factor to displace the rotation rate

    c1 = .5 * max(abs(rotate[1])) / max(abs(rt[0].data))  # Vel in m/s

    plt.figure(figsize=(18, 9))
    plt.subplot2grid((6, 5), (2, 0), colspan=5, rowspan=2)
    plt.plot(time, rt[0].data, color='r', label=ur'Rotation rate')

    plt.plot(time, (0.5 / c1) *
             rotate[1] + fact1, color='k', label=ur'Transversal acceleration')

    plt.xlabel(ur'Time [s]', fontweight='bold', fontsize=13)
    plt.ylabel(
        ur'$\dot{\mathbf{\Omega}}_\mathbf{z}$ [nrad/s] - a$_\mathbf{T}$/2c'
        '[1/s]', fontweight='bold', fontsize=13)
    plt.xlim(0, rt[0].stats.delta * len(rt[0].data))
    plt.ylim(min(rt[0].data), fact1 + max((1. / (2. * c1)) * rotate[1]))
    box_yposition = ((fact1 + max((1. / (2. * c1)) * rotate[1]))
                     - abs(min(rt[0].data)))/2  # box is in middle of figure
    if is_local(baz) == 'non-local':  # gap between annotation and vline
        xgap = 50
    else:
        xgap = 15
    xlim1 = rt[0].stats.delta * len(rt[0].data)  # needed for p-coda later
    bbox_props = dict(boxstyle="square, pad=0.3", fc='white')
    plt.axvline(x=min_pw, linewidth=1)
    plt.annotate('1', xy=(min_pw+xgap, box_yposition), fontsize=14,
                 fontweight='bold', bbox=bbox_props)
    plt.axvline(x=min_sw, linewidth=1)
    plt.annotate('2', xy=(min_sw+xgap, box_yposition), fontsize=14,
                 fontweight='bold', bbox=bbox_props)
    plt.axvline(x=min_lwi, linewidth=1)
    plt.annotate('3', xy=(min_lwi+xgap, box_yposition), fontsize=14,
                 fontweight='bold', bbox=bbox_props)
    plt.axvline(x=min_lwf, linewidth=1)
    plt.annotate('4', xy=(min_lwf+xgap, box_yposition), fontsize=14,
                 fontweight='bold', bbox=bbox_props)
    plt.title(ur'Ring laser and broadband seismometer recordings. Event: %s %s'
              % (startev.date, startev.time))
    plt.grid(True)
    plt.legend(loc=7,shadow=True)
    # P coda
    plt.subplot2grid((6, 5), (0, 0), colspan=2)
    cp = 0.5 * max(abs(pcod_rotate[1][cop_rt[0].stats.sampling_rate * min_pw:
                                      cop_rt[0].stats.sampling_rate * max_pw])
                   ) / max(abs(cop_rt[0].data[cop_rt[0].stats.sampling_rate *
                                              min_pw:cop_rt[0].stats.
                                              sampling_rate * max_pw]))
    minamp1_pcod = min((0.5 / cp) * pcod_rotate[1][rt[0].stats.sampling_rate
                                                   * min_pw:rt[0].stats.
                                                   sampling_rate * max_pw])
    minamp2_pcod = min(cop_rt[0].data[rt[0].stats.sampling_rate * min_pw:rt[0]
                                      .stats.sampling_rate * max_pw])
    maxamp1_pcod = max((0.5 / cp) * pcod_rotate[1][rt[0].stats.sampling_rate
                                                   * min_pw:rt[0].stats.
                                                   sampling_rate * max_pw])
    maxamp2_pcod = max(cop_rt[0].data[rt[0].stats.sampling_rate * min_pw:rt[0]
                                      .stats.sampling_rate * max_pw])
    plt.plot(time, cop_rt[0].data, color='r')
    plt.plot(time, (0.5 / cp) * pcod_rotate[1], color='k')
    plt.xlim(min_pw, max_pw)
    plt.ylim(
        min([minamp1_pcod, minamp2_pcod]),
        max([maxamp1_pcod, maxamp2_pcod]))
    plt.xlabel(ur'Time [s]', fontweight='bold', fontsize=11)
    plt.ylabel(
        ur'$\dot{\mathbf{\Omega}}_\mathbf{z}$ [nrad/s] - a$_\mathbf{T}$/2c'
        '[1/s]', fontweight='bold', fontsize=11)
    plt.title(u'1: P-coda (Highpass, cut-off: %.1f Hz)' % cutoff_pc)
    plt.grid(True)

    # S wave
    plt.subplot2grid((6, 5), (0, 3), colspan=2)
    # colorscale??
    cs = 0.5 * max(abs(rotate[1][rt[0].stats.sampling_rate * min_sw:rt[0].
                                 stats.sampling_rate * max_sw])) /\
        max(abs(rt[0].data[rt[0].stats.sampling_rate * min_sw:rt[0].stats.
                           sampling_rate * max_sw]))
    minamp1_s = min((0.5 / cs) * rotate[1][rt[0].stats.sampling_rate
                                           * min_sw:rt[0].stats.sampling_rate
                                           * max_sw])
    minamp2_s = min(rt[0].data[rt[0].stats.sampling_rate * min_sw:rt[0].stats.
                               sampling_rate * max_sw])
    maxamp1_s = max((0.5 / cs) * rotate[1][rt[0].stats.sampling_rate
                                           * min_sw:rt[0].stats.sampling_rate
                                           * max_sw])
    maxamp2_s = max(rt[0].data[rt[0].stats.sampling_rate * min_sw:rt[0].stats.
                               sampling_rate * max_sw])
    plt.plot(time, rt[0].data, color='r')
    plt.plot(time, (0.5 / cs) * rotate[1], color='k')
    plt.xlim(min_sw, max_sw)
    plt.ylim(
        min([minamp1_s, minamp2_s]),
        max([maxamp1_s, maxamp2_s]))
    plt.xlabel(ur'Time [s]', fontweight='bold', fontsize=11)
    plt.ylabel(
        ur'$\dot{\mathbf{\Omega}}_\mathbf{z}$ [nrad/s] - a$_\mathbf{T}$/2c'
        '[1/s]', fontweight='bold', fontsize=11)
    plt.title(u'2: S-wave (Lowpass, cut-off: %s Hz)' % (cutoff))
    plt.grid(True)

    # surface waves
    plt.subplot2grid((6, 5), (5, 0), colspan=2)
    cl1 = 0.5 * max(abs(rotate[1][rt[0].stats.sampling_rate * min_lwi:rt[0].
                                  stats.sampling_rate * max_lwi])) /\
        max(abs(rt[0].data[rt[0].stats.sampling_rate * min_lwi:rt[0].
            stats.sampling_rate * max_lwi]))
    minamp1_surf = min((0.5 / cl1) * rotate[1][rt[0].stats.sampling_rate *
                                               min_lwi:rt[0].stats.
                                               sampling_rate * max_lwi])
    minamp2_surf = min(rt[0].data[rt[0].stats.sampling_rate * min_lwi:rt[0].
                                  stats.sampling_rate * max_lwi])
    maxamp1_surf = max((0.5 / cl1) * rotate[1][rt[0].stats.sampling_rate
                                               * min_lwi:rt[0].stats.
                                               sampling_rate * max_lwi])
    maxamp2_surf = max(rt[0].data[rt[0].stats.sampling_rate *
                                  min_lwi:rt[0].stats.sampling_rate * max_lwi])
    plt.plot(time, rt[0].data, color='r')
    plt.plot(time, (0.5 / cl1) * rotate[1], color='k')
    plt.xlim(min_lwi, max_lwi)
    plt.ylim(
        min([minamp1_surf, minamp2_surf]),
        max([maxamp1_surf, maxamp2_surf]))
    plt.xlabel(ur'Time [s]', fontweight='bold', fontsize=11)
    plt.ylabel(
        ur'$\dot{\mathbf{\Omega}}_\mathbf{z}$ [rad/s] - a$_\mathbf{T}$/2c'
        '[1/s]', fontweight='bold', fontsize=11)
    plt.title(ur'3: Initial surface waves (Lowpass, cut-off: %s Hz)'
              % (cutoff))
    plt.grid(True)

    # later surface waves
    plt.subplot2grid((6, 5), (5, 3), colspan=2)
    cl2 = 0.5 * max(abs(rotate[1][rt[0].stats.sampling_rate * min_lwf:rt[0].
                                  stats.sampling_rate * max_lwf])) /\
        max(abs(rt[0].data[rt[0].stats.sampling_rate * min_lwf:rt[0].
                stats.sampling_rate * max_lwf]))
    minamp1_lat = min((0.5 / cl2) * rotate[1][rt[0].stats.sampling_rate
                                              * min_lwf:rt[0].stats.
                                              sampling_rate * max_lwf])
    minamp2_lat = min(rt[0].data[rt[0].stats.sampling_rate * min_lwf:rt[0].
                                 stats.sampling_rate * max_lwf])
    maxamp1_lat = max((0.5 / cl2) * rotate[1][rt[0].stats.sampling_rate
                                              * min_lwf:rt[0].stats.
                                              sampling_rate * max_lwf])
    maxamp2_lat = max(rt[0].data[rt[0].stats.sampling_rate * min_lwf:rt[0].
                                 stats.sampling_rate * max_lwf])
    plt.plot(time, rt[0].data, color='r')
    plt.plot(time, (0.5 / cl2) * rotate[1], color='k')
    plt.xlim(min_lwf, max_lwf)
    plt.ylim(
        min([minamp1_lat, minamp2_lat]),
        max([maxamp1_lat, maxamp2_lat]))
    plt.xlabel(ur'Time [s]', fontsize=11, fontweight='bold')
    plt.ylabel(
        ur'$\dot{\mathbf{\Omega}}_\mathbf{z}$ [rad/s] - a$_\mathbf{T}$/2c'
        '[1/s]', fontsize=11, fontweight='bold')
    plt.title(ur'4: Later surface waves (Lowpass, cut-off: %s Hz)'
              % (cutoff))
    plt.grid(True)
    print 'Done.'
    plt.savefig(tag_name + '/' + tag_name[23:] + '_page_2.png')
    plt.close()

    # Cross-correlation analysis
    print 'Obtaining zero-lag correlation coefficients for theoretical' \
        ' backazimuth...'
    corrcoefs, thres = Get_corrcoefs(rt[0], rt[0].data, ac, rotate[1], sec,
                                     station)
    # calculate coefs for different frequency bands
    corrcoefs_band1, thres1 = Get_corrcoefs(rt_band1[0], rt_band1[0].data, ac,
                                            rotate_band1[1], 200, station)
    corrcoefs_band2, thres2 = Get_corrcoefs(rt_band2[0], rt_band2[0].data, ac,
                                            rotate_band2[1], 100, station)
    corrcoefs_band3, thres3 = Get_corrcoefs(rt_band3[0], rt_band3[0].data, ac,
                                            rotate_band3[1], 50, station)
    corrcoefs_band4, thres4 = Get_corrcoefs(rt_band4[0], rt_band4[0].data, ac,
                                            rotate_band4[1], 20, station)
    corrcoefs_band5, thres5 = Get_corrcoefs(rt_band5[0], rt_band5[0].data, ac,
                                            rotate_band5[1], 12, station)
    corrcoefs_band6, thres6 = Get_corrcoefs(rt_band6[0], rt_band6[0].data, ac,
                                            rotate_band6[1], 10, station)
    corrcoefs_band7, thres7 = Get_corrcoefs(rt_band7[0], rt_band7[0].data, ac,
                                            rotate_band7[1], 8, station)
    corrcoefs_band8, thres8 = Get_corrcoefs(rt_band8[0], rt_band8[0].data, ac,
                                            rotate_band8[1], 6, station)

    # zero-lag correlation coefficients for range of backazimuths
    print 'Backazimuth analysis...'
    corrbaz, maxcorr, backas, max_coefs_10deg = \
        backas_analysis(rt[0], rt[0].data, ac, sec, corrcoefs, None, station)

    X, Y = np.meshgrid(np.arange(0, sec * len(corrcoefs), sec), backas)

    # Estimating backazimuth
    print 'Estimating backazimuth...'
    corrsum, backas2, max_ebaz_xcoef, best_ebaz = backas_est(rt, ac, min_sw, max_lwf, station)

    # calculate phase veloc. for windows where corrcoef is good enough (.75)
    # optional TODO: threshold value of corrcoef as a parameter
    # TODO phase velocity estimation by misfit minimization
    print 'Phase velocities...'

    # calculate startindex for phase velocity calculation in frequency bands:
    # starts at the beginning of surface wave arrivals as body waves are
    # not appropriate
    ind_surf = int(min_lwi/sec)

    ind_band = False  # indicator that we are dealing with bands 1-8
    phasv, EBA = phase_vel(rt, sec, corrcoefs, rotate, corrsum, backas2,
                           ind_band, ind_surf)
    # calculates phase velocities for different frequency bands -
    # -> put EBA_bandx here instead of EBA and backas2_bandx,...
    # for frequency dependent EBA
    ind_band = True  # indicator that we are dealing with bands 1-8

    phasv_band1, EBA1 = phase_vel(rt_band1, 200, corrcoefs_band1, rotate_band1,
                                  corrsum, backas2, ind_band, ind_surf)
    phasv_band2, EBA2 = phase_vel(rt_band2, 100, corrcoefs_band2, rotate_band2,
                                  corrsum, backas2, ind_band, ind_surf)
    phasv_band3, EBA3 = phase_vel(rt_band3, 50, corrcoefs_band3, rotate_band3,
                                  corrsum, backas2, ind_band, ind_surf)
    phasv_band4, EBA4 = phase_vel(rt_band4, 20, corrcoefs_band4, rotate_band4,
                                  corrsum, backas2, ind_band, ind_surf)
    phasv_band5, EBA5 = phase_vel(rt_band5, 12, corrcoefs_band5, rotate_band5,
                                  corrsum, backas2, ind_band, ind_surf)
    phasv_band6, EBA6 = phase_vel(rt_band6, 10, corrcoefs_band6, rotate_band6,
                                  corrsum, backas2, ind_band, ind_surf)
    phasv_band7, EBA7 = phase_vel(rt_band7, 8, corrcoefs_band7, rotate_band7,
                                  corrsum, backas2, ind_band, ind_surf)
    phasv_band8, EBA8 = phase_vel(rt_band8, 6, corrcoefs_band8, rotate_band8,
                                  corrsum, backas2, ind_band, ind_surf)
    ind_band = False
    phasv_band1 = phasv_band1[~np.isnan(phasv_band1)]  # filters out all NaN
    phasv_band2 = phasv_band2[~np.isnan(phasv_band2)]
    phasv_band3 = phasv_band3[~np.isnan(phasv_band3)]
    phasv_band4 = phasv_band4[~np.isnan(phasv_band4)]
    phasv_band5 = phasv_band5[~np.isnan(phasv_band5)]
    phasv_band6 = phasv_band6[~np.isnan(phasv_band6)]
    phasv_band7 = phasv_band7[~np.isnan(phasv_band7)]
    phasv_band8 = phasv_band8[~np.isnan(phasv_band8)]

    plt.figure()
    ax = plt.subplot(111)
    lens = len(phasv_band1) + len(phasv_band2) + len(phasv_band3) + \
        len(phasv_band4) + len(phasv_band5) + len(phasv_band6) + \
        len(phasv_band7) + len(phasv_band8)
    arr = np.zeros((lens, 4))
    leni = 0
    for phasv_band in [phasv_band1, phasv_band2, phasv_band3, phasv_band4,
                       phasv_band5, phasv_band6, phasv_band7, phasv_band8]:
        arr[leni:leni + len(phasv_band), 0] = phasv_band
        leni += len(phasv_band)
    arr[:, 2] = [depth]*lens  # source depth
    arr[:, 3] = [0.001 * baz[0]]*lens  # epicentral distance

    l1 = len(phasv_band1)
    l2 = len(phasv_band2)
    l3 = len(phasv_band3)
    l4 = len(phasv_band4)
    l5 = len(phasv_band5)
    l6 = len(phasv_band6)
    l7 = len(phasv_band7)
    l8 = len(phasv_band8)

    arr[:l1, 1] = [1. / 0.015] * len(phasv_band1)
    arr[l1:l1+l2, 1] = [1. / 0.03] * len(phasv_band2)
    arr[l1+l2:l1+l2+l3, 1] = [1./0.07]*len(phasv_band3)
    arr[l1+l2+l3:l1+l2+l3+l4, 1] = [1./0.15]*len(phasv_band4)
    arr[l1+l2+l3+l4:l1+l2+l3+l4+l5, 1] = [1./0.25]*len(phasv_band5)
    arr[l1+l2+l3+l4+l5:l1+l2+l3+l4+l5+l6, 1] = [1./0.35]*len(phasv_band6)
    arr[l1+l2+l3+l4+l5+l6:l1+l2+l3+l4+l5+l6+l7, 1] = [1./0.5]*len(phasv_band7)
    arr[l1+l2+l3+l4+l5+l6+l7:l1+l2+l3+l4+l5+l6+l7+l8, 1] =\
        [1./0.8]*len(phasv_band8)
    np.save(tag_name2, arr)
    phasv_means = []
    phasv_stds = []
    for ksk in [phasv_band1, phasv_band2, phasv_band3, phasv_band4,
                phasv_band5, phasv_band6, phasv_band7, phasv_band8]:
        phasv_means.append(np.mean(ksk))
        # calculates means for single f-bands and appends them to a list
        phasv_stds.append(np.std(ksk))

    print 'Cross-correlation and phase velocity figures...'

    plt.figure(figsize=(18, 9))
    plt.subplot2grid((4, 26), (0, 0), colspan=25)
    plt.plot(time, rt[0].data, color='r', label=ur'Rotation rate')
    plt.plot(time, (1. / (2. * c1)) *
             rotate[1] + fact1, color='k', label=ur'Transversal acceleration')
    plt.ylabel(
        ur'$\dot{\mathbf{\Omega}}_\mathbf{z}$ [nrad/s] - a$_\mathbf{T}$/2c '
        '[1/s]', fontsize=10, fontweight='bold')
    plt.xlim(0, rt[0].stats.delta * len(rt[0].data))
    plt.ylim(min(rt[0].data), fact1 + max((1. / (2. * c1)) * rotate[1]))
    plt.title(ur'Cross-correlation for $\dot\Omega_z$ and a$_T$ in %s seconds '
              'time windows (lowpass, cutoff: %s Hz). Event: %s %s'
              % (sec, cutoff, startev.date, startev.time))
    plt.grid(True)
    plt.legend(loc=7,shadow=True)
    plt.subplot2grid((4, 26), (1, 0), colspan=25)
    plt.scatter(np.arange(0, sec * len(phasv), sec), phasv,
                c=corrcoefs, vmin=0.75, vmax=1, s=35,
                cmap=plt.cm.autumn_r)
    plt.ylabel(ur'Phase velocity [km/s]', fontsize=10, fontweight='bold')
    plt.xlim(0, rt[0].stats.delta * len(rt[0].data))
    plt.ylim(0, 16)
    plt.grid(True)
    fig = plt.subplot2grid((4, 26), (1, 25))
    cmap = mpl.cm.autumn_r
    norm = mpl.colors.Normalize(vmin=0.75, vmax=1)
    cb1 = mpl.colorbar.ColorbarBase(
        fig, cmap=cmap, norm=norm, orientation='vertical',
        ticks=[0.75, 0.8, 0.85, 0.9, 0.95, 1.0])
    cb1.set_label(ur'X-corr. coeff.', fontweight='bold')

    plt.subplot2grid((4, 26), (2, 0), colspan=25)
    plt.plot(np.arange(0, sec * len(corrcoefs), sec), max_coefs_10deg, 'ro-', label='Max. CC for est. BAz', linewidth=1)
    plt.plot(np.arange(0, sec * len(corrcoefs), sec), corrcoefs, 'ko-', label='CC for theo. BAz', linewidth=1)
    plt.plot(np.arange(0, sec * len(corrcoefs) + 1, sec), thres, '--r', lw=2)
    plt.ylabel(ur'X-corr. coeff.', fontsize=10, fontweight='bold')
    plt.text(time[len(time) - 1] + 50, 0.71, ur'0.75', color='red')
    plt.xlim(0, rt[0].stats.delta * len(rt[0].data))
    plt.ylim(0, 1)
    plt.legend(loc=4, shadow=True)
    plt.grid(True)

    plt.subplot2grid((4, 26), (3, 0), colspan=25)
    teobaz = baz[2] * np.ones(len(corrcoefs) + 1)
    plt.pcolor(X, Y, corrbaz, cmap=plt.cm.RdYlGn_r)
    plt.plot(np.arange(0, sec * len(corrcoefs) + 1, sec), teobaz, '--r', lw=2)
    plt.plot(np.arange(0, sec * len(corrcoefs), sec), maxcorr, '.k')
    plt.text(1000, baz[2], str(baz[2])[0:5] + ur'°',
             bbox={'facecolor': 'black', 'alpha': 0.8}, color='r')
    if EBA == backas2[np.asarray(corrsum).argmax()]:
        obsbaz = EBA * np.ones(len(corrcoefs) + 1)
        plt.text(400, EBA, str(EBA)[0:5] + ur'°',
                 bbox={'facecolor': 'black', 'alpha': 0.8}, color='y')
        plt.plot(np.arange(0, sec * len(corrcoefs) + 1, sec),
                 obsbaz, '--y', lw=2)
    plt.xlim(0, rt[0].stats.delta * len(rt[0].data))
    plt.xlabel(ur'Time [s]', fontweight='bold')
    plt.ylabel(ur'BAz [°]', fontsize=10, fontweight='bold')
    plt.ylim([0, 360])
    plt.yticks([0,60,120,180,240,300,360])
    plt.grid(True)
    fig = plt.subplot2grid((4, 26), (3, 25))
    norm = mpl.colors.Normalize(vmin=-1, vmax=1)
    cb1 = mpl.colorbar.ColorbarBase(
        fig, cmap=plt.cm.RdYlGn_r, norm=norm, orientation='vertical')
    cb1.set_label(ur'X-corr. coeff.', fontweight='bold')
    cb1.set_ticks([-1.0,-0.75,-0.5,-0.25,0.0,0.25,0.5,0.75,1.0])
    print 'Done.'
    plt.savefig(tag_name + '/' + tag_name[23:] + '_page_3.png')
    plt.close()

    print 'Analyzing rotations in the P-coda...'
    print 'Obtaining zero-lag correlation coefficients for theoretical' \
        ' backazimuth...'

    if is_local(baz)=="non-local":
        sec_p = 5
    else:
        sec_p = 2

    corrcoefs_p = []
    rt_pcodaxc = frtp[0].data[0: (min_lwi + max_lwi) // 2 * rt_pcoda[0].
                              stats.sampling_rate]
    pcoda_rotatexc = frotate[1][0: (min_lwi + max_lwi) // 2 * facp[0].
                                stats.sampling_rate]

    corrcoefs_p, thres_p = Get_corrcoefs(rt_pcoda[0], rt_pcodaxc, facp,
                                         pcoda_rotatexc, sec_p, station)

    print 'Backazimuth analysis...'

    ind = max_lwi * facp[0].stats.sampling_rate

    corrbazp, maxcorrp, backas, max_coefs_10deg_p =\
        backas_analysis(frtp[0], rt_pcodaxc, facp, sec_p, corrcoefs_p, ind,
                        station)

    Xp, Yp = np.meshgrid(np.arange(0, sec_p * len(corrcoefs_p), sec_p), backas)

    # TAPER
    rt_pcoda.taper(max_percentage=0.05)

    time_p = rt_pcoda[0].stats.delta * np.arange(0, len(rt_pcoda[0].data))
    fact1_p = 2 * max(rt_pcoda[0].data[0: max_lwi * ac_pcoda[0].stats.
                                       sampling_rate])
    c1_p = .5 * max(abs(pcoda_rotate[1][0: max_lwi * ac_pcoda[0].stats.
                                        sampling_rate])) /\
        max(abs(rt_pcoda[0].data[0: max_lwi * ac_pcoda[0].stats.
                sampling_rate]))

    maxcorrp_over50 = []
    for i10 in range(0, len(maxcorrp)):
        if np.max(corrbazp[:, i10]) >= 0.5:
            maxcorrp_over50.append(maxcorrp[i10])
        else:
            maxcorrp_over50.append(0)

    # print maxcorrp_over50
    print 'Cross-correlation figure for the P-coda...'
    plt.figure(figsize=(18, 9))
    plt.subplot2grid((5, 26), (0, 0), colspan=25)
    plt.plot(time_p, ac_pcoda.select(component='Z')[0].data, color='g')
    plt.ylabel(ur'a$_\mathbf{Z}$ [nm/s$^2$]', fontweight='bold', fontsize=11)
    plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
    plt.xlim(0, (min_lwi + max_lwi) // 2)
    plt.ylim(min(ac_pcoda.select(component='Z')[0].
                 data[0: max_lwi * ac_pcoda[0].stats.sampling_rate]),
             max(ac_pcoda.select(component='Z')[0].data[0: max_lwi *
                                                        ac_pcoda[0].stats.
                                                        sampling_rate]))
    plt.title(ur'$\dot{\mathbf{\Omega}}_\mathbf{z}$ and a$_\mathbf{T}$'
              'correlation in the P-coda in a %d seconds time window'
              ' (highpass, cutoff: 1 Hz). Event: %s %s' % (sec_p, startev.date, startev.time))
    plt.axvline(x=min_pw, linewidth=1)
    plt.axvline(x=min_sw, linewidth=1)
    plt.subplot2grid((5, 26), (1, 0), colspan=25, rowspan=2)
    plt.plot(time_p, rt_pcoda[0].data, color='r', label=ur'Rotation rate')
    plt.plot(time_p, (0.5 / c1_p) * pcoda_rotate[1] + fact1_p, color='k',
             label=ur'Transversal acceleration')
    plt.ylabel(ur'$\dot{\mathbf{\Omega}}_\mathbf{z}$ [nrad/s] -'
               'a$_\mathbf{T}$/2c [1/s]', fontweight='bold', fontsize=11)
    plt.xlim(0, (min_lwi + max_lwi) // 2)
    plt.ylim(min(rt_pcoda[0].data[0: max_lwi * ac_pcoda[0].stats.
             sampling_rate]), fact1_p + max((1. / (2. * c1_p)) *
             pcoda_rotate[1][0: max_lwi * ac_pcoda[0].stats.sampling_rate]))
    xlim2 = (min_lwi + max_lwi) // 2
    box_yposition2 = (fact1_p + max((1. / (2. * c1_p)) * pcoda_rotate[1]
                      [0: max_lwi * ac_pcoda[0].stats.sampling_rate]) -
                      np.abs(min(rt_pcoda[0].data[0: max_lwi * ac_pcoda[0].
                             stats.sampling_rate])))/2.
    plt.axvline(x=min_pw, linewidth=1)
    plt.annotate('P-arrival', xy=(min_pw+xgap*float(xlim2/xlim1),
                                  box_yposition2), fontsize=14,
                 fontweight='bold', bbox=bbox_props)
    plt.axvline(x=min_sw, linewidth=1)
    plt.annotate('S-arrival', xy=(min_sw+xgap*float(xlim2/xlim1),
                                  box_yposition2), fontsize=14,
                 fontweight='bold', bbox=bbox_props)
    plt.grid(True)
    plt.legend(loc=6, shadow=True)

    plt.subplot2grid((5, 26), (3, 0), colspan=25)
    plt.plot(np.arange(0, sec_p * len(corrcoefs_p), sec_p), corrcoefs_p, '.k')
    plt.ylabel(ur'X-corr. coeff.', fontweight='bold')
    plt.xlim(0, (min_lwi + max_lwi) // 2)
    plt.ylim(0, 1)
    plt.grid(True)
    plt.subplot2grid((5, 26), (4, 0), colspan=25)
    plt.pcolor(Xp, Yp, corrbazp, cmap=plt.cm.RdYlGn_r)
    plt.plot(np.arange(0, sec_p * len(corrcoefs_p), sec_p),
             maxcorrp_over50, '.k')
    plt.xlim(0, (min_lwi + max_lwi) // 2)
    plt.xlabel(ur'Time [s]', fontweight='bold')
    plt.ylabel(ur'BAz [°]', fontweight='bold')
    plt.ylim([0, 360])
    plt.yticks([0,60,120,180,240,300,360])
    plt.grid(True)

    fig = plt.subplot2grid((5, 26), (4, 25))
    norm = mpl.colors.Normalize(vmin=-1, vmax=1)
    cb1 = mpl.colorbar.ColorbarBase(
        fig, cmap=plt.cm.RdYlGn_r, norm=norm, orientation='vertical')
    cb1.set_label(ur'X-corr. coeff.', fontweight='bold')
    cb1.set_ticks([-1.0,-0.75,-0.5,-0.25,0.0,0.25,0.5,0.75,1.0])
    print 'Done.'
    print 'Saving figures...'
    plt.savefig(tag_name + '/' + tag_name[23:] + '_page_4.png')
    plt.close()

    print 'Storing information about the event...'

    store_info_json(rotate, ac, rt, corrcoefs, baz, arriv_p, EBA, tag_name,
                    station, phasv_means, phasv_stds, startev, event, net_r,
                    net_s, chan1, chan2, chan3, chan4, sta_r, sta_s, source,
                    loc_r, loc_s, event_source, depth, displacement,
                    event.magnitudes[0]['mag'], 0.001*baz[0], max_ebaz_xcoef)

    print 'Done.'

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Comparison of transversal\
        acceleration and vertical rotation rate through direct waveform\
        comparison in different time windows and cross-correlation analysis.')
    parser.add_argument('--station', help='Choice of station: WET or PFO\
        (default is WET)', type=str, default='WET')
    parser.add_argument('--mode', help='Choice of executive mode for WET: \
        fetch event data from GCMT catalog or get it from a link (for plotting\
        beachballs)(default: neries, otherwise: link)', type=str,
                        default='neries')
    parser.add_argument('--filter', help='Choice of filtering: \
        If True, small Mw events at high distances are filtered out according \
        to logarithmic filter function)(default: False, otherwise: True)',
                        type=str, default='False')
    parser.add_argument('--check_files', help='Erase event folders if complete \
        processing failed, e.g. due to unavailable data \
        (default: False, otherwise: True)',type=str, default='False')
    parser.add_argument('--min_magnitude', help='Minimum magnitude for \
        events (default is 3).', type=float or int, default=4.0)
    parser.add_argument('--max_magnitude', help='Maximum magnitude for \
        events (default is 10).', type=float or int, default=10.0)
    parser.add_argument('--min_depth', help='Minimum depth for events in km \
        (default is 0 km). Negative down.', type=float or int, default=0.0)
    parser.add_argument('--max_depth', help='Maximum depth for events in km \
        (default is -1000).', type=float or int, default=1000.)
    parser.add_argument('--min_latitude', help='Minimum latitude for events.\
        Format +/- 90 decimal degrees (default is -90°).', type=float or int,
                        default=-90.0)
    parser.add_argument('--max_latitude', help='Maximum latitude for events \
        (default is 90°).', type=float, default=90.0)
    parser.add_argument('--min_longitude', help='Minimum longitude for \
        events. Format +/- 180 decimal degrees (default is -180°).',
                        type=float or int, default=-180.0)
    parser.add_argument('--max_longitude', help='Maximum longitude for \
        events (default is 180°).', type=float or int, default=180.0)
    parser.add_argument('--min_datetime', help='Earliest date and time for \
        the search. Format is UTC: yyyy-mm-dd-[T hh:mm:ss]. \
        Example: 2010-02-27T05:00', type=str, default=str(datetime.datetime.now()-datetime.timedelta(hours=168)))
    parser.add_argument('--max_datetime', help='Latest date and time for \
        the search (default is today).',type=str, default=str(datetime.datetime.now()))

    args = parser.parse_args()
    station = args.station
    mode = args.mode
    filter_e = args.filter
    check_files = args.check_files
    print 'Downloading Events...'

    # This link leads to the quakeml webpage of a certain event on IRIS.
    # The fetched xml provides a moment tensor data for the beachballs.
    link = 'http://www.iris.edu/spudservice/momenttensor/736631/quakeml'
    #quakeml = '/home/jsalvermoser/waveformCompare/500kmradius_events_qml.xml'
    quakeml = '/import/two-data/salvermoser/waveformCompare/XML-extra_events/extra_events.xml'
    if station == 'WET' and mode == 'qmlfile':
        cat = read_events(quakeml, format='QUAKEML')
        event_source = "ISC_"
        catalog='ISC_'
    elif mode == 'fdsn':
        catalog='GCMT'
        event_source = "IRIS"
        c = fdsnClient(event_source)
        cat = c.get_events(minmagnitude=args.min_magnitude,
                           maxmagnitude=args.max_magnitude,
                           mindepth=args.min_depth, #magnitudetype='Mw',
                           maxdepth=args.max_depth,
                           minlatitude=args.min_latitude,
                           maxlatitude=args.max_latitude,
                           minlongitude=args.min_longitude,
                           maxlongitude=args.max_longitude,
                           starttime=args.min_datetime,
                           endtime=args.max_datetime,
                           catalog=catalog)
    else:
        catalog='GCMT'
        event_source = "GCMT"
        if UTCDateTime(args.min_datetime) < UTCDateTime(2014,1,1):
            cat_all = read_events('NDK_events_before2014.ndk')
        else:
            cat_all = read_events('http://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog/NEW_QUICK/qcmt.ndk')

        cat = cat_all.filter('time > '+str(args.min_datetime), 'time < '+str(args.max_datetime),
                             'magnitude >= '+str(args.min_magnitude), 'magnitude <= '+str(args.max_magnitude),
                             # 'depth >= '+str(args.min_depth), 'depth <= '+str(args.max_depth)
                             'longitude >= '+str(args.min_longitude), 'longitude <= '+str(args.max_longitude),
                             'latitude >= '+str(args.min_latitude), 'latitude <= '+str(args.max_latitude))

    print 'Downloaded %i events. Starting processing...' % len(cat)

    if not os.path.exists('database/static/OUTPUT'): #'database/static/OUTPUT'):
        os.makedirs('database/static/OUTPUT')#'database/static/OUTPUT')
    if not os.path.exists('database/static/OUTPUT/Phase_velocities'):# 'database/static/OUTPUT/Phase_velocities'):
        os.makedirs('database/static/OUTPUT/Phase_velocities') #'database/static/OUTPUT/Phase_velocities')

    contador1 = 0
    contador2 = 0
    contador3 = 0
    for event in cat:
        print '--------------------------------------------------------------'
        print str(event).split('\n')[0]

        print event.preferred_origin().time
        print '--------------------------------------------------------------'
        try:
            tag_name = os.path.join('database/static/OUTPUT', catalog + '_' + str(event.origins[0]['time'].date) + 
                       'T' + str(event.origins[0]['time'].hour) + ':' + 
                       str(event.origins[0]['time'].minute) +'_' + str(event.magnitudes[0]['mag']) + 
                    '_' + str(event.event_descriptions[0]['text'].splitlines()[0].replace(' ', '_').replace(',', '_')))

            if os.path.exists(str(tag_name)):
                print 'This event was already processed...'
                contador1 += 1
            elif not os.path.exists(str(tag_name)):  # mk directory for each event
                os.makedirs(str(tag_name))
                event.write(tag_name + "/" + catalog + '_' + 
                    str(event.origins[0]['time'].date) + 'T' + str(event.origins[0]['time'].hour) + ':' + 
                    str(event.origins[0]['time'].minute) +'_' + str(event.magnitudes[0]['mag']) + 
                    '_' + str(event.event_descriptions[0]['text'].splitlines()[0].replace(' ', '_').replace(',', '_')) + 
                    ".xml", format="QUAKEML")
                if check_files==False:
                    try:
                        plotWaveformComp(event, station, link, mode, filter_e,
                                         event_source)
                        contador2 += 1
                    except Exception as e:
                        contador3 += 1
                        print e
                else:
                    try:
                        plotWaveformComp(event, station, link, mode, filter_e,
                                         event_source)
                        contador2 += 1
                    except Exception as e:
                        contador3 += 1
                        print e
                        print 'Remove incomplete folder ...'
                        shutil.rmtree(tag_name)
                        
        except IndexError:
            print 'No Magnitude picked for this Event'

    print 'Done, no more events to show.'
    print 'From a total of %i event(s):\n %i was/were successfully processed.' \
          '\n %i could not be processed. \n %i already processed.' % (len(cat), contador2, contador3, contador1)

### DEBUGGER
#from IPython.core.debugger import Tracer; Tracer(colors="Linux")()
