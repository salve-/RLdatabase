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


INFORMATION: This script can be used to pick phase velocities by hand.
Phase velocity arrays are not stored.
"""

import matplotlib as mpl
#mpl.use('Agg')
import matplotlib
import matplotlib.pyplot as plt
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
from obspy.imaging.mopad_wrapper import Beach
from obspy.core.utcdatetime import UTCDateTime
from mpl_toolkits.basemap import Basemap
import os
from obspy.signal.filter import highpass
import numpy as np
import shutil
from scipy.signal import detrend
#import urllib2
from xml.dom.minidom import parseString
from collections import OrderedDict
from scipy.fftpack import fft, ifft, fftfreq


# if matplotlib.__version__ < '1.0':  # Matplotlib 1.0 or newer is necessary
#     raise ValueError('I need Matplotlib version 1.0 or newer.')


class RotationalProcessingException(Exception):
    pass


def download_data(origin_time, net, sta, loc, chan):
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

    c = arclinkClient(user='test@obspy.org')
    st = c.get_waveforms(network=net, station=sta, location='', channel=chan,
                         starttime=origin_time-190,
                         endtime=origin_time+3*3600+10)

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
    depth = origin.depth * 0.001  # Depth in km
    if station == 'WET':
        net_r = 'BW'
        net_s = 'GR'
        sta_r = 'RLAS'
        sta_s = 'WET'
        loc_r = ''
        loc_s = ''
        if origin < UTCDateTime(2010,1,1):
            chan1 = 'BAZ'
        else: 
            chan1 = 'BJZ'
        chan2 = 'BHE'
        chan3 = 'BHN'
        chan4 = 'BHZ'
        source = 'http://erde.geophysik.uni-muenchen.de'
        # ringlaser signal
        rt = download_data(startev, net_r, sta_r, loc_r, chan1)
        # broadband station signal
        acE = download_data(startev, net_s, sta_s, loc_s, chan2)
        acN = download_data(startev,  net_s, sta_s, loc_s, chan3)
        acZ = download_data(startev,  net_s, sta_s, loc_s, chan4)
        ac = Stream(traces=[acE[0], acN[0], acZ[0]])

        for ca in [ac[0], ac[1], ac[2], rt[0]]:
            ca.stats.coordinates = AttribDict()
            ca.stats.coordinates['longitude'] = 12.8782
            ca.stats.coordinates['latitude'] = 49.144001
            ca.stats['starttime'] = startev - 180
            ca.stats['sampling_rate'] = 20.
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

    cutoff_pc = 1.0  # Cut-off frequency for the highpass filter in the P-coda
    if is_local == 'local':
        for trr in (rt + ac):
            trr.data = trr.data[0: 3600 * rt[0].stats.sampling_rate]
        rt.decimate(factor=2)
        ac.decimate(factor=2)
        sec = 5
        cutoff = 2.0  # Cut-off freq for the lowpass filter for local events
    elif is_local == 'non-local':
        rt.decimate(factor=4)
        ac.decimate(factor=4)
        sec = 120
        cutoff = 1.0  # Cut-off freq for the lowpass filter for non-loc events
    else:
        for trr in (rt + ac):
            trr.data = trr.data[0: 3600 * rt[0].stats.sampling_rate]
        rt.decimate(factor=2)
        ac.decimate(factor=2)
        sec = 3
        cutoff = 4.0  # Cut-off freq for the lowpass filter for local events
    return rt, ac, sec, cutoff, cutoff_pc


def remove_instr_resp(rt, ac, station, startev):

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
    rt.detrend(type='linear')

    # make sure start and endtimes match for both instruments
    startaim = max([tr.stats.starttime for tr in (ac + rt)])
    endtaim = min([tr.stats.endtime for tr in (ac + rt)])

    ac.trim(startaim, endtaim, nearest_sample=True)
    rt.trim(startaim, endtaim, nearest_sample=True)

    if station == 'WET':
        rt[0].data = rt[0].data * 1. / 6.3191 * 1e-3  # Rotation rate in nrad/s

        ac.detrend(type='linear')
        # TAPER
        taper_percentage = 0.05
        taper = np.blackman(np.int(len(ac[0].data) * taper_percentage))
        taper_left, taper_right = np.array_split(taper, 2)
        taper = np.concatenate([taper_left,
                                np.ones(len(ac[0].data)-len(taper)),
                                taper_right])
        ac[0].data = ac[0].data * taper
        ac[1].data = ac[1].data * taper
        ac[2].data = ac[2].data * taper

        # acceleration in nm/s^2
        # note: single zero to go from velocity to acceleration

        paz_sts2 = {'poles': [(-0.0367429 + 0.036754j),
                              (-0.0367429 - 0.036754j)],
                    'sensitivity': 0.944019640, 'zeros': [0j], 'gain': 1.0}
        ac.simulate(paz_remove=paz_sts2, remove_sensitivity=True)  # nm/s^2


    else:
        rt[0].data = rt[0].data * 1. / 2.5284 * 1e-3  # Rotation rate in nrad/s
        # rate in nrad/s
        ac.detrend(type='linear')

        # TAPER
        taper_percentage = 0.05
        taper = np.blackman(np.int(len(ac[0].data) * taper_percentage))
        taper_left, taper_right = np.array_split(taper, 2)
        taper = np.concatenate([taper_left,
                                np.ones(len(ac[0].data)-len(taper)),
                                taper_right])
        ac[0].data = ac[0].data * taper
        ac[1].data = ac[1].data * taper
        ac[2].data = ac[2].data * taper

        ac.remove_response(output='ACC', pre_filt=(0.005, 0.006, 30., 35.))

        # to nm/s^2
        for traza in (ac):
            traza.data = 1e9 * traza.data

    return rt, ac

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

def filter_and_rotate(ac, rt, baz, cutoff, cutoff_pc,
                      station):
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

    cop_ac = ac.copy()
    cop_rt = rt.copy()
    cop_ac.filter('highpass', freq=cutoff_pc, corners=2, zerophase=True)
    cop_rt.filter('highpass', freq=cutoff_pc, corners=2, zerophase=True)
    ac.filter('lowpass', freq=cutoff, corners=2, zerophase=True)
    rt.filter('lowpass', freq=cutoff, corners=2, zerophase=True)

    # rotate translational signal to theoretical event backazimuth
    rotate = rotate_ne_rt(ac.select(component=compN)[
                          0].data, ac.select(component=compE)[0].data, baz[2])


    rt_bands = []
    rotate_bands = []
    bandw = [0.005, 0.005556, 0.0066, 0.00769, 0.00833, 0.00909, 0.01, 0.0111, 0.0125, 0.0143, 0.016667, 0.02, 0.025, 0.0333, 0.05]
    fc = [0.005, 0.005556, 0.0066, 0.00769, 0.00833, 0.00909, 0.01, 0.0111, 0.0125, 0.0143, 0.016667, 0.02, 0.025, 0.0333, 0.05]
    for ii in range(0,len(bandw)):
        rt_bands.append(gaussianfilter(rt[0][0:32768],0.05,bandw[ii],fc[ii])) ## to flip *(-1)
        rotate_bands.append(gaussianfilter(rotate[1][0:32768],0.05,bandw[ii],fc[ii]))

    ###

    print 'Done'
    return rotate, cop_rt, rt_bands, rotate_bands,fc


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
    TauPy_model = TauPyModel('ak135')
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
    #from IPython.core.debugger import Tracer; Tracer(colors="Linux")()
    compE, compN = station_components(station)
    corrcoefs = []
    sampling_rate = acstr[1].stats.sampling_rate
    for i5 in xrange(0, len(rodat) // (int(sampling_rate) * sec)):
        coeffs = xcorr(rodat[sampling_rate * sec *
                             i5:sampling_rate * sec * (i5 + 1)],
                       rotate_array[sampling_rate * sec * i5:sampling_rate * sec * (i5 + 1)], 0)
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

    return corrbaz, maxcorr, backas


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
    if rt[0].stats.sampling_rate == 5.0:
        sec2 = 30
    else:
        sec2 = 
    step2 = 1
    backas2 = np.linspace(0, 360 - step2, 360 / step2)
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
            corrbaz2.append(corrbazz2[1])
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
    return corrsum, backas2


def phase_vel(sampl_rate, rt, sec, corrcoefs, rotate, corrsum, backas2, ind_band,
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
                phas_v = .001 * 0.5 * max(rotate[1][sampl_rate *
                                          sec * i8:sampl_rate *
                                          sec * (i8 + 1)]) /\
                    max(rt[0].data[sampl_rate * sec *
                        i8:sampl_rate * sec * (i8 + 1)])
            else:
                phas_v = np.NaN
            phasv.append(phas_v)
        phasv = np.asarray(phasv)  # Velocity in km/s
        phasv=[phasv]

    if ind_band:  # dealing with frequency bands
        import matplotlib.pyplot as pl
        pl.figure()
        pl.plot(rt/np.max(rt),'r')
        pl.plot(rotate/np.max(rotate),'k')
        pl.scatter(np.argmax(rt),1,c='r',s=30)
        pl.scatter(np.argmax(rotate),1,c='k',s=30)
        pl.show()
        print np.argmax(rt)
        ok = raw_input('Is the rotation pick ok? ')
        if ok=='y':
            arg_rt = np.argmax(rt)
        else:
            arg_rt = raw_input('insert argument of maximum rotation rate: ')
        ok2 = raw_input('Is the acceleration pick ok? ')
        if ok2 == 'y':
            arg_rotate = np.argmax(rotate)
        else:
            arg_rotate = raw_input('insert argument of maximum transverse acceleration: ')

        phasv1 = .001 * 0.5 * rotate[int(arg_rotate)] / rt[int(arg_rt)]
        ph=[]
        for i in range(arg_rotate-100,arg_rotate+100):
            ph.append(.001 * 0.5 * rotate[i] / rt[i])
        phasv2 = np.mean(ph)
        print 'picked: ',phasv1
        print 'averaged: ',phasv2
        phasv=[phasv1,phasv2]

        # for i8 in xrange(ind_surf, len(corrcoefs)):
        #     if corrcoefs[i8] >= 0.75:
        #         # Velocity in km/s
        #         phas_v = .001 * 0.5 * max(rotate[sampl_rate *
        #                                   sec * i8:sampl_rate *
        #                                   sec * (i8 + 1)]) /\
        #             max(rt[sampl_rate * sec *
        #                 i8:sampl_rate * sec * (i8 + 1)])
        #     else:
        #         phas_v = np.NaN
        #     phasv.append(phas_v)
        # phasv = np.asarray(phasv)  # Velocity in km/s

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
                    station, phasv_means, phasv_aves_est, phasv_stds, startev, event, net_r,
                    net_s, chan1, chan2, chan3, chan4, sta_r, sta_s, source,
                    loc_r, loc_s, event_source, depth, fc, magnitude, distance):
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
    :param tag_name: Name of the .xml file.
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
                    ('fdsn_data_select_source', source)])),
                ('translational', OrderedDict([
                    ('network', net_s),
                    ('station', sta_s),
                    ('loc', loc_s),
                    ('channel_N', chan3),
                    ('channel_E', chan2),
                    ('channel_Z', chan4),
                    ('fdsn_data_select_source', source)]))
                ])),
            ('event_id', event.resource_id.id),
            ('fdsn_event_source', event_source),
            ('starttime', str(startev-180)),
            ('endtime', str(startev+3*3600)),
            ('magnitude', magnitude),
            ('depth', depth),
            ('depth_unit', 'km'),
            ('epicentral_distance', distance),
            ('epicentral_distance_unit', 'km'),
            ('peak_transverse_acceleration', PAT),
            ('peak_transverse_acceleration_unit', 'nm/s^2'),
            ('theoretical_backazimuth', TBA),
            ('theoretical_backazimuth_unit', 'degree'),
            ('estimated_backazimuth', EBA),
            ('estimated_backazimuth_unit', 'degree'),
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
                    ('freq_center', fc[0]),
                    ('freq_unit', 'Hz'),
                    ('mean_phase_vel', phasv_means[0]),
                    ('averaged_phase_vel', phasv_aves_est[0]),
                    ('vel_std', phasv_stds[0]),
                    ('vel_unit', 'km/s')])),
                ('band_2', OrderedDict([
                    ('freq_center', fc[1]),
                    ('freq_unit', 'Hz'),
                    ('mean_phase_vel', phasv_means[1]),
                    ('averaged_phase_vel', phasv_aves_est[1]),
                    ('vel_std', phasv_stds[1]),
                    ('vel_unit', 'km/s')])),
                ('band_3', OrderedDict([
                    ('freq_center', fc[2]),
                    ('freq_unit', 'Hz'),
                    ('mean_phase_vel', phasv_means[2]),
                    ('averaged_phase_vel', phasv_aves_est[2]),
                    ('vel_std', phasv_stds[2]),
                    ('vel_unit', 'km/s')])),
                ('band_4', OrderedDict([
                    ('freq_center', fc[3]),
                    ('freq_unit', 'Hz'),
                    ('mean_phase_vel', phasv_means[3]),
                    ('averaged_phase_vel', phasv_aves_est[3]),
                    ('vel_std', phasv_stds[3]),
                    ('vel_unit', 'km/s')])),
                ('band_5', OrderedDict([
                    ('freq_center', fc[4]),
                    ('freq_unit', 'Hz'),
                    ('mean_phase_vel', phasv_means[4]),
                    ('averaged_phase_vel', phasv_aves_est[4]),
                    ('vel_std', phasv_stds[4]),
                    ('vel_unit', 'km/s')])),
                ('band_6', OrderedDict([
                    ('freq_center', fc[5]),
                    ('freq_unit', 'Hz'),
                    ('mean_phase_vel', phasv_means[5]),
                    ('averaged_phase_vel', phasv_aves_est[5]),
                    ('vel_std', phasv_stds[5]),
                    ('vel_unit', 'km/s')])),
                ('band_7', OrderedDict([
                    ('freq_center', fc[6]),
                    ('freq_unit', 'Hz'),
                    ('mean_phase_vel', phasv_means[6]),
                    ('averaged_phase_vel', phasv_aves_est[6]),
                    ('vel_std', phasv_stds[6]),
                    ('vel_unit', 'km/s')])),
                ('band_8', OrderedDict([
                    ('freq_center', fc[7]),
                    ('freq_unit', 'Hz'),
                    ('mean_phase_vel', phasv_means[7]),
                    ('averaged_phase_vel', phasv_aves_est[7]),
                    ('vel_std', phasv_stds[7]),
                    ('vel_unit', 'km/s')])),
                ('band_9', OrderedDict([
                    ('freq_center', fc[8]),
                    ('freq_unit', 'Hz'),
                    ('mean_phase_vel', phasv_means[8]),
                    ('averaged_phase_vel', phasv_aves_est[8]),
                    ('vel_std', phasv_stds[8]),
                    ('vel_unit', 'km/s')])),
                ('band_10', OrderedDict([
                    ('freq_center', fc[9]),
                    ('freq_unit', 'Hz'),
                    ('mean_phase_vel', phasv_means[9]),
                    ('averaged_phase_vel', phasv_aves_est[9]),
                    ('vel_std', phasv_stds[9]),
                    ('vel_unit', 'km/s')])),
                ('band_11', OrderedDict([
                    ('freq_center', fc[10]),
                    ('freq_unit', 'Hz'),
                    ('mean_phase_vel', phasv_means[10]),
                    ('averaged_phase_vel', phasv_aves_est[10]),
                    ('vel_std', phasv_stds[10]),
                    ('vel_unit', 'km/s')])),
                ('band_12', OrderedDict([
                    ('freq_center', fc[11]),
                    ('freq_unit', 'Hz'),
                    ('mean_phase_vel', phasv_means[11]),
                    ('averaged_phase_vel', phasv_aves_est[11]),
                    ('vel_std', phasv_stds[11]),
                    ('vel_unit', 'km/s')])),
                ('band_13', OrderedDict([
                    ('freq_center', fc[12]),
                    ('freq_unit', 'Hz'),
                    ('mean_phase_vel', phasv_means[12]),
                    ('averaged_phase_vel', phasv_aves_est[12]),
                    ('vel_std', phasv_stds[12]),
                    ('vel_unit', 'km/s')])),
                ('band_14', OrderedDict([
                    ('freq_center', fc[13]),
                    ('freq_unit', 'Hz'),
                    ('mean_phase_vel', phasv_means[13]),
                    ('averaged_phase_vel', phasv_aves_est[13]),
                    ('vel_std', phasv_stds[13]),
                    ('vel_unit', 'km/s')])),
                ('band_15', OrderedDict([
                    ('freq_center', fc[14]),
                    ('freq_unit', 'Hz'),
                    ('mean_phase_vel', phasv_means[14]),
                    ('averaged_phase_vel', phasv_aves_est[14]),
                    ('vel_std', phasv_stds[14]),
                    ('vel_unit', 'km/s')]))
                ]))
            ])
    outfile = open(tag_name + '/' + tag_name[7:] + '.json', 'wt')
    json.dump(dic, outfile, indent=4)

    outfile.close()


def plotWaveformComp(event, station, mode, filter_e, event_source):

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
    """

    tag_name = os.path.join('OUTPUT', catalog + '_' + str(event.origins[0]['time'].date) +
     'T' + str(event.origins[0]['time'].hour) + ':' + str(event.origins[0]['time'].minute) +
     '_' + str(event.magnitudes[0]['mag']) + '_' + 
     str(event.event_descriptions[0]['text'].splitlines()[0].replace(' ', '_').replace(',', '_')))
    tag_name2 = os.path.join('OUTPUT/Phase_velocities', str(event)
                             .splitlines()[0].replace('\t', '')
                             .replace(' ', '_').replace('|', '_')
                             .replace(':', '_'))

    # event information:
    latter, lonter, depth, startev, rt, ac, baz, gcdist, net_r, net_s, chan1,\
        chan2, chan3, chan4, sta_r, sta_s, loc_r, loc_s, source =\
        event_info_data(event, station, mode)

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

    

    # Preprocesing of rotational and translational signal

    # check if event is local
    # + resample signal accordingly
    # + Different cut-off frequency for the lowpass filter,
    #   depending on the epicentral distance
    rt, ac, sec, cutoff, cutoff_pc =\
        resample(is_local(baz), baz, rt, ac)

    print 'Removing instrument response...'
    rt, ac = remove_instr_resp(rt, ac, station, startev)

    print 'Filtering and rotating...'
    rotate, cop_rt, rt_bands, rotate_bands, fc =\
        filter_and_rotate(ac, rt, baz, cutoff, cutoff_pc, station)

    print 'Getting arrival times...'
    init_sec = startev - ac[0].stats.starttime
    # When the event starts in the fetched data

    arriv_p, arriv_s = ps_arrival_times(baz[0], depth, init_sec)

    min_pw, max_pw, min_sw, max_sw, min_lwi, max_lwi, min_lwf, max_lwf =\
        time_windows(baz, arriv_p, arriv_s, init_sec, is_local(baz))

    rt.taper(max_percentage=0.05)


    # Cross-correlation analysis
    print 'Obtaining zero-lag correlation coefficients for theoretical' \
        ' backazimuth...'
    corrcoefs, thres = Get_corrcoefs(rt[0], rt[0].data, ac, rotate[1], sec,
                                     station)
    corrcoefs_bands = []
    secs = [200, 100, 50, 20, 12, 10, 8, 6, 6 ,6 ,6 ,6 ,6 ,6 , 6]
    # from IPython.core.debugger import Tracer; Tracer(colors="Linux")()

    for i_ in range(0,len(rt_bands)):
        corrcoefs_tmp, thres_tmp = Get_corrcoefs(rt_bands[i_], rt_bands[i_], ac,
                                            rotate_bands[i_], secs[i_], station)
        corrcoefs_bands.append(corrcoefs_tmp)

    # zero-lag correlation coefficients for range of backazimuths
    print 'Backazimuth analysis...'
    corrbaz, maxcorr, backas = \
        backas_analysis(rt[0], rt[0].data, ac, sec, corrcoefs, None, station)

    # Estimating backazimuth
    print 'Estimating backazimuth...'
    corrsum, backas2 = backas_est(rt, ac, min_sw, max_lwf, station)

    # calculate phase veloc. for windows where corrcoef is good enough (.75)
    # optional TODO: threshold value of corrcoef as a parameter
    # TODO phase velocity estimation by misfit minimization
    print 'Phase velocities...'

    # calculate startindex for phase velocity calculation in frequency bands:
    # starts at the beginning of surface wave arrivals as body waves are
    # not appropriate
    ind_surf = int(min_lwi/sec)

    ind_band = False  # indicator that we are dealing with bands 1-8
    phasv, EBA = phase_vel(int(rt[0].stats.sampling_rate), rt, sec, corrcoefs, rotate, corrsum, backas2,
                           ind_band, ind_surf)
    # calculates phase velocities for different frequency bands -
    # -> put EBA_bandx here instead of EBA and backas2_bandx,...
    # for frequency dependent EBA
    ind_band = True  # indicator that we are dealing with bands 1-8
    phasv_means = []
    phasv_stds = []
    phasv_aves_est = [] #averaged over several values not only maxima

    for ii_ in range(0,len(rt_bands)):
        phasv_tmp, EBA_tmp = phase_vel(int(rt[0].stats.sampling_rate), rt_bands[ii_], secs[ii_], corrcoefs_bands[ii_], rotate_bands[ii_],
                                  corrsum, backas2, ind_band, ind_surf)
        # phasv_tmp = phasv_tmp[~np.isnan(phasv_tmp)] # filters out all NaN
        phasv_means.append(phasv_tmp[0])
        phasv_aves_est.append(phasv_tmp[1])
        phasv_stds.append('not defined')
    print phasv_means

    #from IPython.core.debugger import Tracer; Tracer(colors="Linux")()
    ind_band = False

    print 'Storing information about the event...'


    store_info_json(rotate, ac, rt, corrcoefs, baz, arriv_p, EBA, tag_name,
                    station, phasv_means, phasv_aves_est, phasv_stds, startev, event, net_r,
                    net_s, chan1, chan2, chan3, chan4, sta_r, sta_s, source,
                    loc_r, loc_s, event_source, depth, fc,
                    event.magnitudes[0]['mag'], 0.001*baz[0])

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
    print 'Downloading Events...'

    catalog='GCMT'
    event_source = "IRIS"


    if UTCDateTime(args.min_datetime) < UTCDateTime(2014,1,1):
            cat_all = read_events('/home/jsalvermoser/waveformCompare/NDK_events_before2014.ndk')
    else:
        cat_all = read_events('http://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog/NEW_QUICK/qcmt.ndk')

    print cat_all
    cat = cat_all.filter('time > '+str(args.min_datetime), 'time < '+str(args.max_datetime),
                             'magnitude >= '+str(args.min_magnitude), 'magnitude <= '+str(args.max_magnitude),
                             # 'depth >= '+str(args.min_depth), 'depth <= '+str(args.max_depth)
                             'longitude >= '+str(args.min_longitude), 'longitude <= '+str(args.max_longitude),
                             'latitude >= '+str(args.min_latitude), 'latitude <= '+str(args.max_latitude))

    print 'Downloaded %i events. Starting processing...' % len(cat)

    if not os.path.exists('OUTPUT'):
        os.makedirs('OUTPUT')
    if not os.path.exists('OUTPUT/Phase_velocities'):
        os.makedirs('OUTPUT/Phase_velocities')    

    contador1 = 0
    contador2 = 0
    for event in cat:
        print '--------------------------------------------------------------'
        print str(event).split('\n')[0]

        print event.preferred_origin().time
        print '--------------------------------------------------------------'
        tag_name = os.path.join('OUTPUT', catalog + '_' + str(event.origins[0]['time'].date) + 
                   'T' + str(event.origins[0]['time'].hour) + ':' + 
                   str(event.origins[0]['time'].minute) +'_' + str(event.magnitudes[0]['mag']) + 
                '_' + str(event.event_descriptions[0]['text'].splitlines()[0].replace(' ', '_').replace(',', '_')))

        if os.path.exists(str(tag_name)):
            print 'This event was already processed...'
            contador1 = contador1 + 1

        elif not os.path.exists(str(tag_name)):  # mk directory for each event
            os.makedirs(str(tag_name))
            event.write(tag_name + "/" + catalog + '_' + 
                str(event.origins[0]['time'].date) + 'T' + str(event.origins[0]['time'].hour) 
                + ':' + str(event.origins[0]['time'].minute) +'_' + str(event.magnitudes[0]['mag']) + 
                '_' + str(event.event_descriptions[0]['text'].splitlines()[0].replace(' ', '_').replace(',', '_')) 
                + ".xml", format="QUAKEML")
            try:
                plotWaveformComp(event, station, mode, filter_e,
                                 event_source)
                contador2 = contador2 + 1
            except Exception as e:
                contador2 = contador2-1
                print e

    print 'Done, no more events to show.'
    print 'From a total of %i event(s), %i was/were successfully processed. ' \
          '%i already processed.' % (len(cat), contador2, contador1)
