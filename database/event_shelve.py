#!/usr/bin/env python
# -*- coding: utf-8 -*-

from collections import OrderedDict
from flask import Response
import fnmatch
import obspy
import os
import re
import json
import glob
import shelve
import warnings
import uuid
import numpy as np
import geojson
from geojson import Feature, FeatureCollection, Point
from obspy.core import UTCDateTime


class EventShelveException(Exception):
    """
    Exception raised by this module.
    """
    pass


class EventShelveWarning(UserWarning):
    """
    Warning raised by this module.
    """
    pass

def locations2degrees(lat1, long1, lat2, long2):
    lat1=float(lat1); long1=float(long1); lat2=float(lat2); long2=float(long2);
    print lat1, type(lat1)
    lat1 = np.radians(lat1)
    lat2 = np.radians(lat2)
    long1 = np.radians(long1)
    long2 = np.radians(long2)
    print 'hi'
    gds = []
    for i in range(0,len(long2)):
        long_diff = long2[i] - long1
        gd = np.degrees(
            np.arctan2(
                np.sqrt((
                    np.cos(lat2[i]) * np.sin(long_diff)) ** 2 +
                    (np.cos(lat1) * np.sin(lat2[i]) - np.sin(lat1) *
                        np.cos(lat2[i]) * np.cos(long_diff)) ** 2),
                np.sin(lat1) * np.sin(lat2[i]) + np.cos(lat1) * np.cos(lat2[i]) *
                np.cos(long_diff)))
        gds.append(gd)
    print 'hi'
    return gds

def locations2degrees_single(lat1, long1, lat2, long2):
    lat1=float(lat1); long1=float(long1); lat2=float(lat2); long2=float(long2);
    lat1 = np.radians(lat1)
    lat2 = np.radians(lat2)
    long1 = np.radians(long1)
    long2 = np.radians(long2)
    long_diff = long2 - long1
    gd = np.degrees(
        np.arctan2(
            np.sqrt((
                np.cos(lat2) * np.sin(long_diff)) ** 2 +
                (np.cos(lat1) * np.sin(lat2) - np.sin(lat1) *
                    np.cos(lat2) * np.cos(long_diff)) ** 2),
            np.sin(lat1) * np.sin(lat2) + np.cos(lat1) * np.cos(lat2) *
            np.cos(long_diff)))
    return gd

class EventShelve(object):
    def __init__(self, shelve_path, root_folder, quakeml_glob_expr,
                 regex_expr=None):
        """
        Initializes the EventShelve object.
        :param shelve_path:
        :param root_folder:
        :param quakeml_glob_expr:
        :param regex_expr:
        """
        self._s = shelve.open(shelve_path)

        # First step is to get all files.
        quakeml_filenames = []
        for root, _, filenames in os.walk(root_folder):
            for filename in fnmatch.filter(filenames, quakeml_glob_expr):
                quakeml_filenames.append(os.path.normpath(os.path.relpath(
                    os.path.join(root, filename))))

        quakeml_filenames = set(quakeml_filenames)
        filenames_in_shelve = set(self._s.keys())




        # Delete all files no longer available.
        to_be_removed = filenames_in_shelve - quakeml_filenames
        for filename in to_be_removed:
            del self._s[filename]

        filenames_in_shelve = set(self._s.keys())
        print(filenames_in_shelve)

        # Find files that need to be added.
        to_be_added = quakeml_filenames - filenames_in_shelve
        for _i, filename in enumerate(to_be_added):
            print("Indexing file %i of %i: %s ..." % (_i + 1, len(to_be_added), filename))
            cat = obspy.readEvents(filename)
            if len(cat) == 0:
                continue
            elif len(cat) > 1:
                msg = ("File '%s' contains %i events. Only one event per "
                       "file is supported. Will be skipped." %
                       (filename, len(cat)))
                warnings.warn(msg)
                continue

            ev = cat[0]
            with open(filename[:-3]+'json') as data_file:
                data = json.load(data_file)
            peak_rot = data['peak_vertical_rotation_rate']
            max_xcoef = data['peak_correlation_coefficient']
            theo_baz = data['theoretical_backazimuth']
            SNR_rot = data['vertical_rotation_rate_SNR']


            # Get the event id used for that event.
            event_id = None
            if regex_expr is not None:
                match = re.match(regex_expr, ev.resource_id.id)
                if match:
                    try:
                        event_id = match.group(1)
                    except IndexError:
                        pass
            if not event_id:
                event_id = str(uuid.uuid4())

            origin = ev.preferred_origin() or ev.origins[0]
            magnitude = ev.preferred_magnitude() or ev.magnitudes[0]

            event_info = {
                "event_id": event_id,
                "latitude": origin.latitude,
                "longitude": origin.longitude,
                "time": origin.time,
                "depth_in_km": origin.depth / 1000.0,
                "magnitude": magnitude.mag,
                "magnitude_type": magnitude.magnitude_type,
                "peak_vertical_rotation_rate": peak_rot,
                "max_correlation": max_xcoef,
                "theoretical_BAZ": theo_baz,
                "SNR_rotation_rate": SNR_rot
            }

            self._s[filename] = event_info

        # Copy to in memory dictionary.
        self.events = OrderedDict(self._s)

        # Close shelve.
        self._s.close()

    def query(self, starttime=None, endtime=None, minlatitude=None,
              maxlatitude=None, minlongitude=None, maxlongitude=None,
              latitude=None, longitude=None, maxradius=None, minradius=None,
              mindepth=None, maxdepth=None, minmagnitude=None,
              maxmagnitude=None, mincor=None, maxcor=None, minpeakrot=None, minSNR=None, 
              limit=None, offset=1, orderby="time", event_id=None, lat_circ=None, 
              lon_circ=None, circ_dist=None, query_id=None,format='quakeml', **kwargs):
        """
        FDSN event service like queries.
        """
        counter = 0
        actually_used_counter = 0

        found_events = {}
        liste=[]
        # Find all events according to the query.
        for filename, event in self.events.iteritems():
            if lon_circ == None: # need to set default values for lat lon here to calculate dist
                lon_circ = 12.88
            if lat_circ == None:
                lat_circ = 49.15
            dist = locations2degrees_single(lat_circ, lon_circ, event["latitude"], event["longitude"])
            # if format == None: #set default to quakeml to make the fdsn service work!
            #     format='quakeml'
            if (event_id is None or event["event_id"] == event_id) and \
                    (starttime is None or event["time"] >= starttime) and \
                    (endtime is None or event["time"] <= endtime) and \
                    (minlatitude is None or
                     event["latitude"] >= float(minlatitude)) and \
                    (maxlatitude is None or
                     event["latitude"] <= float(maxlatitude)) and \
                    (minlongitude is None or
                     event["longitude"] >= float(minlongitude)) and \
                    (maxlongitude is None or
                     event["longitude"] <= float(maxlongitude)) and \
                    (circ_dist is None or
                     dist <= float(circ_dist)) and \
                    (mindepth is None or
                     event["depth_in_km"] >= float(mindepth)) and \
                    (maxdepth is None or
                     event["depth_in_km"] <= float(maxdepth)) and \
                    (mincor is None or
                     event["max_correlation"] >= float(mincor)) and \
                    (maxcor is None or
                     event["max_correlation"] <= float(maxcor)) and \
                    (minpeakrot is None or
                     event["peak_vertical_rotation_rate"] >= float(minpeakrot)) and \
                    (minSNR is None or
                     event["SNR_rotation_rate"] >= float(minSNR)) and \
                    (minmagnitude is None or
                     event["magnitude"] >= float(minmagnitude)) and \
                    (maxmagnitude is None or
                     event["magnitude"] <= float(maxmagnitude)):
                # counter += 1
                # if counter <= offset:
                #     continue
                actually_used_counter += 1
                if limit is not None and limit >= actually_used_counter:
                    break
                found_events[filename] = event


        print "Found events:", len(found_events)

        # this is just used as a test feature collection: if the catalog does not contain
        # events for the specified parameters, we create a dummy collection that contains
        # more than the allowed number of events to being able to filter in the index.html
        if format=='map':
            if len(found_events) == 0:
                ft = Feature(geometry=Point([0,0]))
                fts = [ft]*2501
                Feature_Coll = FeatureCollection(fts) 
            elif len(found_events)>2500: # this gives a variable that is only used for error messages
                Feature_Coll = FeatureCollection(Feature(geometry=Point([0,0])))
            else:
                if len(found_events.keys()) == 0:
                    msg = ("Could not find events for the specified parameters")
                    warnings.warn(EventShelveWarning)  
                else:
                    features=[]
                    for ev in range(0,len(found_events.keys())):
                        URL_P1= str(found_events.keys()[ev])[:-4] + '_page_1.png'
                        URL_P2= str(found_events.keys()[ev])[:-4] + '_page_2.png'
                        URL_P3= str(found_events.keys()[ev])[:-4] + '_page_3.png'
                        URL_P4= str(found_events.keys()[ev])[:-4] + '_page_4.png'
                        lon = found_events[str(found_events.keys()[ev])]['longitude']
                        lat = found_events[str(found_events.keys()[ev])]['latitude']
                        time = found_events[str(found_events.keys()[ev])]['time']
                        depth = found_events[str(found_events.keys()[ev])]['depth_in_km']
                        magnitude = found_events[str(found_events.keys()[ev])]['magnitude']
                        #region = found_evstr(found_events.keys()[ev])]['text']
                        p = Point([lon,lat])
                        f=Feature(geometry=p, properties={"time" : str(time), "magnitude" : str(magnitude),
                         "depth" : str(depth),
                         "page1_URL": URL_P1,"page2_URL" : URL_P2, "page3_URL" : URL_P3,
                         "page4_URL" : URL_P4})
                        features.append(f)

                    Feature_Coll = FeatureCollection(features)
            return Feature_Coll
        elif format=='quakeml':
            if query_id is None:
                query_id = "smi:local/%s" % str(uuid.uuid4())
            else:
                query_id = "smi:" + query_id.replace("http://", "")

            cat_str = ("<?xml version='1.0' encoding='utf-8'?>\n"
                   '<ns0:quakeml xmlns:ns0="http://quakeml.org/xmlns/quakeml/'
                   '1.2" xmlns:rotational_seismology_database="http://www.rotational-seismology.org" xmlns="http://quakeml.org/xmlns/bed/1.2">\n'
                   '  <eventParameters publicID="%s">\n'
                   "    {events}\n"
                   "  </eventParameters>\n"
                   "</ns0:quakeml>" % query_id)

            pattern = re.compile(r"<event\s.*<\/event>", re.DOTALL)
            event_strings = []
            for filename in found_events.iterkeys():
                with open(filename, "rt") as fh:
                    event_str = fh.read()
                    event_str = re.findall(pattern, event_str)[0]
                    if event_str is None:
                        msg = ("Could not extract event string from event '%'. "
                               "Will be skipped." % filename)
                        warnings.warn(EventShelveWarning)
                        continue
                    event_strings.append(event_str)
            cat_str = cat_str.format(events="\n    ".join(event_strings))
            return Response(cat_str, content_type='text/xml; charset=utf-8')


    def __del__(self):
        self._s.close()
