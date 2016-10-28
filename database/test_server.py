#!/usr/bin/env python
# -*- coding: utf-8 -*-

import flask
from flask import render_template,request, redirect, url_for, jsonify, Response
from flask_flatpages import FlatPages
from flask.ext.cache import Cache
from wtforms import Form, TextField, validators


import inspect
import obspy
import os
import sys

import config
from event_shelve_test import EventShelve
from obspy.fdsn import Client

DEBUG = False
FLATPAGES_AUTO_RELOAD = DEBUG
FLATPAGES_EXTENSION = '.md'

ROOT_URL = "/"

PATH = os.path.dirname(os.path.abspath(inspect.getfile(
    inspect.currentframe())))


app = flask.Flask("FDSNEventService")
cache = Cache(app, config={"CACHE_TYPE": "simple"})
pages = FlatPages(app)

print("Initializing event shelve...")
# Init the event shelve.
event_shelve = EventShelve(
    shelve_path=config.SHELVE_DB_PATH,
    root_folder=config.QUAKEML_ROOT_DIR,
    quakeml_glob_expr=config.QUAKEML_FILES_GLOB,
    regex_expr=config.REGEX_FOR_EVENT_ID)
print("Done initializing event shelve...")

### CLASSES
class QueryForm(Form):
    station = TextField('Station', [validators.Length(min=0, max=25)])
    start = TextField('Start Time', [validators.Length(min=0, max=17)])
    end = TextField('End Time', [validators.Length(min=0, max=17)])
    minmag = TextField('Min. Magnitude', [validators.Length(min=0, max=5)])
    maxmag = TextField('Max. Magnitude', [validators.Length(min=0, max=5)])
    mindepth = TextField('Min. Depth', [validators.Length(min=0, max=5)])
    maxdepth = TextField('Max. Depth', [validators.Length(min=0, max=5)])
    minlon = TextField('Min. Longitude', [validators.Length(min=0, max=9)])
    maxlon = TextField('Max. Longitude', [validators.Length(min=0, max=9)])
    minlat = TextField('Min. Latitude', [validators.Length(min=0, max=9)])
    maxlat = TextField('Max. Latitude', [validators.Length(min=0, max=9)])
    mincor = TextField('Min. Correlation', [validators.Length(min=0, max=5)])
    maxcor = TextField('Max. Correlation', [validators.Length(min=0, max=5)])
    minpeakrot = TextField('Min. Peak Rotation', [validators.Length(min=0, max=5)])
    minSNR = TextField('Min. Rotation Rate SNR', [validators.Length(min=0, max=5)])
    lat_circ = TextField('Latitude (circ)', [validators.Length(min=0, max=5)])
    lon_circ = TextField('Longitude (circ)', [validators.Length(min=0, max=5)])
    circ_dist = TextField('Circle radius', [validators.Length(min=0, max=5)])
    format = TextField('Format', [validators.Length(min=0, max=5)])

@app.route('/', methods=['GET', 'POST'])
def index():
    form = QueryForm(request.form)
    if request.method == 'POST' and form.validate():
        Station = form.station.data
        Start = form.start.data
        End = form.end.data
        Minmagnitude = form.minmag.data
        Maxmagnitude = form.maxmag.data
        Mindepth = form.mindepth.data
        Maxdepth = form.maxdepth.data
        Minlongitude = form.minlon.data
        Maxlongitude = form.maxlon.data
        Minlatitude = form.minlat.data
        Maxlatitude = form.maxlat.data
        Mincorrelation = form.mincor.data
        Maxcorrelation = form.maxcor.data
        Minpeakrotation = form.minpeakrot.data
        Lat_circ = form.lat_circ.data
        Lon_circ = form.lon_circ.data
        Circ_dist = form.Circ_dist.data
        MinSNR = form.minSNR.data
        Format = form.format.data
    return render_template('index_test.html', pages=pages, form=form)

@app.route('/<path:path>/')
def page(path):
    page = pages.get_or_404(path)
    return render_template('page.html', page=page)


@app.route('/tag/<string:tag>/')
def tag(tag):
    tagged = [p for p in pages if tag in p.meta.get('tags', [])]
    return render_template('tag.html', pages=tagged, tag=tag)


@app.route(ROOT_URL + "version")
def version():
    """
    Return the version string of the webservice.
    """
    return "0.0.1"


@app.route(ROOT_URL + "fdsnws/event/1/application.wadl")
@cache.cached()
def wadl():
    """
    Return the WADL file.
    """
    with open(os.path.join(PATH, "application.wadl"), "rb") as fh:
        wadl_string = fh.read()
    return Response(wadl_string, mimetype='text/xml')


@app.route("/fdsnws/event/1/query", methods=['GET', 'POST'])
def query():
    """
    The actual query route.
    """

    arguments = {key: value for key, value in flask.request.args.items()}

    #form = QueryForm(request.form)

    # Map short to long arguments.
    mappings = {
        "station": "station",
        "start": "starttime",
        "end": "endtime",
        "minlat": "minlatitude",
        "maxlat": "maxlatitude",
        "minlon": "minlongitude",
        "maxlon": "maxlongitude",
        "lat": "latitude",
        "lon": "longitude",
        "minmag": "minmagnitude",
        "maxmag": "maxmagnitude",
        "mindepth": "mindepth",
        "maxdepth": "maxdepth",
        "mincor": "mincor",
        "maxcor": "maxcor",
        "minpeakrot": "minpeakrot",
        "minSNR": "minSNR",
        "circ_dist": "circ_dist",
        "lat_circ": "lat_circ",
        "lon_circ": "lon_circ",
        "format": "format"
    }
    for key, value in mappings.items():
        if key in arguments:
            arguments[value] = arguments[key]

    # Convert times.
    if "starttime" in arguments:
        arguments["starttime"] = obspy.UTCDateTime(arguments["starttime"])
    if "endtime" in arguments:
        arguments["endtime"] = obspy.UTCDateTime(arguments["endtime"])
    arguments["query_id"] = flask.request.base_url

    try:
        if "format" not in arguments:
            arguments["format"] = "quakeml"
        cat = event_shelve.query(**arguments)
    except Exception as e:
        return str(e), 500, {}


    if cat is None:
        return ("Request was properly formatted and submitted but no data "
                "matches the selection", 204, {})
    if arguments["format"]=="map":
        return str(cat)
    elif arguments["format"]=="quakeml":
        return cat

if __name__ == "__main__":
    if config.PUBLIC is True:
        app.run(host="0.0.0.0", port=config.PORT)
    else:
        app.run(port=config.PORT,debug=True)
   
