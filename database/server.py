#!/usr/bin/env python
# -*- coding: utf-8 -*-

import flask
from flask import render_template,request, redirect, url_for, jsonify
from flask_flatpages import FlatPages
from flask.ext.cache import Cache
from flask_frozen import Freezer
from wtforms import Form, BooleanField, TextField, PasswordField, validators


import inspect
import obspy
import os
import sys

import config
from event_shelve import EventShelve
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
freezer = Freezer(app)

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
    start = TextField('Start Time', [validators.Length(min=0, max=17)])
    end = TextField('End Time', [validators.Length(min=0, max=17)])
    minmag = TextField('Min. Magnitude', [validators.Length(min=0, max=5)])
    maxmag = TextField('Max. Magnitude', [validators.Length(min=0, max=5)])
    minlon = TextField('Min. Longitude', [validators.Length(min=0, max=9)])
    maxlon = TextField('Max. Longitude', [validators.Length(min=0, max=9)])
    minlat = TextField('Min. Latitude', [validators.Length(min=0, max=9)])
    maxlat = TextField('Max. Latitude', [validators.Length(min=0, max=9)])

@app.route('/', methods=['GET', 'POST'])
def index():
    form = QueryForm(request.form)
    if request.method == 'POST' and form.validate():
        Start = form.start.data
        End = form.end.data
        Minmagnitude = form.minmag.data
        Maxmagnitude = form.maxmag.data
        Minlongitude = form.minlon.data
        Maxlongitude = form.maxlon.data
        Minlatitude = form.minlat.data
        Maxlatitude = form.maxlat.data
    return render_template('index.html', pages=pages, form=form)


# @app.route('/welcome')
# def welcome():
#     return render_template('welcome.html')

# @app.route('/files')
# def show_files():

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


@app.route(ROOT_URL + "application.wadl")
@cache.cached()
def wadl():
    """
    Return the WADL file.
    """
    with open(os.path.join(PATH, "application.wadl"), "rb") as fh:
        wadl_string = fh.read()
    return wadl_string


@app.route("/query", methods=['GET', 'POST'])
def query():
    """
    The actual query route.
    """

    arguments = {key: value for key, value in flask.request.args.items()}

    #form = QueryForm(request.form)

    # Map short to long arguments.
    mappings = {
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
        cat = event_shelve.query(**arguments)
    except Exception as e:
        return str(e), 500, {}


    if cat is None:
        return ("Request was properly formatted and submitted but no data "
                "matches the selection", 204, {})

    return str(cat)

if __name__ == "__main__":
    if config.PUBLIC is True:
        app.run(host="0.0.0.0", port=config.PORT)
    else:
        app.run(port=config.PORT,debug=True)
   
