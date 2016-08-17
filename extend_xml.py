from obspy import Catalog, UTCDateTime
from obspy.core.event import read_events
from obspy.core.util.attribdict import AttribDict
import json
import glob

def parse_json(filename,parameter1,parameter2,parameter3,parameter4,parameter5,parameter6,parameter7,parameter8,parameter9):
	with open(filename) as data_file:
		data = json.load(data_file)
	return data[parameter1], data[parameter2], data[parameter3], data[parameter4],\
	data[parameter5], data[parameter6],data[parameter7],data[parameter8],data[parameter9]

path = "/import/two-data/salvermoser/waveformCompare"

filenames_json = glob.glob(path + "/database/static/OUTPUT/*/*.json")
filenames_xml = glob.glob(path + "/database/static/OUTPUT/*/*.xml")

counter=0

for i_ in range(0,len(filenames_json)):
	d,pcc,tSNR,rSNR,tba,eba,peak_tra,freq, p_rot = parse_json(filenames_json[i_],'epicentral_distance',
		'peak_correlation_coefficient','transverse_acceleration_SNR','vertical_rotation_rate_SNR',
		'theoretical_backazimuth','estimated_backazimuth','peak_transverse_acceleration', 
		'frequency_at_peak_vertical_rotation_rate', 'peak_vertical_rotation_rate')

	ns = 'http://www.rotational-seismology.org'

	params = AttribDict()
	params.namespace = ns
	params.value = AttribDict()

	params.value.epicentral_distance = AttribDict()
	params.value.epicentral_distance.namespace = ns
	params.value.epicentral_distance.value = d
	params.value.epicentral_distance.attrib = {'unit':"km"}

	params.value.transverse_acceleration_SNR= AttribDict()
	params.value.transverse_acceleration_SNR.namespace = ns
	params.value.transverse_acceleration_SNR.value = tSNR

	params.value.vertical_rotation_rate_SNR = AttribDict()
	params.value.vertical_rotation_rate_SNR.namespace = ns
	params.value.vertical_rotation_rate_SNR.value = rSNR

	params.value.theoretical_backazimuth = AttribDict()
	params.value.theoretical_backazimuth.namespace = ns
	params.value.theoretical_backazimuth.value = tba
	params.value.theoretical_backazimuth.attrib = {'unit':"degree"}

	params.value.estimated_backazimuth = AttribDict()
	params.value.estimated_backazimuth.namespace = ns
	params.value.estimated_backazimuth.value = eba
	params.value.estimated_backazimuth.attrib = {'unit':"degree"}

	params.value.peak_transverse_acceleration = AttribDict()
	params.value.peak_transverse_acceleration.namespace = ns
	params.value.peak_transverse_acceleration.value = peak_tra
	params.value.peak_transverse_acceleration.attrib = {'unit':"nm/s^2"}

	params.value.frequency_at_peak_rotation = AttribDict()
	params.value.frequency_at_peak_rotation.namespace = ns
	params.value.frequency_at_peak_rotation.value = freq
	params.value.frequency_at_peak_rotation.attrib = {'unit':"Hz"}

	params.value.peak_vertical_rotation_rate = AttribDict()
	params.value.peak_vertical_rotation_rate.namespace = ns
	params.value.peak_vertical_rotation_rate.value = p_rot
	params.value.peak_vertical_rotation_rate.attrib = {'unit':"nrad/s"}

	params.value.peak_correlation_coefficient = AttribDict()
	params.value.peak_correlation_coefficient.namespace = ns
	params.value.peak_correlation_coefficient.value = pcc

	cat = read_events(pathname_or_url=filenames_xml[i_], format='QUAKEML')

	cat[0].extra = AttribDict()
	cat[0].extra.params = params
	cat.write(filenames_xml[i_], "QUAKEML",
		nsmap={html="rotational_seismology_database": ns})
	counter+=1
	print counter