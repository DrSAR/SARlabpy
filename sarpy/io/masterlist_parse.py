import sarpy.helpers
import os

def masterlist_parse(BrukerObject, getScanLabels = False):

	'''parse the file and:

		- get all the patients (easy: file.keys())
		- get all the studies abbreviations 
		- get all the scan labels'''

	masterlist = BrukerObject._masterlist

	# First strip away the General Section

	masterlist.pop('General',None)

	# Then get the patients

	patients = masterlist.keys()

	studies_dict = {}
	scanlabels_dict = {}
	scanlabels = []


	for pat in patients:

		#
		studies_dict[pat] = masterlist[pat].keys()
		scanlabels_dict[pat] = {}

	# Then get all the scan labels

	for k,v in studies_dict.iteritems():

		for stud in v:
			scanlabels_dict[k][stud] = {}

			try:
				scanlabels.append(masterlist[k][stud]['scanlabels'].keys())
				for lbl in masterlist[k][stud]['scanlabels'].keys():
					scanlabels_dict[k][stud][lbl] = masterlist[k][stud]['scanlabels'][lbl]

			except KeyError: # Skip studies that don't exist
				pass

	if getScanLabels:
		# Flatten and uniquify all the lists
		scanlabels = list(set(sarpy.helpers.flatten_list(scanlabels)))
		return scanlabels_dict,scanlabels
	else:
		return scanlabels_dict

def get_scans_from_masterlist(BrukerObject,scan_label_list):

	''' Helper function to get a dictionary of scan labels and scans 

		Input is any bruker object, and a list of scan_labels
	'''


	# First help the user out a bit by taking the input and turning it into a list 
	if type(scan_label_list) is str:

		tmp = []
		tmp.append(scan_label_list)
		scan_label_list = tmp

		print('String converted to list of size 1 \n Next time, please supply a proper list')

	a,b,scanlabels_dict,d = masterlist_parse(BrukerObject)
	scan_dict = {}

	for scan_label in scan_label_list:

		scan_dict[scan_label] = []

		for k,v in scanlabels_dict.iteritems():
			for stud in v:
				if scan_label in v[stud].keys():
					scan_dict[scan_label].append(os.path.join(k+'.'+stud.strip('study ')+'/',v[stud][scan_label]))	

	return scan_dict

