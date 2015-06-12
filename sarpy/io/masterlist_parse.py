import sarpy.helpers

def masterlist_parse(BrukerObject,
					 ):

	'''parse the file and:

		- get all the patients (easy: file.keys())
		- get all the studies abbreviations 
		- get all the scan labels'''

	masterlist = BrukerObject._masterlist

	# First strip away the General Section

	masterlist.pop('General',None)

	# Then get the patients

	patients = masterlist.keys()

	# Then get the studies abbreviations

	studies = []
	studies_dict = {}

	scanlabels = []
	scanlabels_dict = {}

	for pat in patients:

		studies.append(masterlist[pat].keys())

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

	# Flatten and uniquify all the lists

	scanlabels = list(set(sarpy.helpers.flatten_list(scanlabels)))
	studies = list(set(sarpy.helpers.flatten_list(studies)))


	return studies, studies_dict, scanlabels_dict,scanlabels

