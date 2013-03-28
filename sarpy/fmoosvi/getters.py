# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 14:41:54 2013

@author: fmoosvi
"""

import sarpy
import json
import collections

def get_num_slices(scan_object, pdata_num = 0):
    
    """
    Returns the number of slices

    :param object scan_object: scan object from a study
    :param integer pdata_num: reconstruction number, according to python numbering.
            default reconstruction is pdata_num = 0.
    :return: integer num_slices: number of slices, detrmined from visu_pars
    """     
    
       ## Ridiculouly long (but definitely complete and correct) method of 
       #  obtaining the third dimension from visu_pars

    if scan_object.pdata[pdata_num].visu_pars.VisuCoreDim == 2:
        num_slices = 1 # we hope this will get updated below in the case of
              # multi-slice 2D data
    else:
        num_slices = scan_object.pdata[pdata_num].visu_pars.VisuCoreSize[2]
    try:
        VisuFGOrderDesc = scan_object.pdata[pdata_num].visu_pars.VisuFGOrderDesc
    except KeyError:
        # data missing, we should have all we need anyway
        pass
    else:
        for VisuFGOrderDesc_element in VisuFGOrderDesc:
            if VisuFGOrderDesc_element[1].strip()=='<FG_SLICE>':
                num_slices = int(VisuFGOrderDesc_element[0])
                
    return num_slices


def get_patients_from_experiment(Experiment_name, verbose = False):
    
    subject_list = get_unique_list_elements(sarpy.Experiment(Experiment_name).get_SUBJECT_id())

    Patients = []  
    
    for subject in subject_list:
        curr_patient = sarpy.Patient(subject)
        Patients.append(curr_patient)        
        
        if verbose:
            print('Patient {0} has {1} sessions.'.format(subject,len(curr_patient.studies)))
    
    return Patients
    
def get_unique_list_elements(list, idfun=None):

    #Get uniue items in list (maintaining order)
    # From: http://www.peterbe.com/plog/uniqifiers-benchmark
    
   # order preserving
   if idfun is None:
       def idfun(x): return x
   seen = {}
   result = []
   for item in list:
       marker = idfun(item)
       if marker in seen: continue
       seen[marker] = 1
       result.append(item)
   return result         

def magic_tuple(scn):
    return (scn.shortdirname, scn.acqp.ACQ_protocol_name)

if __name__ == '__main__':
#    pat_list = get_patients_from_experiment('NecS3', verbose = True)
#    
#    master_sheet = collections.OrderedDict()
#    for pat in pat_list:
#        nm = pat.get_SUBJECT_id()
#        assert len(set(nm)) ==1
#        master_sheet[nm[0]]={}
#        
#        if pat.get_SUBJECT_id()[0] == 'NecS3Hs12':
#            stdy1= pat.studies[0]
#            stdy2= pat.studies[3]
#        elif pat.get_SUBJECT_id()[0] == 'NecS3Hs04':
#            stdy1= pat.studies[0]
#            stdy2= pat.studies[2]
#        elif len(pat.studies) != 2:
#            print 'not 2 studies: %s ' % pat.get_SUBJECT_id()
#            continue
#        else:
#            stdy1 = pat.studies[0]
#            stdy2 = pat.studies[1]
#
#        master_sheet[nm[0]]['0h'] = stdy1.shortdirname
#        master_sheet[nm[0]]['24h'] = stdy2.shortdirname
#        
#        #LL
#        scn = stdy1.find_scan_by_protocol('04_')
#        if len(scn) > 1:
#            print 'more than one LL on day 1 %s' % len(scn) + scn[0].shortdirname
#        master_sheet[nm[0]]['0h-LL']= scn[0].shortdirname
#        scn = stdy2.find_scan_by_protocol('04_')
#        if len(scn) > 1:
#            print 'more than one LL on day 2 %s' % len(scn) + scn[0].shortdirname
#        master_sheet[nm[0]]['24h-LL']= scn[0].shortdirname
#
#        #DCE 1 & 2         
#        scn = stdy1.find_scan_by_protocol('06_')
#        if len(scn) >= 2:
#            master_sheet[nm[0]]['0h-DCE1']= magic_tuple(scn[-2])
#            master_sheet[nm[0]]['0h-DCE2']= magic_tuple(scn[-1])
#        elif len(scn) == 1:
#            print '1 DCE on day 1 %s ' % stdy1.shortdirname
#            
#        scn = stdy2.find_scan_by_protocol('06_')
#        if len(scn) >= 2:
#            master_sheet[nm[0]]['24h-DCE1']= magic_tuple(scn[-2])
#            master_sheet[nm[0]]['24h-DCE2']= magic_tuple(scn[-1])
#        elif len(scn) ==1:
#            print '1 DCE on day 2 %s ' % stdy1.shortdirname
#          
#         # RARE
#        scn = stdy1.find_scan_by_protocol('05_.*RARE')
#        if len(scn) == 2:
#            master_sheet[nm[0]]['0h-IR_A']= magic_tuple(scn[0])
#            master_sheet[nm[0]]['0h-IR_B']= magic_tuple(scn[1])
#        elif len(scn) == 1:
#            master_sheet[nm[0]]['0h-IR_A']= magic_tuple(scn[0])
#        else:
#            print 'RARE not clear, day1 %s ' % stdy1.shortdirname
#            
#        scn = stdy2.find_scan_by_protocol('05_.*RARE')
#        if len(scn) == 2:
#            master_sheet[nm[0]]['24h-IR_A']= magic_tuple(scn[0])
#            master_sheet[nm[0]]['24h-IR_B']= magic_tuple(scn[1])
#        elif len(scn) == 1:
#            master_sheet[nm[0]]['24h-IR_A']= magic_tuple(scn[0])
#        else:
#            print 'RARE not clear, day2 %s ' % stdy2.shortdirname
#          
#          
#         
#        
#    with open('/tmp/master_sheet.json','wb') as outfile:
#        json.dump(master_sheet, outfile, indent=4)

    
    with open('/tmp/master_sheet.json','r') as infile:
        x = json.load(infile)
        
        
    for k,v in x.iteritems():
        try:
            print k, x[k]['24h-DCE2'][0]
        except:
            print k, 'too bad'
    