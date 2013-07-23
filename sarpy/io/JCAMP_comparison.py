# -*- coding: utf-8 -*-
"""
Copyright: SARlab members, UBC, Vancouver, 2013
"""
from __future__ import print_function

import re
import numpy
import sarpy

# these are the different methods to compare JCAMP parameters by
def regex_comp(a,b): return re.match(a,b) # useful for strings
def arr_comp(a,b): return (a==b).all()    # numpy provides an array comparison
def int_comp(a,b): return a==b # this should work for int, float and lists

default_param = frozenset(['default'])
dictionary = {frozenset([
'ACQ_BF_enable', 'ACQ_calib_date', 'ACQ_coil_config_file', 'ACQ_completed', 'ACQ_contrast_agent', 'ACQ_DS_enabled', 'ACQ_echo_descr', 'ACQ_experiment_mode', 'ACQ_fit_function_name', 'ACQ_flipback', 'ACQ_institution', 'ACQ_method', 'ACQ_movie_descr', 'ACQ_operation_mode', 'ACQ_operator', 'ACQ_patient_pos', 'ACQ_pipe_status', 'ACQ_protocol_location', 'ACQ_protocol_name', 'ACQ_Routing', 'ACQ_Routing_base', 'ACQ_scan_name', 'ACQ_scan_size', 'ACQ_slice_orient', 'ACQ_slice_sepn_mode', 'ACQ_station', 'ACQ_status', 'ACQ_sw_version', 'ACQ_switch_pll_enabled', 'ACQ_time', 'ACQ_transmitter_coil', 'ACQ_trigger_enable', 'ACQ_trigger_reference', 'ACQ_user_filter', 'ACQ_word_size', 'AcqMode', 'AQ_mod', 'BSPulseEnum', 'BYTORDA', 'CalcSpoiler', 'CalculateGradDelay', 'CalculateRPDelay', 'constNEchoes', 'DATATYPE', 'DATE', 'DIGMOD', 'DIGTYP', 'DQDMODE', 'DSPFIRM', 'EchoTimeMode', 'EdcOnOff', 'ExcPulseAuto', 'ExcPulseEnum', 'FILE_LOCATION', 'filename', 'FitFunctionName', 'FQ1LIST', 'FQ2LIST', 'FQ3LIST', 'FQ8LIST', 'FWInphase', 'GO_block_size', 'GO_data_save', 'GO_disp_update', 'GO_init_files', 'GO_macro', 'GO_online_reco', 'GO_raw_data_format', 'GO_reco_display', 'GO_reco_each_nr', 'GO_time_est', 'GO_use_macro', 'GradientDelayYesNo', 'GradOff', 'GradSync', 'GRDPROG', 'GS_continue', 'GS_disp_update', 'GS_get_info_points', 'GS_image_type', 'GS_info_dig_filling', 'GS_info_max_point', 'GS_info_normalized_area', 'GS_online_reco', 'GS_reco_display', 'GS_typ', 'HalfAcquisition', 'HPPRGN', 'INSTRUM', 'InvPulseEnum', 'LOCNUC', 'MeasMode', 'Method', 'NUC1', 'NUC2', 'NUC3', 'NUC4', 'NUC5', 'NUC6', 'NUC7', 'NUC8', 'NUCLEUS', 'OnlineReco', 'OPT_EDCOnOff', 'OPT_ManAdjustment', 'OPT_NavKeepData', 'OPT_RFLOnOff', 'ORIGIN', 'OWNER', 'PAPS', 'POWMOD', 'ProcessingMacro', 'PULPROG', 'PVM_BlBloodOnOff', 'PVM_ChPulse1Enum', 'PVM_ChPulse2Enum', 'PVM_ChPulse3Enum', 'PVM_DecOnOff', 'PVM_DeriveGains', 'PVM_DiffPrepMode', 'PVM_DigAutSet', 'PVM_DigFilter', 'PVM_DigQuad', 'PVM_DwDgSwitch', 'PVM_DwDirectScale', 'PVM_DwMeasMode', 'PVM_DwRfcPulseEnum', 'PVM_DwVisiblePars', 'PVM_EncOrder1', 'PVM_EncOrder2', 'PVM_EncUseMultiRec', 'PVM_EpiAutoGhost', 'PVM_EpiBlipsOff', 'PVM_EpiCombine', 'PVM_EpiDoubleShotAdj', 'PVM_EpiDriftCorr', 'PVM_EpiEchoTimeShifting', 'PVM_EpiGradSync', 'PVM_EpiGrappaSegAdj', 'PVM_EpiNavigatorMode', 'PVM_EpiPrefixNavYes', 'PVM_EpiRampComp', 'PVM_EpiRampForm', 'PVM_EpiRampMode', 'PVM_EpiTrajAdjAutomatic', 'PVM_EpiTrajAdjComp', 'PVM_EpiTrajAdjMeasured', 'PVM_EpiTrajAdjRampform', 'PVM_EpiTrajAdjRampmode', 'PVM_EpiTrajAdjYesNo', 'PVM_EvolutionOnOff', 'PVM_FatSupDeriveGainMode', 'PVM_FatSupOnOff', 'PVM_FatSupprPulseEnum', 'PVM_FlipBackOnOff', 'PVM_FlowSatDeriveGainMode', 'PVM_FovSatOnOff', 'PVM_GeoMode', 'PVM_InFlowSatOnOff', 'PVM_Isotropic', 'PVM_MagTransOnOff', 'PVM_MajSliceOri', 'PVM_MotionSupOnOff', 'PVM_MovieOnOff', 'PVM_NoeOnOff', 'PVM_Nucleus1', 'PVM_Nucleus1Enum', 'PVM_Nucleus2', 'PVM_Nucleus2Enum', 'PVM_Nucleus3Enum', 'PVM_Nucleus4Enum', 'PVM_Nucleus5Enum', 'PVM_Nucleus6Enum', 'PVM_Nucleus7Enum', 'PVM_Nucleus8Enum', 'PVM_ObjOrderScheme', 'PVM_OperationMode', 'PVM_OvsConflict', 'PVM_OvsDeriveGainMode', 'PVM_OvsOnOff', 'PVM_OvsPulseEnum', 'PVM_ppgFlag1', 'PVM_ppgFlag2', 'PVM_ppgFlag3', 'PVM_ppgFlag4', 'PVM_PreemphasisFileEnum', 'PVM_PreemphasisSpecial', 'PVM_RefAttMod1', 'PVM_RefAttStat1', 'PVM_RefScanPCYN', 'PVM_RefScanYN', 'PVM_SatSlicesDeriveGainMode', 'PVM_SatSlicesPulseEnum', 'PVM_ScanTimeStr', 'PVM_SelIrDeriveGainMode', 'PVM_SelIrOnOff', 'PVM_SelIrPulseEnum', 'PVM_ShimCoeffCalVersion', 'PVM_ShimCoeffHwMode', 'PVM_ShimCoeffHwStatus', 'PVM_SpatDimEnum', 'PVM_SpecDimEnum', 'PVM_TaggingDeriveGainMode', 'PVM_TaggingDir', 'PVM_TaggingMode', 'PVM_TaggingOnOff', 'PVM_TriggerMode', 'PVM_TriggerModule', 'PVM_TriggerOutOnOff', 'PVM_TuneShimAdjFreq', 'PVM_TuneShimForceSubset', 'PVM_TuneShimIncIter', 'PVM_TuneShimSubset', 'PVM_UserType', 'PVM_VoxExcOrder', 'PVM_VoxMethodType', 'PVM_VpPulse1Enum', 'PVM_VpPulse2Enum', 'PVM_WsCalcSpoiler', 'PVM_WsCalcTiming', 'PVM_WsDeriveGainMode', 'PVM_WsMode', 'PVM_WsOnOff', 'RampMode', 'RECO_map_mode', 'RECO_wordtype', 'RecoMethMode', 'RecoMode', 'RecoOptimise', 'RecoOptimizeMemory', 'RecoRegridNSetDefaults', 'RecoUseOrigin', 'RefPulseEnum', 'Reorder', 'ResearchMode', 'RFSpoilerOnOff', 'RFSpoilOnOff', 'seg_mode', 'ShimHwInfoAvail', 'ShowAllPars', 'SignalType', 'SlabSel', 'SlabSelect', 'SliceAlter', 'SOLVENT', 'SPNAM0', 'SPNAM1', 'SPNAM2', 'SPNAM3', 'SPNAM4', 'SPNAM5', 'SPNAM6', 'SPNAM7', 'SPNAM8', 'SPNAM9', 'SPNAM10', 'SPNAM11', 'SPNAM12', 'SPNAM13', 'SPNAM14', 'SPNAM15', 'TITLE', 'TraAdjMode', 'TraAutomaticParSetting', 'TraExcPulseEnum', 'TrajectoryMeasured', 'TraRecoMode', 'Trigger', 'UseOrigin', 'USERNAME', 'VoxPulse1Enum', 'VoxPulse2Enum', 'VoxPulse3Enum', 'WeightingMode', 'WSPulseEnum', 'YesNoMinEchoTime', ]): regex_comp,
frozenset([
'ACQ_grad_matrix', 'ACQ_trim', 'PVM_ChInterPulseDelay', 'PVM_ChSpoilerOnDuration', 'PVM_ChSpoilerStrength', 'PVM_DwBMat', 'PVM_DwDir', 'PVM_DwGradVec', 'PVM_DwSpDir', 'PVM_EpiPhaseCorrection', 'PVM_FovSatSliceVec', 'PVM_RefScanPC', 'PVM_ShimCoeffCoeff', 'PVM_ShimCoeffCoeffSD', 'PVM_SPackArrGradOrient', 'PVM_VoxArrGradOrient', 'PVM_VoxArrPosition', 'PVM_VoxArrSize', 'PVM_VpInterPulseDelay', 'PVM_VpSpoilerOnDuration', 'PVM_VpSpoilerStrength', 'SpoilerDuration', 'SpoilerStrengthArr', ]): arr_comp,
default_param: int_comp}





def list_param_by_type(type_=str, exp_names=['NecS3'], fname = None):
    '''
    create part of the soure code in this module:
    sort JCAMP parameter names into types (str, array and other) and create a
    dictionary that associate frozensets of these parameters with the
    corresponding compare function
    '''

    set_dict = {}
    for typ in [float, numpy.ndarray, str, list, int]:
        set_dict[typ]=set([])

    exp_all = None
    for x in exp_names:
        exp_i = sarpy.Experiment(x)
        if exp_all is None:
            exp_all = exp_i
        else:
            for stdy in exp_i.studies:
                exp_all.add_study(study=stdy)

    for stdy in exp_all.studies:
        for scn in stdy.scans:
            print(scn)
            for t in scn.acqp.__dict__:
                set_dict[type(scn.acqp[t])].add(t)
            for t in scn.method.__dict__:
                set_dict[type(scn.method[t])].add(t)

    if fname:
        with open(fname, 'w') as f:
            print('{frozenset([', file=f)
            for a in sarpy.natural_sort(set_dict[str]):
                print("'%s', " % a, file=f, end='')
            print("]): regex_comp,", file=f)
            print('frozenset([', file=f)
            for a in sarpy.natural_sort(set_dict[numpy.ndarray]):
                print("'%s', " % a, file=f, end='')
            print("]): arr_comp,", file=f)
            print("default_param: int_comp}", file=f)
    else:
        return {frozenset(set_dict[str]): regex_comp,
                frozenset(set_dict[numpy.ndarray]): arr_comp,
                default_param: int_comp,
               }

if __name__ == "__main__":
    list_param_by_type(exp_names=['readfid', 'NecS3Hs15', 'DiLL', 'HPG'],
                             fname = 'JCAMP_dict.py')