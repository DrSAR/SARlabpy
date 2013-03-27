/*
 *******************************************************************
 *
 * $Source: /bscl/CvsTree/bscl/gen/config/proto.head,v $
 *
 * Copyright (c) 1995
 * BRUKER ANALYTISCHE MESSTECHNIK GMBH
 * D-76287 Rheinstetten, Germany
 *
 * All Rights Reserved
 *
 *
 * $State: Exp $
 *
 *******************************************************************
 */

#ifndef _P_
#	if defined(HAS_PROTO) || defined(__STDC__) || defined(__cplusplus)
#		define _P_(s) s
#	else
#		define _P_(s) ()
#	endif
#endif

/* /opt/PV5.1/prog/parx/src/BSB1mapFLASH/initMeth.c */
void initMeth _P_((void));
/* /opt/PV5.1/prog/parx/src/BSB1mapFLASH/loadMeth.c */
void loadMeth _P_((const char *));
/* /opt/PV5.1/prog/parx/src/BSB1mapFLASH/parsRelations.c */
void backbone _P_((void));
void dsRelations _P_((void));
void dsRange _P_((void));
void ExcPulseAngleRelation _P_((void));
void ExcPulseEnumRelation _P_((void));
void ExcPulseRelation _P_((void));
void ExcPulseRange _P_((void));
void echoTimeRels _P_((void));
void SliceSegDurRels _P_((void));
double minLoopRepetitionTime _P_((void));
void repetitionTimeRels _P_((void));
void LocalGeometryMinimaRels _P_((double *, double *));
void LocalGradientStrengthRels _P_((void));
void LocalFrequencyOffsetRels _P_((void));
void Local_NAveragesRange _P_((void));
void Local_NAveragesHandler _P_((void));
void localHandleSliceGeometry _P_((void));
void LocalEchoTime1Relation _P_((void));
void Local_NMovieFramesRange _P_((void));
void Local_MovieOnOffRange _P_((void));
void Local_NMovieFramesRels _P_((void));
void Local_MovieOnOffRels _P_((void));
void Local_UpdateMovie _P_((void));
void EchoTimeModeRels _P_((void));
void EchoTimeModeRange _P_((void));
void EffSWhRange _P_((void));
void EffSWhRel _P_((void));
void localInversionRel _P_((void));
void ControlGradientLimits _P_((YesNo));
void RecoMethModeRel _P_((void));
void RecoMethModeVisPar _P_((void));
void MaskModeRel _P_((void));
void GaussBroadRange _P_((void));
void MaskWeightRange _P_((void));
void BSPulseEnumRelation _P_((void));
void BSPulseRelation _P_((void));
void BSPulseRange _P_((void));
/* /opt/PV5.1/prog/parx/src/BSB1mapFLASH/BaseLevelRelations.c */
void SetBaseLevelParam _P_((void));
void SetBasicParameters _P_((void));
void SetFrequencyParameters _P_((void));
void SetGradientParameters _P_((void));
void SetInfoParameters _P_((void));
void SetMachineParameters _P_((void));
void SetPpgParameters _P_((void));
void SetACQ_obj_orderForMovie _P_((void));
/* /opt/PV5.1/prog/parx/src/BSB1mapFLASH/RecoRelations.c */
void SetRecoParam _P_((void));
void SetNewRecoParam _P_((void));
void RecoUserUpdateRel _P_((void));
void SWIRecoRel _P_((void));
int PowerOfTwo _P_((int));
/* /opt/PV5.1/prog/parx/src/BSB1mapFLASH/deriveVisu.c */
void deriveVisu _P_((void));
