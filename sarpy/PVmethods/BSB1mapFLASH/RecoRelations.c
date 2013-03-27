/****************************************************************
 *
 * $Source: /pv/CvsTree/pv/gen/src/prg/methods/FLASH/RecoRelations.c,v $
 *
 * Copyright (c) 2002
 * Bruker BioSpin MRI GmbH
 * D-76275 Ettlingen, Germany
 *
 * All Rights Reserved
 *
 * $Id: RecoRelations.c,v 1.15.2.4 2009/08/14 14:25:05 sako Exp $
 *
 ****************************************************************/


static const char resid[] = "$Id: RecoRelations.c,v 1.15.2.4 2009/08/14 14:25:05 sako Exp $ (C) 2002 Bruker BioSpin MRI GmbH";

#define DEBUG           0
#define DB_MODULE       0
#define DB_LINE_NR      0



#include "method.h"

void SetRecoParam( void )

{
  
  int dim,i,movieframes,size;
  
  DB_MSG(("-->SetRecoParam\n"));
  
  /* set baselevel reconstruction parameter */
  /* default initialization of reco based on acqp pars allready set */
  
  ATB_InitDefaultReco();
  
  /* select method specific reconstruction method */
  SetNewRecoParam();

 
  if(PVM_MovieOnOff == On)
  {
    movieframes = PVM_NMovieFrames;
  }
  else
  {
    movieframes =1; 
  }  

  /* set reco rotate according to phase offsets     */

  dim = PTB_GetSpatDim();

  /* set reco sizes and ft_mode (dim 2&3) */
  /* (dim 1 is kept as it was set by ATB_InitDefaultReco) */
  for(i=0; i<dim; i++)
  {
    size = (int)(PVM_Matrix[i]*PVM_AntiAlias[i]);
    RECO_ft_mode[i] = (size == PowerOfTwo(size)) ?  COMPLEX_FFT:COMPLEX_FT;
    RECO_ft_size[i] = size;
    RECO_size[i] = PVM_Matrix[i];
  }  
 
  ParxRelsParRelations("RECO_ft_mode",Yes);
  ParxRelsParRelations("RECO_ft_size",Yes);
  ParxRelsParRelations("RECO_size",Yes);

  
  ATB_SetRecoRotate(PVM_EffPhase1Offset,
                    PVM_Fov[1]*PVM_AntiAlias[1],
                    NSLICES,     
                    movieframes,           /* 1 echo image    */
                    1) ;         /* phase1 direction*/
  
  if(dim==3)
  {
    
    ATB_SetRecoRotate(PVM_EffPhase2Offset, 
		      PVM_Fov[2]*PVM_AntiAlias[2],
		      NSLICES,     
		      movieframes,           /* 1 echo image    */
		      2) ;         /* phase2 direction*/
    
  }
  
  
  /* set reco offset */
  
  ATB_SetRecoOffset(RECO_ft_size,
                    PVM_AntiAlias,
                    NI,              
                    dim);
  for(i=0;i<dim;i++)
    RECO_fov[i]= PVM_FovCm[i];
  
  ParxRelsParRelations("RECO_fov",Yes);
  
  ATB_SetRecoTranspositionFromLoops(PtrType3x3 ACQ_grad_matrix[0],
				    NSLICES,
				    1, /* no loops within slice loop */
				    NI,
				    ACQ_obj_order);
   

  
  DB_MSG(("<--SetRecoParam\n"));
  
}

void SetNewRecoParam(void)
{
  DB_MSG(("-->SetNewRecoParam\n"));

  if(PVM_EncUseMultiRec == Yes || PVM_EncPftAccel1 > 1.0 || RecoMethMode==SWI)
  {
    int i, k, ftSize[3];

    RECO_mode = USER_MODE;
    ParxRelsParRelations("RECO_mode",Yes);

    for(i=0; i<PTB_GetSpatDim(); i++)
    {
      ftSize[i] = (int)(PVM_Matrix[i]*PVM_AntiAlias[i]);  
      RECO_ft_mode[i] = (ftSize[i] == PowerOfTwo(ftSize[i])) ?  COMPLEX_FFT:COMPLEX_FT;
      RECO_ft_size[i] = ftSize[i];
      RECO_size[i] = PVM_Matrix[i];
    }
  
  
    ATB_InitUserModeReco(ACQ_dim,
			 PVM_EncMatrix,
			 ftSize,
			 PVM_EncSteps1,
			 PVM_EncSteps2,
			 PVM_EncNReceivers,
			 PVM_EncPpiAccel1,
			 PVM_EncPpiRefLines1,
			 NI,
			 ACQ_obj_order,
			 ACQ_phase_factor,
			 PVM_EchoPosition);
  
  
    /* set scaling values for phased array coils */
    for(k=0; k<PVM_EncNReceivers;k++)
      RecoScaleChan[k] = PVM_EncChanScaling[k];
    
  }
  else
    RECO_mode = FT_MODE;
  
  DB_MSG(("<--SetNewRecoParam\n"));
}

void RecoUserUpdateRel (void)
{
  DB_MSG(("-->RecoUserUpdateRel\n"));
  
  if(RecoMethMode==SWI)
    SWIRecoRel(); 
  else
    ParxRelsParRelations("RecoUserUpdate",Yes);
  
  DB_MSG(("<--RecoUserUpdateRel\n"));
}

void SWIRecoRel (void)
{
  DB_MSG(("-->SWIRecoRel"));

  if (RecoUserUpdate == Yes)
  {
    int i, n, dim;
    char arg1[128], arg2[128], arg3[128], arg4[128];
    int nrReceivers;
    double delta[3];
    double tau[3];
    double gamma[3];    
    
    dim = PTB_GetSpatDim();
    nrReceivers = RecoNrActiveReceivers(); /* GetNrActiveReceivers */

    if(WeightingMode==phase_image)
    {
      RecoCombineMode = AddImages;
      RECO_image_type = PHASE_IMAGE;
    }
    else
    {
      RecoCombineMode = SumOfSquares;
      RECO_image_type = MAGNITUDE_IMAGE;  
    }
    ParxRelsParRelations("RECO_image_type",Yes);

    /* parameters for RecoGaussWinMultFilter: exp(-delta*t - tau*t^2 -gamma) 
       tau = (N*pi*broadening/(2*FOV))^2
       delta = -2*tau*t0 
       gamma = tau * t0^2 */
 
    for(i=0;i<dim;i++)
    {
      tau[i] = pow((PVM_Matrix[i]*M_PI*GaussBroadening/(2*PVM_Fov[i])),2);  
    }

    delta[0]=-2*tau[0]*(PVM_EchoPosition/100)/PVM_EncZfRead;
    delta[1]=-2*tau[1]*0.5;
    delta[2]=-2*tau[2]*0.5;

    gamma[0]=tau[0]*pow(((PVM_EchoPosition/100)/PVM_EncZfRead),2);
    gamma[1]=tau[1]*pow(0.5,2);
    gamma[2]=tau[2]*pow(0.5,2);

    sprintf(arg1,"winDirection=0;delta=%f;tau=%f;gamma=%f",delta[0],tau[0],gamma[0]);
    sprintf(arg2,"winDirection=1;delta=%f;tau=%f;gamma=%f",delta[1],tau[1],gamma[1]);
    sprintf(arg3,"winDirection=2;delta=%f;tau=%f;gamma=%f",delta[2],tau[2],gamma[2]);

    if(WeightingMode == positive_mask)
      sprintf(arg4,"mask=1;weighting=%f",MaskWeighting);
    else if (WeightingMode == negative_mask)
      sprintf(arg4,"mask=2;weighting=%f",MaskWeighting);      
    else if (WeightingMode == magnitude_mask)
      sprintf(arg4,"mask=3;weighting=%f",MaskWeighting);  
    else if (WeightingMode == phase_image)
      sprintf(arg4,"mask=4;weighting=%f",MaskWeighting);

    RecoDeriveInputProcess();
    RecoDeriveComputeProcess();
  
    if (ACQ_scan_type == Scan_Experiment)
    {
      /* apply gauss filtering and phase weighting to data */ 
      if(dim==2)
      {  
	for (n = 0; n < nrReceivers; n++)
	{
	  RecoComputeAppendStage("PREPPASS",n,"Q","RecoTeeFilter","TEE","");
	  RecoComputeAppendStage("PREPPASS",n,"TEE","RecoPhaseWeightingFilter","SWI",arg4);
	  RecoComputeAddStage("PREPPASS",n,"RecoFTFilter","FTINV1","direction=0;exponent=1");
	  RecoComputeAddStage("PREPPASS",n,"RecoFTFilter","FTINV2","direction=1;exponent=1");
	  RecoComputeAddStage("PREPPASS",n,"RecoGaussWinMultFilter","GAUSS1",arg1);	
	  RecoComputeAddStage("PREPPASS",n,"RecoGaussWinMultFilter","GAUSS2",arg2);
	  RecoComputeAddStage("PREPPASS",n,"RecoFTFilter","FTINV3","direction=0;exponent=-1");
	  RecoComputeAddStage("PREPPASS",n,"RecoFTFilter","FTINV4","direction=1;exponent=-1");
	  RecoComputeConnectStages("PREPPASS",n,"TEE","FTINV1");
	  RecoComputeConnectStages("PREPPASS",n,"FTINV1","FTINV2");
	  RecoComputeConnectStages("PREPPASS",n,"FTINV2","GAUSS1");
	  RecoComputeConnectStages("PREPPASS",n,"GAUSS1","GAUSS2");
	  RecoComputeConnectStages("PREPPASS",n,"GAUSS2","FTINV3");
	  RecoComputeConnectStages("PREPPASS",n,"FTINV3","FTINV4");
	  RecoComputeConnectStages("PREPPASS",n,"FTINV4","SWI.mask");
	}
      }
      
      if(dim == 3)
      {
	for (n = 0; n < nrReceivers; n++)
	{
	  RecoComputeAppendStage("PREPPASS",n,"Q","RecoTeeFilter","TEE","");
	  RecoComputeAppendStage("PREPPASS",n,"TEE","RecoPhaseWeightingFilter","SWI",arg4);
	  RecoComputeAddStage("PREPPASS",n,"RecoFTFilter","FTINV1","direction=0;exponent=1");
	  RecoComputeAddStage("PREPPASS",n,"RecoFTFilter","FTINV2","direction=1;exponent=1");
	  if(PVM_EncPpiAccel1>1)
	    RecoComputeAddStage("PREPPASS",n,"RecoFTFilter","FTINV3","direction=2;exponent=-1");
	  else
	    RecoComputeAddStage("PREPPASS",n,"RecoFTFilter","FTINV3","direction=2;exponent=1");
	  RecoComputeAddStage("PREPPASS",n,"RecoGaussWinMultFilter","GAUSS1",arg1);	
	  RecoComputeAddStage("PREPPASS",n,"RecoGaussWinMultFilter","GAUSS2",arg2);
	  RecoComputeAddStage("PREPPASS",n,"RecoGaussWinMultFilter","GAUSS3",arg3);
	  RecoComputeAddStage("PREPPASS",n,"RecoFTFilter","FTINV4","direction=0;exponent=-1");
	  RecoComputeAddStage("PREPPASS",n,"RecoFTFilter","FTINV5","direction=1;exponent=-1");
	  if(PVM_EncPpiAccel1>1)
	    RecoComputeAddStage("PREPPASS",n,"RecoFTFilter","FTINV6","direction=2;exponent=1");
	  else
	    RecoComputeAddStage("PREPPASS",n,"RecoFTFilter","FTINV6","direction=2;exponent=-1");
	  RecoComputeConnectStages("PREPPASS",n,"TEE","FTINV1");
	  RecoComputeConnectStages("PREPPASS",n,"FTINV1","FTINV2");
	  RecoComputeConnectStages("PREPPASS",n,"FTINV2","FTINV3");
	  RecoComputeConnectStages("PREPPASS",n,"FTINV3","GAUSS1");
	  RecoComputeConnectStages("PREPPASS",n,"GAUSS1","GAUSS2");
	  RecoComputeConnectStages("PREPPASS",n,"GAUSS2","GAUSS3");
	  RecoComputeConnectStages("PREPPASS",n,"GAUSS3","FTINV4");
	  RecoComputeConnectStages("PREPPASS",n,"FTINV4","FTINV5");
	  RecoComputeConnectStages("PREPPASS",n,"FTINV5","FTINV6");
	  RecoComputeConnectStages("PREPPASS",n,"FTINV6","SWI.mask");
	}
      }
    }

    RecoDeriveOutputProcess();
  }

  DB_MSG(("<--SWIRecoRel"));
}

int PowerOfTwo(int x)
/* returns next power of two */
{

 return (1<<(int)ceil(log((double)x)/log(2.0)));
}
