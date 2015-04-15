#include "type.h"
#include <stdio.h>
#include <math.h>


void
InitMelFBwindows(MelFB_Window * FirstWin, float StFreq, float SmplFreq, int FFTLength, int NumChannels, BOOLEAN normalize)
{
  int i, j, k;
  int *centrFreq;
  float freq;
  float start_mel;
  float fs_per_2_mel;
  float normFBconst;
  MelFB_Window *p1, *p2;

  centrFreq = (int *) malloc(sizeof(int) * NumChannels);

  /* Constants for calculation */
  start_mel = 2595.0 * log10(1.0 + (float) StFreq / 700.0);
  fs_per_2_mel = 2595.0 * log10(1.0 + (SmplFreq / 2) / 700.0);

  for (i = 0; i < NumChannels; i++) {
    freq = 700 * (pow(10, (start_mel + (float) i / (NumChannels - 1) * (fs_per_2_mel - start_mel)) / 2595.0) - 1.0);
    centrFreq[i] = (int) (FFTLength * freq / SmplFreq + 0.5);
  }
  p1 = FirstWin;

  /*----------------------
   * first freq. window 0
   *----------------------*/
  p1->StartingPoint = centrFreq[0];
  p1->Length = centrFreq[1] - centrFreq[0];
  p1->Data = (float *) malloc(sizeof(float) * p1->Length);
  if (p1->Data == NULL) {
    fprintf(stderr, "ERROR:   Memory allocation error occured!\r\n");
    exit(0);
  }
  normFBconst = 0.0;
  for (j = 0; j < p1->Length; j++) {
    p1->Data[j] = 1.0 - (float) j / (float) p1->Length;
    normFBconst += p1->Data[j];
  }
  if (normalize)
    for (j = 0; j < p1->Length; j++)
      p1->Data[j] /= normFBconst;
  p2 = (MelFB_Window *) malloc(sizeof(MelFB_Window));
  if (p2 == NULL) {
    fprintf(stderr, "ERROR:   Memory allocation error occured!\r\n");
    exit(0);
  }
  p1->Next = p2;
  p1 = p2;

  /*----------------------------------
   * freq. windows 1->NumChannels - 2
   *----------------------------------*/
  for (i = 1; i < NumChannels - 1; i++) {
    p1->StartingPoint = centrFreq[i - 1] + 1;
    p1->Length = centrFreq[i + 1] - centrFreq[i - 1] - 1;
    p1->Data = (float *) malloc(sizeof(float) * p1->Length);
    if (p1->Data == NULL) {
      fprintf(stderr, "ERROR:   Memory allocation error occured!\r\n");
      exit(0);
    }
    normFBconst = 0.0;
    for (j = 0; j < (centrFreq[i] - centrFreq[i - 1]); j++) {
      p1->Data[j] = (float) (j + 1) / (float) (centrFreq[i] - centrFreq[i - 1]);
      normFBconst += p1->Data[j];
    }
    for (j = (centrFreq[i] - centrFreq[i - 1]), k = 0; j < p1->Length; j++, k++) {
      p1->Data[j] = 1.0 - (k + 1) / (float) (centrFreq[i + 1] - centrFreq[i]);
      normFBconst += p1->Data[j];
    }
    if (normalize)
      for (j = 0; j < p1->Length; j++)
	p1->Data[j] /= normFBconst;

    p2 = (MelFB_Window *) malloc(sizeof(MelFB_Window));
    if (p2 == NULL) {
      fprintf(stderr, "ERROR:   Memory allocation error occured!\r\n");
      exit(0);
    }
    p1->Next = p2;
    p1 = p2;
  }

  /*-----------------------------------
   * last freq. window NumChannels - 1
   *-----------------------------------*/
  p1->StartingPoint = centrFreq[NumChannels - 2] + 1;
  p1->Length = centrFreq[NumChannels - 1] - centrFreq[NumChannels - 2];
  p1->Data = (float *) malloc(sizeof(float) * p1->Length);
  normFBconst = 0.0;
  for (j = 0; j < p1->Length; j++) {
    p1->Data[j] = (float) (j + 1) / (float) p1->Length;
    normFBconst += p1->Data[j];
  }
  if (normalize)
    for (j = 0; j < p1->Length; j++)
      p1->Data[j] /= normFBconst;

  p1->Next = NULL;

  free(centrFreq);

  return;
}


void
InitMelIDCTbasis(float **melIDCTbasis, MelFB_Window * FirstWin, short melorder, int sampFreq, int FFTLength)
{
  MelFB_Window *p1;
  float centrFreq[WF_MEL_ORDER];
  float deltaFreq[WF_MEL_ORDER];
  float startingFreq;
  float sum;
  float linStep;
  int i, j;

  linStep = sampFreq / (float) FFTLength;

  /*-------------------------------------
   * calculating centrFreq and deltaFreq
   *-------------------------------------*/
  p1 = FirstWin;
  for (j = 0; j < melorder; j++) {
    if (j == 0) {		/* freq. window 0 */
      centrFreq[j] = p1->StartingPoint * linStep;
      p1 = p1->Next;
    } else {
      if (j == (melorder - 1)) {/* freq. window melorder- */
	centrFreq[j] = (p1->StartingPoint + p1->Length - 1) * linStep;
      } else {			/* freq. windows 1->(melorder-2) */
	startingFreq = p1->StartingPoint * linStep;
	sum = 0.0;
	centrFreq[j] = 0.0;
	for (i = 0; i < p1->Length; i++) {
	  centrFreq[j] += p1->Data[i] * (startingFreq + i * linStep);
	  sum += p1->Data[i];
	}
	centrFreq[j] /= sum;
	p1 = p1->Next;
      }
    }
  }
  for (j = 0; j < melorder; j++) {	/* first and last deltaFreq is only half */
    if (j == 0)
      deltaFreq[j] = (centrFreq[j + 1] - centrFreq[j]) / sampFreq;
    else {
      if (j == (melorder - 1))
	deltaFreq[j] = (centrFreq[j] - centrFreq[j - 1]) / sampFreq;
      else
	deltaFreq[j] = (centrFreq[j + 1] - centrFreq[j - 1]) / sampFreq;
    }
  }

  /*------
   * IDCT
   *------*/
  for (i = 0; i < melorder; i++)/* time axis */
    for (j = 0; j < melorder; j++)	/* frequency axis */
      melIDCTbasis[i][j] = deltaFreq[j] * cos(PIx2 * i * centrFreq[j] / sampFreq);
}



MelFB_Window *
CMelFBAlloc()
{
  MelFB_Window *FirstWin = calloc(1, sizeof(MelFB_Window));

  if (FirstWin == NULL)
    return NULL;
  else
    return FirstWin;
}



void NoiseReInit(NoiseSupStructX *NSX)
{
  X_INT16 i, j;

  NSX->nsVar.SampFreq = 8000;


  /*---------
   * buffers
   *---------*/
  for (i=0 ; i<NS_BUFFER_SIZE ; i++)
    {
      NSX->nsVar.buffers.FirstStageInFloatBuffer [i]  = 0.0;
      NSX->nsVar.buffers.SecondStageInFloatBuffer [i] = 0.0;
    }

  NSX->nsVar.buffers.nbFramesInFirstStage   = 0;
  NSX->nsVar.buffers.nbFramesInSecondStage  = 0;
  NSX->nsVar.buffers.nbFramesOutSecondStage = 0;

  /*-----------
   * dc_filter
   *-----------*/
  NSX->nsVar.prevSamples.lastSampleIn = 0.0;
  NSX->nsVar.prevSamples.lastDCOut    = 0.0;

  /*-----------
   * gain_fact
   *-----------*/
  for (i=0 ; i<3 ; i++) NSX->nsVar.gainFact.denEn1 [i] = 0.0;
  NSX->nsVar.gainFact.lowSNRtrack = 0.0;
  NSX->nsVar.gainFact.alfaGF = 0.8;

    /*-----------
   * agc_data_ns
   *-----------*/
  NSX->nsVar.agc_data_ns.agc_gain = 0.0;
  NSX->nsVar.agc_data_ns.Last_Spower = 0.0;
  NSX->nsVar.agc_data_ns.Spower = 0.0;
  NSX->nsVar.agc_data_ns.nb_adapt = 0;
  /*-------------
   * vad_data_ns
   *-------------*/
  NSX->nsVar.vadNS.nbFrame [0]    = NSX->nsVar.vadNS.nbFrame [1] = 0;
  NSX->nsVar.vadNS.hangOver       = 0;
  NSX->nsVar.vadNS.nbSpeechFrames = 0;
  NSX->nsVar.vadNS.meanEn         = 0.0;
  NSX->nsVar.vadNS.flagVAD        = 0;

  NSX->nsVar.vadCA.hangOver       = 0;
  NSX->nsVar.vadCA.nbSpeechFrames = 0;
  NSX->nsVar.vadCA.meanCA         = 0.0;
  NSX->nsVar.vadCA.flagVAD        = 0;

  /*-------------
   * vad_data_fd
   *-------------*/
  NSX->nsVar.vadFD.MelMean = 0.0;
  NSX->nsVar.vadFD.VarMean = 0.0;
  NSX->nsVar.vadFD.AccTest = 0.0;
  NSX->nsVar.vadFD.AccTest2 = 0.0;
  NSX->nsVar.vadFD.SpecMean = 0.0;
  NSX->nsVar.vadFD.MelValues[0] = NSX->nsVar.vadFD.MelValues[1] = 0.0;
  NSX->nsVar.vadFD.SpecValues = 0.0;
  NSX->nsVar.vadFD.SpeechInVADQ = 0.0;
  NSX->nsVar.vadFD.SpeechInVADQ2 = 0.0;

  /*----------
   * spectrum
   *----------*/
  for (i=0 ; i<NS_SPEC_ORDER ; i++)
    {
      NSX->nsVar.spectrum.denSigSE1[i]   = NSX->nsVar.spectrum.denSigSE2[i]   = 0.0;
      NSX->nsVar.spectrum.nSigSE1[i]     = NSX->nsVar.spectrum.nSigSE2[i]     = 0.0;
      NSX->nsVar.spectrum.noiseSE1[i]    = NSX->nsVar.spectrum.noiseSE2[i]    = NS_EPS;
    }

  for (i=0 ; i<NS_SPEC_ORDER; i++)
    {
      for (j=0 ; j<NS_PSD_MEAN_ORDER ; j++)
	NSX->nsVar.spectrum.PSDMeanBuffer1[i][j] = 0.0;
    }

  for (i=0 ; i<NS_SPEC_ORDER ; i++)
    {
      for (j=0 ; j<NS_PSD_MEAN_ORDER ; j++)
	NSX->nsVar.spectrum.PSDMeanBuffer2[i][j] = 0.0;
    }

  NSX->nsVar.spectrum.indexBuffer1 = 0;
  NSX->nsVar.spectrum.indexBuffer2 = 0;

  /*--------
   * ns_tmp
   *--------*/
  // Hanning window for spectrum estimation
  for (i=0 ; i<NS_FRAME_LENGTH ; i++)
    NSX->nsTmp.sigWindow[i] = 0.5 - 0.5 * cos ((PIx2 * ((X_FLOAT32) i + 0.5)) / (X_FLOAT32) (NS_FRAME_LENGTH));

  // Hanning window for impulse response windowing
  for (i=0 ; i<NS_FILTER_LENGTH ; i++)
    NSX->nsTmp.IRWindow[i] = 0.5 - 0.5 * cos ((PIx2 * ((X_FLOAT32) i + 0.5)) / (X_FLOAT32) (NS_FILTER_LENGTH));

  // mel FB windows
  NSX->nsTmp.FirstWindow = CMelFBAlloc ();
  if (NSX->nsTmp.FirstWindow == NULL)
    {
      fprintf (stderr, "ERROR:   Memory allocation error occured!\r\n");
      exit(0);
    }
  InitMelFBwindows (NSX->nsTmp.FirstWindow, 0.0, (X_FLOAT32)NSX->nsVar.SampFreq, 2*(NS_SPEC_ORDER-1), WF_MEL_ORDER, 1);

  // mel IDCT
  NSX->nsTmp.melIDCTbasis = (X_FLOAT32 **) malloc (sizeof (X_FLOAT32 *) * WF_MEL_ORDER);
  if (NSX->nsTmp.melIDCTbasis == NULL)
    {
      fprintf (stderr, "ERROR:   Memory allocation error occured!\r\n");
      exit(0);
    }

  for (i=0 ; i<WF_MEL_ORDER ; i++)
    {
      NSX->nsTmp.melIDCTbasis[i] = (X_FLOAT32 *) malloc (sizeof (X_FLOAT32) * WF_MEL_ORDER);
      if (NSX->nsTmp.melIDCTbasis[i] == NULL)
	{
	  fprintf (stderr, "ERROR:   Memory allocation error occured!\r\n");
	  exit(0);
	}
    }

  InitMelIDCTbasis (NSX->nsTmp.melIDCTbasis, NSX->nsTmp.FirstWindow, WF_MEL_ORDER, NSX->nsVar.SampFreq, 2*(NS_SPEC_ORDER-1));
}


void logMMSE_Init(LOGMMSE_VAR *logMMSEVar)
{
  X_INT16 i, j;

  logMMSEVar->SampFreq = samplRate;
  logMMSEVar->nbFrame  = 0;

  for (i=0 ; i < FFT_LEN ; i++)
    {
      logMMSEVar->Buf2denoise[i]  = 0.0;
      logMMSEVar->FFTIn[i].real  = 0.0;
      logMMSEVar->FFTIn[i].imag  = 0.0;
      logMMSEVar->spect[i]       = 0.0;
      logMMSEVar->noiseMean[i]   = 0.0;
      logMMSEVar->Xk_prev[i]     = 0.0;
      logMMSEVar->sig[i]         = 0.0;
      logMMSEVar->gammak[i]      = 0.0;
      logMMSEVar->ksi[i]         = 0.0;
      logMMSEVar->log_sigma_k[i] = 0.0;
      logMMSEVar->A[i]           = 0.0;
      logMMSEVar->vk[i]          = 0.0;
      logMMSEVar->ei_vk[i]       = 0.0;
      logMMSEVar->hw[i]          = 0.0;
    }

  for (i=0 ; i < FRAME_SHIFT ; i++)
    logMMSEVar->old_denosebuf[i] = 0.0;
    /*-----------
   * agc_data_ns
   *-----------*/
  logMMSEVar->agc_data_ns.agc_gain = 0.0;
  logMMSEVar->agc_data_ns.Last_Spower = 0.0;
  logMMSEVar->agc_data_ns.Spower = 0.0;
  logMMSEVar->agc_data_ns.nb_adapt = 0;
  /*-------------
   * vad_data_ns
   *-------------*/
  logMMSEVar->vadNS.nbFrame [0]    =logMMSEVar->vadNS.nbFrame [1] = 0;
  logMMSEVar->vadNS.hangOver       = 0;
  logMMSEVar->vadNS.nbSpeechFrames = 0;
  logMMSEVar->vadNS.meanEn         = 0.0;
  logMMSEVar->vadNS.flagVAD        = 0;



  /*--------
   * ns_tmp
   *--------*/
  // Hamming window for spectrum estimation
  for (i = 0 ; i < FRAME_LEN ; i++)
    logMMSEVar->sigWindow[i] = 0.54 - 0.46 * cos ((PIx2 *  i ) / (X_FLOAT32) (FRAME_LEN));

}

