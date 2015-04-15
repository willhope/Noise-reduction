#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "type.h"
#include "rfft.h"
#define max(a,b) ((a>b)?(a):(b))



static void DoSigWindowing (X_FLOAT32 *Data, X_FLOAT32 *window, X_INT16 frameLength, X_INT16 FFTLength)
{
  X_INT16 i;

  // windowing
  for (i=0 ; i<frameLength ; i++)
    Data [i] = Data[i] * window [i];

  // zero padding
  for (i=frameLength ; i<FFTLength ; i++)
    Data [i] = 0.0;

  return;
}


static void FFTtoPSD (X_FLOAT32 *FFTIn, X_FLOAT32 *PSDOut, X_INT16 FFTLength)
{
  X_INT16 i, j;

  FFTIn[0] = FFTIn[0] * FFTIn[0];

  for (i=1,j=FFTLength-1; i<FFTLength/2 ; i++,j--)
    {
      FFTIn[i] = (FFTIn[i] * FFTIn[i] + FFTIn[j] * FFTIn[j]);
    }

  FFTIn[i] = FFTIn[i] * FFTIn[i];

  // from 129 to 65 FB
  for (i=0,j=0 ; i<NS_SPEC_ORDER-1 ; i++,j+=2)
    {
      PSDOut[i] = (FFTIn[j] + FFTIn[j+1]) / 2.0;
    }
  PSDOut[i] = FFTIn[j];

  return;
}


static void PSDMean (X_INT16 *indexBuffer, X_FLOAT32 *PSDIn, X_FLOAT32 *PSDOut, X_FLOAT32 *PSDBuffer)
{
  X_INT16 i, index;

  *indexBuffer = 1 - *indexBuffer;

  for (i=0 ; i<NS_SPEC_ORDER ; i++)
    {
      index = i * NS_PSD_MEAN_ORDER + *indexBuffer;
      PSDBuffer[index] = PSDIn[i];

      index += 1 - 2 * *indexBuffer;
      PSDOut[i] = (PSDBuffer[index] + PSDIn[i]) / NS_PSD_MEAN_ORDER;
    }
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: ApplyWF
 *
 * PURPOSE:       Convolving input data with denoised filter impulse response
 *
 * INPUT:
 *   *data        Pointer to input noisy signal
 *   *predata     Pointer to input noisy signal
 *   *filter      Pointer to denoised filter impulse response
 *   frameShift   Frame shift
 *   melOrder     Impulse response length
 *
 * OUTPUT:
 *   result       Output denoised data
 *
 * RETURN VALUE:
 *   none
 *
 *---------------------------------------------------------------------------*/
static void ApplyWF (X_FLOAT32 *data, X_FLOAT32 *predata, X_FLOAT32 *filter, X_FLOAT32 *result, X_INT16 frameShift, X_INT16 melOrder)
{
  X_INT16 i, j;

  for (i=0 ; i<frameShift ; i++)
    {
      result[i] = 0.0;
      for (j=-melOrder ; j<=(i<=melOrder ? i : melOrder) ; j++)
	{
	  result[i] += (filter[j + melOrder] * data[i-j]);
	}
      for (j=i+1 ; j<=melOrder ; j++)
	{
	  result[i] += (filter[j + melOrder] * predata[frameShift - j + i]);
	}
    }
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: VAD
 *
 * PURPOSE:       Voice Ativity Detection
 *
 * INPUT:
 *   fstage       Stage of two stage noise suppression
 *   *vadNS       Pointer to VAD structure
 *   *newShiftFrame Pointer to new input frame
 *
 * OUTPUT:        VAD for noise suppression structure updated
 *
 *
 * RETURN VALUE:
 *   none
 *
 *---------------------------------------------------------------------------*/
static void VAD (X_INT16 fstage, VAD_DATA_NS *vadNS, const X_FLOAT32 *newShiftFrame)
{
  X_INT16 i;
  X_INT16 flagVAD;
  X_INT16 hangOver;
  X_INT16 nbSpeechFrames;
  X_INT32 nbFrame;
  X_FLOAT32 meanEn;
  X_FLOAT32 frameEn;
  X_FLOAT32 lambdaLTE;

  nbFrame        = vadNS->nbFrame [fstage];
  meanEn         = vadNS->meanEn;
  flagVAD        = vadNS->flagVAD;
  hangOver       = vadNS->hangOver;
  nbSpeechFrames = vadNS->nbSpeechFrames;

  if (nbFrame < 2147483647) nbFrame++;
  vadNS->nbFrame [fstage] = nbFrame;

  if (fstage == 1) return;

  if (nbFrame < NS_NB_FRAME_THRESHOLD_LTE)
    lambdaLTE = 1 - 1 / (X_FLOAT32) nbFrame;
  else
    lambdaLTE = NS_LAMBDA_LTE_LOWER_E;

  frameEn = 64.0;

  for (i=0 ; i<NS_FRAME_SHIFT ; i++)
    frameEn += newShiftFrame[i] * newShiftFrame[i];

  frameEn = (X_FLOAT32) (0.5 + (log (frameEn / 64.0) / log(2.0)) * 16.0);

  if (((frameEn - meanEn) < NS_SNR_THRESHOLD_UPD_LTE) || (nbFrame < NS_MIN_FRAME))
    {
      if ((frameEn < meanEn) || (nbFrame < NS_MIN_FRAME))
	meanEn += (1 - lambdaLTE) * (frameEn - meanEn);
      else
	meanEn += (1 - NS_LAMBDA_LTE_HIGHER_E) * (frameEn - meanEn);

      if (meanEn < NS_ENERGY_FLOOR)
	meanEn = NS_ENERGY_FLOOR;
    }
  if (nbFrame > 4)
    {
      if ((frameEn - meanEn) > NS_SNR_THRESHOLD_VAD)
	{
	  flagVAD = 1;
	  nbSpeechFrames++;
	}
      else
	{
	  if (nbSpeechFrames > NS_MIN_SPEECH_FRAME_HANGOVER)
	    hangOver = NS_HANGOVER;
	  nbSpeechFrames = 0;
	  if (hangOver != 0)
	    {
	      hangOver--;
	      flagVAD = 1;
	    }
	  else
	    flagVAD = 0;
	}
    }
  vadNS->meanEn         = meanEn;
  vadNS->flagVAD        = flagVAD;
  vadNS->hangOver       = hangOver;
  vadNS->nbSpeechFrames = nbSpeechFrames;

  return;
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: FilterCalc
 *
 * PURPOSE:       Computation of noise suppression filter in the frequency domain
 *
 * INPUT:
 *   fstage       Stage of two stage noise reduction
 *   *NSX         Pointer to noise suppression structure
 *   *PSDMeaned   Pointer to smoothed PSD
 *
 * OUTPUT:
 *   W            Noise suppression filter
 *
 * RETURN VALUE:
 *   none
 *
 *---------------------------------------------------------------------------*/
static void FilterCalc (X_INT16 fstage, NoiseSupStructX *NSX, X_FLOAT32 *PSDMeaned, X_FLOAT32 *W)
{
  NS_VAR    nsVar = NSX->nsVar;

  X_INT16   i;
  X_INT16   flagVAD;
  X_INT16   nbFrame;

  X_FLOAT32 SNRprio;
  X_FLOAT32 SNRpost;
  X_FLOAT32 lambdaNSE;
  X_FLOAT32 *denSigSE;
  X_FLOAT32 *nSigSE;
  X_FLOAT32 *noiseSE;

  nSigSE    = NSX->nsVar.spectrum.nSigSE1;
  noiseSE   = NSX->nsVar.spectrum.noiseSE1;
  denSigSE  = NSX->nsVar.spectrum.denSigSE1;

  if (fstage == 1)
    {
      nSigSE    = NSX->nsVar.spectrum.nSigSE2;
      noiseSE   = NSX->nsVar.spectrum.noiseSE2;
      denSigSE  = NSX->nsVar.spectrum.denSigSE2;
    }

  nbFrame   = nsVar.vadNS.nbFrame [fstage];
  flagVAD   = nsVar.vadNS.flagVAD;

  /*-------------------------------------------------------
   * Choice of the Noise Estimation according to 2WF stage
   * VAD based NE only at 1st stage of 2WF
   *-------------------------------------------------------*/

  /*--------------------------------
   * non VAD based Noise Estimation
   *--------------------------------*/
  if (fstage == 1)
    {
      // noise estimation in energy
      for (i=0 ; i<NS_SPEC_ORDER ; i++)
	noiseSE[i] *= noiseSE[i];

      if (nbFrame < 11)
	{
	  lambdaNSE = 1 - 1 / (X_FLOAT32) nbFrame;
	  for (i=0 ; i<NS_SPEC_ORDER ; i++)
	    {
	      noiseSE[i] = lambdaNSE * noiseSE[i] + (1 - lambdaNSE) * PSDMeaned[i];
	    }
	}
      else
	{
	  X_FLOAT32 upDate;
	  for (i=0 ; i<NS_SPEC_ORDER ; i++)
	    {
	      upDate = 0.9 + 0.1 * (PSDMeaned [i] / (PSDMeaned [i] + noiseSE [i])) *
		(1.0 + 1.0 / (1.0 + 0.1 * (PSDMeaned[i] / noiseSE[i])));
	      noiseSE [i] *= upDate;
	    }
	}

      // store noise estimation values in magnitude
      for (i=0 ; i<NS_SPEC_ORDER ; i++)
	{
	  noiseSE [i] = (X_FLOAT32) sqrt (noiseSE[i]);
	  if (noiseSE[i] < NS_EPS) noiseSE[i] = NS_EPS;
	}
    }

  /*----------------------------------------------------------------------------
   * Spectral estimations of noisy signal are now stored in magnitude in nSigSE
   *----------------------------------------------------------------------------*/
  for (i=0 ; i<NS_SPEC_ORDER ; i++)
    {
      nSigSE[i]    = (X_FLOAT32) sqrt (nSigSE[i]);
      PSDMeaned[i] = (X_FLOAT32) sqrt (PSDMeaned[i]);
    }

  /*-------------------------------------------
   * VAD based Noise Estimation (in magnitude)
   *-------------------------------------------*/
  if (fstage == 0)
    {
      if (nbFrame < NS_NB_FRAME_THRESHOLD_NSE)
	lambdaNSE = 1 - 1 / (X_FLOAT32) nbFrame;
      else
	lambdaNSE =  NS_LAMBDA_NSE;

      if (flagVAD == 0)
	{
	  for (i=0 ; i<NS_SPEC_ORDER ; i++)
	    {
	      noiseSE [i] = lambdaNSE * noiseSE [i] + (1 - lambdaNSE) * PSDMeaned[i];
	      if (noiseSE[i] < NS_EPS) noiseSE[i] = NS_EPS;
	    }
	}
    }

  /*--------------------------------------
   * noise suppression filter calculation
   *--------------------------------------*/
  for (i=0 ; i<NS_SPEC_ORDER ; i++)
    {
      SNRpost = (PSDMeaned[i] / noiseSE[i]) - 1;
      SNRprio = NS_BETA * (denSigSE [i] / noiseSE [i]) + (1 - NS_BETA) * max (0, SNRpost);
      W[i] = SNRprio / (1 + SNRprio);
      SNRprio = W[i] * PSDMeaned[i] / noiseSE[i];
      SNRprio = max (SNRprio, NS_RSB_MIN);
      W[i] = SNRprio / (1 + SNRprio);
      denSigSE [i] = W[i] * nSigSE[i];
    }

  return;
}


void
DoMelFB(float *SigFFT, MelFB_Window * FirstWin)
{
  MelFB_Window *p1;
  float Sum[WF_MEL_ORDER];
  int i, j;

  p1 = FirstWin;
  j = 0;
  while (p1) {
    Sum[j] = 0.0;
    for (i = 0; i < p1->Length; i++)
      Sum[j] += SigFFT[p1->StartingPoint + i] * p1->Data[i];

    j++;
    p1 = p1->Next;
  }

  for (j--; j >= 0; j--)
    SigFFT[j] = Sum[j];

  return;
}


static void DoGainFact (X_INT16 fstage, NoiseSupStructX *NSX, X_FLOAT32 *W)
{
  X_INT16 i;

  X_FLOAT32 averSNR;
  X_FLOAT32 noiseEn;
  X_FLOAT32 lambdaSNR;
  X_FLOAT32 *noiseSE  = NSX->nsVar.spectrum.noiseSE2;
  X_FLOAT32 *denSigSE = NSX->nsVar.spectrum.denSigSE1;

  GAIN_FACT *gainFact = &(NSX->nsVar.gainFact);

  if (fstage == 0)
    {
      gainFact->denEn1 [0] = gainFact->denEn1 [1];
      gainFact->denEn1 [1] = gainFact->denEn1 [2];
      gainFact->denEn1 [2] = 0.0;
      for (i=0 ; i<NS_SPEC_ORDER ; i++) gainFact->denEn1 [2] += denSigSE [i]; // new denEn1
    }
  else
    {
      noiseEn = 0.0;
      for (i=0 ; i<NS_SPEC_ORDER ; i++) noiseEn += noiseSE [i];
      averSNR = (gainFact->denEn1 [0] * gainFact->denEn1 [1] * gainFact->denEn1 [2]) / (noiseEn * noiseEn * noiseEn);

      if (averSNR > 0.00001)
	averSNR = (20 * log10 (averSNR)) / 3.0;
      else
	averSNR = -100.0 / 3.0;

      if ( ((averSNR - gainFact->lowSNRtrack) < 10.0) || (NSX->nsVar.vadNS.nbFrame[fstage] < NS_MIN_FRAME) )
	{
	  if (NSX->nsVar.vadNS.nbFrame[fstage] < NS_MIN_FRAME)
	    lambdaSNR = 1.0 - 1.0 / (X_FLOAT32) NSX->nsVar.vadNS.nbFrame[fstage];
	  else
	    {
	      if (averSNR < gainFact->lowSNRtrack)
		lambdaSNR = 0.95;
	      else
		lambdaSNR = 0.99;
	    }
	  gainFact->lowSNRtrack += (1.0 - lambdaSNR) * (averSNR - gainFact->lowSNRtrack);
	}

      if (gainFact->denEn1 [2] > 100) // no change if very low signal
	{
	  if (averSNR < (gainFact->lowSNRtrack + 3.5))
	    {
	      gainFact->alfaGF += 0.15;
	      if (gainFact->alfaGF > 0.8) gainFact->alfaGF = 0.8;
	    }
	  else
	    {
	      gainFact->alfaGF -= 0.3;
	      if (gainFact->alfaGF < 0.1) gainFact->alfaGF = 0.1;
	    }
	}

      for (i=0 ; i<WF_MEL_ORDER ; i++)
	W[i] = gainFact->alfaGF * W[i] + (1.0 - gainFact->alfaGF) * 1.0;
    }
}

/*----------------------------------------------------------------------------
 * FUNCTION NAME: DoFilterWindowing
 *
 * PURPOSE:       Filter windowing
 *
 * INPUT:
 *   *filterIRIn  Pointer to filter impulse response
 *   *hanningWindow Pointer to Hanning window
 *
 * OUTPUT:
 *   filterIROut  Output windowed filter impulse response
 *
 * RETURN VALUE:
 *   none
 *
 *---------------------------------------------------------------------------*/
static void DoFilterWindowing (X_FLOAT32 *filterIRIn, X_FLOAT32 *hanningWindow, X_FLOAT32 *filterIROut)
{
  X_INT16 i, j;

  for (i=NS_HALF_FILTER_LENGTH,j=0 ; i<NS_FILTER_LENGTH ; i++,j++)
    {
      filterIROut [NS_HALF_FILTER_LENGTH+j] = filterIRIn [j] * hanningWindow [i];
      filterIROut [NS_HALF_FILTER_LENGTH-j] = filterIROut [NS_HALF_FILTER_LENGTH+j];
    }
}




static void DCOffsetFil (X_FLOAT32 *Data, DC_FILTER *prevSamples, X_INT16 DataLength)
{
  X_INT16 i;
  X_FLOAT32 aux;
  X_FLOAT32 *Prev_x = &(prevSamples->lastSampleIn);
  X_FLOAT32 *Prev_y = &(prevSamples->lastDCOut);

  // y[n] = x[n] - x[n-1] + 0.9990234375 * y[n-1]

  for (i=0 ; i<DataLength ; i++)
    {
      aux = Data[i];
      Data[i] = Data[i] - *Prev_x + 0.9990234375 * *Prev_y;
      *Prev_x = aux;
      *Prev_y = Data[i];
    }
}

void
DoMelIDCT(float *inData, float **melIDCTbasis, int melOrder, int timeLength)
{
  int t, f;
  float *outData;

  outData = (float *) malloc(sizeof(float) * timeLength);
  if (outData == NULL) {
    fprintf(stderr, "ERROR:   Memory allocation error occured!\r\n");
    exit(0);
  }
  for (t = 0; t < timeLength; t++) {
    outData[t] = 0.0;
    for (f = 0; f < melOrder; f++) {
      outData[t] += inData[f] * melIDCTbasis[t][f];
    }
  }
  for (t = 0; t < timeLength; t++)
    inData[t] = outData[t];

  free(outData);
}


int DoNoiseSup(X_INT16 *SigBuf,X_FLOAT32 *OutBuf,NoiseSupStructX *NSX)
{
    X_INT16 i;
    X_INT16 fstage;

	X_FLOAT32 *FirstStageInFloatBuffer  = NSX->nsVar.buffers.FirstStageInFloatBuffer;
    X_FLOAT32 *SecondStageInFloatBuffer = NSX->nsVar.buffers.SecondStageInFloatBuffer;

    X_INT32 *nbFramesInFirstStage     = &(NSX->nsVar.buffers.nbFramesInFirstStage);
    X_INT32 *nbFramesInSecondStage    = &(NSX->nsVar.buffers.nbFramesInSecondStage);
    X_INT32 *nbFramesOutSecondStage   = &(NSX->nsVar.buffers.nbFramesOutSecondStage);

    X_INT16 *indexBuffer = NULL;     //

    X_FLOAT32 *prvFrame = NULL;      //
    X_FLOAT32 *curFrame = NULL;      //
    X_FLOAT32 *nSigSE = NULL;        //
    X_FLOAT32 *PSDMeanBuffer = NULL; //



    X_FLOAT32 *W = NSX->nsTmp.tmpMem + NS_SPEC_ORDER;               // scratch memory
    X_FLOAT32 *filterIR = NSX->nsTmp.tmpMem;                  // scratch memory
    X_FLOAT32 *signalIn = NSX->nsTmp.tmpMem;                  // scratch memory
    X_FLOAT32 *signalOut = NSX->nsTmp.tmpMem + NS_SPEC_ORDER; // scratch memory
    X_FLOAT32 *PSDMeaned = NSX->nsTmp.tmpMem;                 // scratch memory


//    X_INT16 bufData16kSize;          //

//    X_FLOAT32 fb16k [NS_SPEC_ORDER]; //
//    X_FLOAT32 *bufferData16k;        //
//    X_FLOAT32 *CodeForBands16k;      //
//    X_FLOAT32 *tmpWorkSpace16k;      // temporary work space for 16k processing
//    X_FLOAT32 *BandsForCoding16k;    //

    //
//    MelFB_Window *FirstWindow16k;    //

    BOOLEAN dataToProcess;           // data to process ?

	for (i=0 ; i<NS_FRAME_SHIFT ; i++)
    FirstStageInFloatBuffer [NS_DATA_IN_BUFFER + i] = (float)SigBuf[i];


	(*nbFramesInFirstStage)++;


    for (fstage=0 ; fstage < NS_stage ; fstage++)
    {
      dataToProcess = FALSE;

	  if ((fstage == 0) && ((*nbFramesInFirstStage - *nbFramesInSecondStage) > NS_NB_FRAMES_LATENCY))
	  {
	  /*-------------------
	   * process 1st stage
	   *-------------------*/
	  dataToProcess = TRUE;

	  nSigSE = NSX->nsVar.spectrum.nSigSE1;
	  PSDMeanBuffer = &(NSX->nsVar.spectrum.PSDMeanBuffer1 [0][0]);
	  indexBuffer = &(NSX->nsVar.spectrum.indexBuffer1);

	  for (i=0 ; i<NS_FRAME_LENGTH ; i++)
	    signalIn[i] = FirstStageInFloatBuffer [NS_ANALYSIS_WINDOW_8K + i];

	  prvFrame = &(NSX->nsVar.buffers.FirstStageInFloatBuffer [NS_PRV_FRAME]);
	  curFrame = &(NSX->nsVar.buffers.FirstStageInFloatBuffer [NS_CUR_FRAME]);
	  }

	  else if ((fstage == 1) && ((*nbFramesInSecondStage - *nbFramesOutSecondStage) > NS_NB_FRAMES_LATENCY))
	  {
	    /*-------------------
	     * process 2nd stage
	     *-------------------*/
	  dataToProcess = TRUE;

	  nSigSE = NSX->nsVar.spectrum.nSigSE2;
	  PSDMeanBuffer = &(NSX->nsVar.spectrum.PSDMeanBuffer2 [0][0]);
	  indexBuffer = &(NSX->nsVar.spectrum.indexBuffer2);

	  for (i=0 ; i<NS_FRAME_LENGTH ; i++)
	    signalIn[i] = SecondStageInFloatBuffer [NS_ANALYSIS_WINDOW_8K + i];

	  prvFrame = &(NSX->nsVar.buffers.SecondStageInFloatBuffer [NS_PRV_FRAME]);
	  curFrame = &(NSX->nsVar.buffers.SecondStageInFloatBuffer [NS_CUR_FRAME]);
	 }

	 if (!dataToProcess) continue; // no processing required

	 DoSigWindowing (signalIn, NSX->nsTmp.sigWindow, NS_FRAME_LENGTH, NS_FFT_LENGTH);

     rfft (signalIn, NS_FFT_LENGTH, NS_FFT_ORDER);

      /*-------------------------------------------------------------
       * FFT spectrum (signalIn) -> power spectrum (NSX->nSigSE)
       *-------------------------------------------------------------*/
     FFTtoPSD (signalIn, nSigSE, NS_FFT_LENGTH);

     PSDMean (indexBuffer, nSigSE, PSDMeaned, PSDMeanBuffer);

      /*---------------------------
       * VAD for Noise Suppression
       *---------------------------*/
     VAD (fstage, &(NSX->nsVar.vadNS), curFrame);

      /*--------------------------
       * filter gains calculation
       *--------------------------*/
     FilterCalc (fstage, NSX, PSDMeaned, W);

     DoMelFB (W, NSX->nsTmp.FirstWindow);


	 DoGainFact (fstage, NSX, W);

      /*-----------------
       * mel inverse DCT
       *-----------------*/
     DoMelIDCT (W, NSX->nsTmp.melIDCTbasis, WF_MEL_ORDER, WF_MEL_ORDER);
      for (i=1 ; i<WF_MEL_ORDER ; i++) W [2 * WF_MEL_ORDER - 1 - i] = W [i];

      /*------------------
       * filter windowing
       *------------------*/
     DoFilterWindowing (W, NSX->nsTmp.IRWindow, filterIR);

      /*--------------------------------------------------------------
       * apply WF to noisy signal, output samples stored in signalOut
       *--------------------------------------------------------------*/
     ApplyWF (curFrame, prvFrame, filterIR, signalOut, NS_FRAME_SHIFT, NS_HALF_FILTER_LENGTH);

	 if (fstage == 0)
	 {
	  for (i=0 ; i<NS_FRAME_SHIFT ; i++)
      {
          SecondStageInFloatBuffer[NS_DATA_IN_BUFFER + i] = signalOut[i];
          if(NS_stage == 1)
          {
              DCOffsetFil (signalOut, &(NSX->nsVar.prevSamples), NS_FRAME_SHIFT);
              for (i=0 ; i<NS_FRAME_SHIFT ; i++)
              {
                  OutBuf[i] = signalOut[i];

              }
          }

      }


	  NSX->nsVar.buffers.nbFramesInSecondStage++;
	 }

	 else if (fstage == 1)
	 {
       DCOffsetFil (signalOut, &(NSX->nsVar.prevSamples), NS_FRAME_SHIFT);
	   for (i=0 ; i<NS_FRAME_SHIFT ; i++)
	   {
		   OutBuf[i] = signalOut[i];
	   }

	   NSX->nsVar.buffers.nbFramesOutSecondStage++;

	 }

    }

     if (NSX->nsVar.buffers.nbFramesInFirstStage)
	 {
      // copy the content of FirstStageInFloatBuffer [80..319] to [0..239]
      for (i=0 ; i<NS_DATA_IN_BUFFER ; i++)
	     FirstStageInFloatBuffer [i] = FirstStageInFloatBuffer [i + NS_FRAME_SHIFT];
     }

  /*---------------------------------------------
   * shift data in NSX->SecondStageInFloatBuffer
   *---------------------------------------------*/
     if (NSX->nsVar.buffers.nbFramesInSecondStage)
       for (i=0 ; i<NS_DATA_IN_BUFFER ; i++)
         SecondStageInFloatBuffer [i] = SecondStageInFloatBuffer [i + NS_FRAME_SHIFT];


   return TRUE;


}
