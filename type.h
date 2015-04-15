#ifndef _TYPE_H
#define _TYPE_H

#include "fft.h"


#define NS_stage                      (X_INT16)(2)
#define NS_FILTER_LENGTH              (X_INT16)(17)
#define NS_HALF_FILTER_LENGTH         (X_INT16)(8)
#define NS_NB_FRAME_THRESHOLD_NSE     (X_INT16)(100)

#define NS_EPS                        (X_FLOAT32)(exp ((double) -10.0))
#define NS_BETA                       (X_FLOAT32)(0.98)
#define NS_RSB_MIN                    (X_FLOAT32)(0.079432823)
#define NS_LAMBDA_NSE                 (X_FLOAT32)(0.99)
#define NS_SPEC_FLOOR                 (X_FLOAT32)(exp ((double) -10.0))
#define NS_LOG_SPEC_FLOOR             (X_FLOAT32)(-10.0)

/*----------------------------------------------------------------------------
 *                            logMMSE
 *----------------------------------------------------------------------------*/
#define samplRate                     (X_INT16)(8000)
#define FRAME_LEN                     (X_INT16)(20*samplRate/1000)
#define OVERLAP_RATE                  (X_FLOAT32) (0.5)
#define FRAME_SHIFT                   (X_INT16)(FRAME_LEN * OVERLAP_RATE)
#define FFT_LEN                       (X_INT16)(256)
#define FFT_ORDER                     (X_INT16)(8)

/*----------------------------------------------------------------------------
 *                                     FFT
 *----------------------------------------------------------------------------*/
#define NS_SPEC_ORDER                 (X_INT16)(65)
#define NS_FFT_ORDER                  (X_INT16)(8)
#define NS_FFT_LENGTH                 (X_INT16)(256)
#define NS_FREQUENCY_BINS             (X_INT16)(129)

/*----------------------------------------------------------------------------
 *                                  BUFFERING
 *----------------------------------------------------------------------------*/
#define NS_PRV_FRAME                  (X_INT16)(0)
#define NS_CUR_FRAME                  (X_INT16)(80)
#define NS_FRAME_SHIFT                (X_INT16)(80)
#define NS_BUFFER_SIZE                (X_INT16)(320)
#define NS_FRAME_LENGTH               (X_INT16)(200)
#define NS_DATA_IN_BUFFER             (X_INT16)(240)
#define NS_SCRATCH_MEM_SIZE           (NS_FFT_LENGTH)
#define NS_NB_FRAMES_LATENCY          (X_INT16)(2)
#define NS_ANALYSIS_WINDOW_8K         (X_INT16)(60)
#define NS_ANALYSIS_WINDOW_16K        (X_INT16)(80)
#define NS_NB_FRAMES_IN_BUFFER        (X_INT16)(4)

/*----------------------------------------------------------------------------
 *                                  PSD MEAN
 *----------------------------------------------------------------------------*/
#define NS_PSD_MEAN_ORDER             (X_INT16)(2)

/*----------------------------------------------------------------------------
 *                                     VAD
 *----------------------------------------------------------------------------*/
#define NS_HANGOVER                   (X_INT16)(15)
#define NS_MIN_FRAME                  (X_INT16)(10)
#define NS_SNR_THRESHOLD_VAD          (X_INT16)(15)
#define NS_SNR_THRESHOLD_UPD_LTE      (X_INT16)(20)
#define NS_NB_FRAME_THRESHOLD_LTE     (X_INT16)(10)
#define NS_MIN_SPEECH_FRAME_HANGOVER  (X_INT16)(4)

#define NS_ENERGY_FLOOR               (X_FLOAT32)(80.0)
#define NS_LAMBDA_LTE_LOWER_E         (X_FLOAT32)(0.97)
#define NS_LAMBDA_LTE_HIGHER_E        (X_FLOAT32)(0.99)

#define WF_MEL_ORDER                  25
#define PIx2                          6.28318530717958647692


#define FALSE 0
#define TRUE (!FALSE)


#define BOOLEAN int
#define X_FLOAT32 float
#define X_INT16 short int
#define X_INT32 int

typedef struct vad_data_ns VAD_DATA_NS;
typedef struct vad_data_fd VAD_DATA_FD;
typedef struct vad_data_ca VAD_DATA_CA;
typedef struct ns_var      NS_VAR;
typedef struct ns_tmp      NS_TMP;
typedef struct dc_filter   DC_FILTER;
typedef struct gain_fact   GAIN_FACT;
typedef struct agc_data_ns AGC_DATA_NS;
typedef struct buffers     BUFFERS;
typedef struct spectrum    SPECTRUM;
typedef struct MelFB_Window MelFB_Window;
typedef struct NoiseSupStructX NoiseSupStructX;
typedef struct LOGMMSE_VAR LOGMMSE_VAR;

struct dc_filter
{
  X_FLOAT32 lastSampleIn;                  // last input sample of DC offset compensation
  X_FLOAT32 lastDCOut;                     // last output sample of DC offset compensation
} ;

struct vad_data_ns
{
  X_INT32   nbFrame [2];                   // frame number
  X_INT16   flagVAD;                       // VAD flag (1 == SPEECH, 0 == NON SPEECH)
  X_INT16   hangOver;                      // hangover
  X_INT16   nbSpeechFrames;                // nb speech frames (used to set hangover)
  X_FLOAT32 meanEn;                        // mean energy
};

struct vad_data_ca
{
  X_INT16   flagVAD;                       // VAD flag (1 == SPEECH, 0 == NON SPEECH)
  X_INT16   hangOver;                      // hangover
  X_INT16   nbSpeechFrames;                // nb speech frames (used to set hangover)
  X_FLOAT32 meanCA;                        // mean energy
};

struct vad_data_fd
{
  X_FLOAT32 MelMean;			   //
  X_FLOAT32 VarMean;			   //
  X_FLOAT32 AccTest;			   //
  X_FLOAT32 AccTest2;			   //
  X_FLOAT32 SpecMean;			   //
  X_FLOAT32 MelValues[2];		   //
  X_FLOAT32 SpecValues;		           //
  X_FLOAT32 SpeechInVADQ;		   //
  X_FLOAT32 SpeechInVADQ2;		   //
};

struct gain_fact
{
  X_FLOAT32 denEn1 [3];                    // previous denoised frames energies
  X_FLOAT32 lowSNRtrack;                   // low SNR track
  X_FLOAT32 alfaGF;                        // gain factor applied in 2nd stage
};

struct agc_data_ns
{
  X_FLOAT32 Last_Spower;
  X_FLOAT32 Spower;
  X_FLOAT32 agc_gain;
  X_INT16   nb_adapt;

};

struct buffers
{
  X_INT32 nbFramesInFirstStage;            // nb frames in first stage
  X_INT32 nbFramesInSecondStage;           // nb frames in second stage
  X_INT32 nbFramesOutSecondStage;          // nb frames out of second stage

  X_FLOAT32 FirstStageInFloatBuffer  [NS_BUFFER_SIZE]; // first stage buffer
  X_FLOAT32 SecondStageInFloatBuffer [NS_BUFFER_SIZE]; // second stage buffer
};

struct spectrum
{
  X_FLOAT32 nSigSE1 [NS_SPEC_ORDER];       // 1st stage noisy signal spectrum estimation
  X_FLOAT32 nSigSE2 [NS_SPEC_ORDER];       // 2nd stage noisy signal spectrum estimation
  X_FLOAT32 noiseSE1 [NS_SPEC_ORDER];      // 1st stage noise spectrum estimation
  X_FLOAT32 noiseSE2 [NS_SPEC_ORDER];      // 2nd stage noise spectrum estimation
  X_FLOAT32 denSigSE1 [NS_SPEC_ORDER];     // 1st stage denoised signal spectrum estimation
  X_FLOAT32 denSigSE2 [NS_SPEC_ORDER];     // 2nd stage denoised signal spectrum estimation

  X_INT16   indexBuffer1;                  // where to enter new PSD, alternatively 0 and 1
  X_INT16   indexBuffer2;                  // where to enter new PSD, alternatively 0 and 1
  X_FLOAT32 PSDMeanBuffer1 [NS_SPEC_ORDER][NS_PSD_MEAN_ORDER]; // 1st stage PSD Mean buffer
  X_FLOAT32 PSDMeanBuffer2 [NS_SPEC_ORDER][NS_PSD_MEAN_ORDER]; // 2nd stage PSD Mean buffer
};

struct MelFB_Window
{
  int StartingPoint;
  int Length;
  float *Data;
  struct MelFB_Window *Next;
};

struct ns_var
{
  X_INT16     SampFreq;                    // sampling frequency
  BOOLEAN     Do16kHzProc;                 // do 16 kHz processing ?
  BUFFERS     buffers;                     // signal buffers
  DC_FILTER   prevSamples;                 // previous samples for DC offset removal
  GAIN_FACT   gainFact;                    // gain factorization variables
  AGC_DATA_NS agc_data_ns;
  VAD_DATA_NS vadNS;                       // VAD for noise suppression data
  VAD_DATA_FD vadFD;                       // VAD for frame dropping data
  VAD_DATA_CA vadCA;                       // VAD for frame dropping data
  SPECTRUM    spectrum;                    // spectrum data
};


struct ns_tmp
{
  MelFB_Window *FirstWindow;               // chained list for Mel filtering
  X_FLOAT32 **melIDCTbasis;                // mel-frequency inverse DCT basis
  X_FLOAT32 IRWindow [NS_FILTER_LENGTH];   // filter impulse response window
  X_FLOAT32 sigWindow [NS_FRAME_LENGTH];   // signal window
  X_FLOAT32 tmpMem [NS_SCRATCH_MEM_SIZE];  // scratch memory
};



struct NoiseSupStructX
{
  NS_VAR    nsVar;                         // non sharable data
  NS_TMP    nsTmp;                         // sharable data
};

struct LOGMMSE_VAR
{
  X_INT16     SampFreq;                    // sampling frequency
  X_INT32     nbFrame;
  X_FLOAT32   Buf2denoise[FFT_LEN];
  X_FLOAT32   old_denosebuf[FRAME_SHIFT];
  Complex     FFTIn[FFT_LEN];
  X_FLOAT32   spect[FFT_LEN];
  X_FLOAT32   noiseMean[FFT_LEN];
  X_FLOAT32   sig[FFT_LEN];
  X_FLOAT32   gammak[FFT_LEN];
  X_FLOAT32   ksi[FFT_LEN];
  X_FLOAT32   Xk_prev[FFT_LEN];
  X_FLOAT32   log_sigma_k[FFT_LEN];
  X_FLOAT32   A[FFT_LEN];
  X_FLOAT32   vk[FFT_LEN];
  X_FLOAT32   ei_vk[FFT_LEN];
  X_FLOAT32   hw[FFT_LEN];
  X_FLOAT32   sigWindow [FRAME_LEN];
  AGC_DATA_NS agc_data_ns;
  VAD_DATA_NS vadNS;                       // VAD for noise suppression data
};

int DoNoiseSup(X_INT16 *SigBuf,X_FLOAT32 *OutData,NoiseSupStructX *NSX);

void NoiseReInit(NoiseSupStructX *NSX);

void logMMSE_Init(LOGMMSE_VAR *NSX);

#endif
