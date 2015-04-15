

#include <stdio.h>
#include "type.h"



static int
ReadBufWave(FILE * fp_in, X_INT16 * buf, int nSamples)
{
  X_INT16 s;
  int i;

  for (i = 0; i < nSamples; i++)
    {
    if (fread(&s, sizeof(X_INT16), 1, fp_in) != 1)
      return FALSE;
    buf[i] = (X_INT16) s;
    }
  return TRUE;
}



int main(int argc, char *argv[])
{
    int i;
	FILE *fp_in,*fp_out;
	X_INT16 *SigBuf,*SigOut;
	X_FLOAT32 *OutBuf;
	int bpcount = 0;
	char  temp[44];

	long FrameCounter=0;
	NoiseSupStructX *NSX;
	LOGMMSE_VAR *logMMSE;

	/*NSX=(NoiseSupStructX *)malloc(sizeof(NoiseSupStructX ));
	NoiseReInit(NSX);*/
	logMMSE = (LOGMMSE_VAR *)malloc(sizeof(LOGMMSE_VAR ));
    logMMSE_Init(logMMSE);
	fp_in=fopen("test_wav/sp02_babble_sn01.wav","rb");
	fp_out=fopen("test_wav/out_logMMSE.pcm","wb");

	fread(temp,sizeof(char),22*2,fp_in);
	if(!fp_in)
		printf("Unable to read!");
	SigBuf = (X_INT16 *)malloc(sizeof(X_INT16)*FRAME_LEN);
	/*SigOut = (X_INT16 *)malloc(sizeof(X_INT16)*FRAME_LEN);
	OutBuf = (X_FLOAT32 *)malloc(sizeof(X_FLOAT32)*FRAME_LEN);*/

    while (ReadBufWave(fp_in, SigBuf, FRAME_LEN))
	{
      FrameCounter++;
	  noise_estimate(SigBuf,logMMSE);
	  if(FrameCounter == 7)
        break;
	}
    fclose(fp_in);
	free(SigBuf);

    FrameCounter = 0;
    SigBuf = (X_INT16 *)malloc(sizeof(X_INT16)*FRAME_SHIFT);
	SigOut = (X_INT16 *)malloc(sizeof(X_INT16)*FRAME_SHIFT);
	OutBuf = (X_FLOAT32 *)malloc(sizeof(X_FLOAT32)*FRAME_SHIFT);
	fp_in=fopen("test_wav/sp02_babble_sn01.wav","rb");
	if(!fp_in)
		printf("Unable to read!");
    fread(temp,sizeof(char),22*2,fp_in);

    /*GtAGC *agc = malloc(sizeof(GtAGC));

    SetGtagc(agc,10,-10,-30,1);*/

	while (ReadBufWave(fp_in, SigBuf, FRAME_SHIFT))
	{
      FrameCounter++;
	  //DoNoiseSup(SigBuf,OutBuf,NSX);

      logMMSE_denosie(SigBuf, OutBuf, logMMSE);

	  //DoGtagc(agc,OutBuf,NS_FRAME_SHIFT,1,bpcount,1);
	  //AGC_compute(OutBuf,NS_FRAME_SHIFT,NSX);
	  if(FrameCounter > 1)
      {
         for(i = 0; i < FRAME_SHIFT; i++)
         SigOut[i] = (X_INT16)OutBuf[i];
         fwrite(SigOut,sizeof(X_INT16),FRAME_SHIFT,fp_out);
      }

	}

	printf("Frame number:%d \n",FrameCounter);
	fclose(fp_in);
	fclose(fp_out);
	free(logMMSE);
    free(OutBuf);
	free(SigOut);
	free(SigBuf);


	return TRUE;

}
