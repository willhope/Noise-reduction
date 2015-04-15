#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "type.h"
#include "fft.h"
#define max(a,b) ((a>b)?(a):(b))
#define min(a,b) ((a<b)?(a):(b))


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

double expp(double x)
{
    int m,i,j;
    double s,p,ep,h,aa,bb,w,xx,g,r,q;
    static double t[5]={-0.9061798459,-0.5384693101,0.0,
                         0.5384693101,0.9061798459};
    static double c[5]={0.2369268851,0.4786286705,0.5688888889,
                        0.4786286705,0.2369268851};
    m=1;
    if (x==0) x=1.0e-10;
    if (x<0.0) x=-x;
    r=0.57721566490153286060651;
    q=r+log(x);
    h=x; s=fabs(0.0001*h);
    p=1.0e+35; ep=0.000001; g=0.0;
    while ((ep>=0.0000001)&&(fabs(h)>s))
      { g=0.0;
        for (i=1;i<=m;i++)
          { aa=(i-1.0)*h; bb=i*h;
            w=0.0;
            for (j=0;j<=4;j++)
              { xx=((bb-aa)*t[j]+(bb+aa))/2.0;
                w=w+(exp(-xx)-1.0)/xx*c[j];
              }
            g=g+w;
          }
        g=g*h/2.0;
        ep=fabs(g-p)/(1.0+fabs(g));
        p=g; m=m+1; h=x/m;
      }
    g=q+g;
    return(g);
}

int logMMSE_denosie(X_INT16 *sigBuf, X_FLOAT32 *OutBuf, LOGMMSE_VAR *logMMSE)
{
    X_INT16 i = 0;
    X_INT16 fstage;
    X_FLOAT32 tmp_var       = 0.0;
    X_FLOAT32 vad_decision  = 0.0;
    X_FLOAT32 log_sigma_sum = 0.0;


    X_FLOAT32 *Buf2denoise  = &(logMMSE->Buf2denoise);
    X_FLOAT32 *old_denosebuf= &(logMMSE->old_denosebuf);
    Complex   *FFTIn        = &(logMMSE->FFTIn);
    X_FLOAT32 *spect        = &(logMMSE->spect);
    X_FLOAT32 *Xk_prev      = &(logMMSE->Xk_prev);
    X_FLOAT32 *noiseMean    = &(logMMSE->noiseMean);
    X_FLOAT32 *sig          = &(logMMSE->sig);
    X_FLOAT32 *gammak       = &(logMMSE->gammak);
    X_FLOAT32 *ksi          = &(logMMSE->ksi);
    X_FLOAT32 *log_sigma_k  = &(logMMSE->log_sigma_k);
    X_FLOAT32 *A            = &(logMMSE->A);
    X_FLOAT32 *vk           = &(logMMSE->vk);
    X_FLOAT32 *ei_vk        = &(logMMSE->ei_vk);
    X_FLOAT32 *hw           = &(logMMSE->hw);
    X_INT32 *Frame_count   = &(logMMSE->nbFrame);

    X_FLOAT32 aa  = 0.98;
    X_FLOAT32 mu  = 0.98;
    X_FLOAT32 eta = 0.15;
    X_FLOAT32 ksi_min = pow(10,-2.5);

    (*Frame_count) ++;

    for(fstage = 0; fstage < 1; fstage ++)
    {
      if((*Frame_count) == 1)
        {
            for(i = 0; i < FRAME_SHIFT; i ++)
            {
                Buf2denoise[i] = sigBuf[i];
            }
            continue;
        }

        for(i = 0; i < FRAME_SHIFT; i ++)
        {
            Buf2denoise[i + FRAME_SHIFT] = sigBuf[i];
        }

        DoSigWindowing (Buf2denoise, logMMSE->sigWindow, FRAME_LEN, FFT_LEN);
        for(i = 0; i < FFT_LEN; i ++)
        {
            FFTIn[i].real = Buf2denoise[i];
            FFTIn[i].imag = 0.0;
        }



        //fft4(FFTIn, 4);
        fft(FFTIn,FFT_LEN);

        for(i = 0; i < FFT_LEN; i ++)
        {
           spect[i] = FFTIn[i].real * FFTIn[i].real + FFTIn[i].imag * FFTIn[i].imag;
        }



        for(i = 0; i < FFT_LEN; i ++)
        {

            sig[i] = spect[i];
            gammak[i] = min((sig[i]/noiseMean[i]),40);
            //printf("%f ",gammak[i]);
            if((*Frame_count == 2))
            {
                ksi[i] = aa + (1 - aa) * max(gammak[i]-1,0);
            }
            else
            {
                ksi[i] = aa * Xk_prev[i] / noiseMean[i] + (1-aa) * max(gammak[i]-1,0);
                ksi[i] = max(ksi_min,ksi[i]);
                //printf("%f ",Xk_prev[i]);
            }

            log_sigma_k[i] = gammak[i] * ksi[i] / (1 + ksi[i]) - log(1 + ksi[i]);

            log_sigma_sum += log_sigma_k[i];
        }

        for(i = 0; i < FFT_LEN; i ++)
        {
            //printf("%f ",log_sigma_k[i]);
        }
       /* printf("\n");
        printf("\n");*/
        vad_decision = log_sigma_sum / FRAME_LEN;
        if(vad_decision < eta)
        {
            for(i = 0; i < FFT_LEN; i ++)
                noiseMean[i] = mu * noiseMean[i] + (1 - mu) * sig[i];
        }
        for(i = 0; i < FFT_LEN; i ++)
        {
            A[i]       = ksi[i] / (1 + ksi[i]);
            vk[i]      = A[i] * gammak[i];
            //printf("%f ",vk[i]);
            ei_vk[i]   = 0.5 * (-expp(vk[i]));

            hw[i]      = A[i] * exp(ei_vk[i]);
            tmp_var    = sqrt(spect[i]) * hw[i];
            Xk_prev[i] = tmp_var * tmp_var;
            FFTIn[i].real = FFTIn[i].real * hw[i];
            FFTIn[i].imag = FFTIn[i].imag * hw[i];
        }

        //ifft4(FFTIn,4);
        ifft(FFTIn,FFT_LEN);
        for(i = 0; i < FRAME_SHIFT; i ++)
        {
            OutBuf[i] = old_denosebuf[i] + FFTIn[i].real;
            //printf("%f ",OutBuf[i]);
            old_denosebuf[i] = FFTIn[FRAME_SHIFT + i].real;
        }


        for(i = 0; i < FRAME_SHIFT; i ++)
        {
            Buf2denoise[i] = sigBuf[i];
        }
    }





}



int noise_estimate(X_INT16 *sigBuf, LOGMMSE_VAR *logMMSE)
{
    X_INT16 i = 0;

    X_FLOAT32 tmp_var       = 0.0;
    X_FLOAT32 vad_decision  = 0.0;


    X_FLOAT32 *Buf2denoise   = &(logMMSE->Buf2denoise);
    Complex   *FFTIn        = &(logMMSE->FFTIn);
    X_FLOAT32 *spect        = &(logMMSE->spect);
    X_FLOAT32 *noiseMean    = &(logMMSE->noiseMean);
    X_INT32 *Frame_count   = &(logMMSE->nbFrame);

    (*Frame_count) ++;


    for(i = 0; i < FRAME_LEN; i ++)
    {
        Buf2denoise[i] = sigBuf[i];
    }
    DoSigWindowing (Buf2denoise, logMMSE->sigWindow, FRAME_LEN, FFT_LEN);
    for(i = 0; i < FFT_LEN; i ++)
    {
        FFTIn[i].real = Buf2denoise[i];
        FFTIn[i].imag = 0.0;
    }

    fft4(FFTIn, 4);
    //fft(FFTIn,128)

    for(i = 0; i < FFT_LEN; i ++)
    {
        spect[i] = FFTIn[i].real * FFTIn[i].real + FFTIn[i].imag * FFTIn[i].imag;
    }

    if((*Frame_count) < 7)
    {
        for(i = 0; i < FFT_LEN; i ++)
        {
            noiseMean[i] += sqrt(spect[i]);
        }

       //(*Frame_count) = 0;

    }

    else
    {
        for(i = 0; i < FFT_LEN; i ++)
        {
            noiseMean[i] /= 6;
            noiseMean[i]  = noiseMean[i] * noiseMean[i];
            //printf("%f ",noiseMean[i]);
        }
       /* printf("\n");
        printf("\n");*/
        (*Frame_count) = 0;
    }


    return 0;



}
