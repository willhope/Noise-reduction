/*===============================================================================*/
/*      ETSI ES 202 050   Distributed Speech Recognition                         */
/*      Advanced Front-End Feature Extraction Algorithm & Compression Algorithm  */
/*      C-language software implementation                                       */
/*      Version 1.1.3   October, 2003                                            */
/*===============================================================================*/
/*-------------------------------------------------------------------------------
 *
 * FILE NAME: rfft.c
 * PURPOSE: Real valued, in-place split-radix FFT.
 *
 *-------------------------------------------------------------------------------*/
/*-----------------
 * File Inclusions
 *-----------------*/
#include <math.h>

#include "rfft.h"

/*---------------------------------------------------------------------------
 * FUNCTION NAME: rfft
 *
 * PURPOSE:       Real valued, in-place split-radix FFT
 *
 * INPUT:
 *   x            Pointer to input and output array
 *   n            Length of FFT, must be power of 2
 *
 * OUTPUT         Output order
 *                  Re(0), Re(1), ..., Re(n/2), Im(N/2-1), ..., Im(1)
 *
 * RETURN VALUE
 *   none
 *
 * DESIGN REFERENCE:
 *                IEEE Transactions on Acoustic, Speech, and Signal Processing,
 *                Vol. ASSP-35. No. 6, June 1987, pp. 849-863.
 *
 *                Subroutine adapted from fortran routine pp. 858-859.
 *                Note corrected printing errors on page 859:
 *                    SS1 = SIN(A3) -> should be SS1 = SIN(A);
 *                    CC3 = COS(3)  -> should be CC3 = COS(A3)
 *
 *---------------------------------------------------------------------------*/
void 
rfft(float *x, int n, int m)
{
  int j, i, k, is, id;
  int i0, i1, i2, i3, i4, i5, i6, i7, i8;
  int n2, n4, n8;
  float xt, a0, e, a, a3;
  float t1, t2, t3, t4, t5, t6;
  float cc1, ss1, cc3, ss3;
  float *r0;

  /* Digit reverse counter */

  j = 0;
  r0 = x;

  for (i = 0; i < n - 1; i++) {

    if (i < j) {
      xt = x[j];
      x[j] = *r0;
      *r0 = xt;
    }
    r0++;

    k = n >> 1;

    while (k <= j) {
      j = j - k;
      k >>= 1;
    }
    j += k;
  }

  /* Length two butterflies */
  is = 0;
  id = 4;

  while (is < n - 1) {

    for (i0 = is; i0 < n; i0 += id) {
      i1 = i0 + 1;
      a0 = x[i0];
      x[i0] += x[i1];
      x[i1] = a0 - x[i1];
    }

    is = (id << 1) - 2;
    id <<= 2;
  }

  /* L shaped butterflies */
  n2 = 2;
  for (k = 1; k < m; k++) {
    n2 <<= 1;
    n4 = n2 >> 2;
    n8 = n2 >> 3;
    e = (M_PI * 2) / n2;
    is = 0;
    id = n2 << 1;
    while (is < n) {
      for (i = is; i <= n - 1; i += id) {
	i1 = i;
	i2 = i1 + n4;
	i3 = i2 + n4;
	i4 = i3 + n4;
	t1 = x[i4] + x[i3];
	x[i4] -= x[i3];
	x[i3] = x[i1] - t1;
	x[i1] += t1;

	if (n4 != 1) {
	  i1 += n8;
	  i2 += n8;
	  i3 += n8;
	  i4 += n8;
	  t1 = (x[i3] + x[i4]) / M_SQRT2;
	  t2 = (x[i3] - x[i4]) / M_SQRT2;
	  x[i4] = x[i2] - t1;
	  x[i3] = -x[i2] - t1;
	  x[i2] = x[i1] - t2;
	  x[i1] = x[i1] + t2;
	}
      }
      is = (id << 1) - n2;
      id <<= 2;
    }

    for (j = 1; j < n8; j++) {
      a = j * e;
      a3 = 3 * a;
      cc1 = cos(a);
      ss1 = sin(a);
      cc3 = cos(a3);
      ss3 = sin(a3);

      is = 0;
      id = n2 << 1;

      while (is < n) {
	for (i = is; i <= n - 1; i += id) {
	  i1 = i + j;
	  i2 = i1 + n4;
	  i3 = i2 + n4;
	  i4 = i3 + n4;
	  i5 = i + n4 - j;
	  i6 = i5 + n4;
	  i7 = i6 + n4;
	  i8 = i7 + n4;
	  t1 = x[i3] * cc1 + x[i7] * ss1;
	  t2 = x[i7] * cc1 - x[i3] * ss1;
	  t3 = x[i4] * cc3 + x[i8] * ss3;
	  t4 = x[i8] * cc3 - x[i4] * ss3;
	  t5 = t1 + t3;
	  t6 = t2 + t4;
	  t3 = t1 - t3;
	  t4 = t2 - t4;
	  t2 = x[i6] + t6;
	  x[i3] = t6 - x[i6];
	  x[i8] = t2;
	  t2 = x[i2] - t3;
	  x[i7] = -x[i2] - t3;
	  x[i4] = t2;
	  t1 = x[i1] + t5;
	  x[i6] = x[i1] - t5;
	  x[i1] = t1;
	  t1 = x[i5] + t4;
	  x[i5] = x[i5] - t4;
	  x[i2] = t1;
	}
	is = (id << 1) - n2;
	id <<= 2;
      }
    }
  }
}
