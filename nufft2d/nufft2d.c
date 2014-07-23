/* --------------------------------------------------- */
/* Automatically generated by mwrap                    */
/* --------------------------------------------------- */

/* Code generated by mwrap */
/*
  Copyright statement for mwrap:

  mwrap -- MEX file generation for MATLAB and Octave
  Copyright (c) 2007-2008 David Bindel

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in
  all copies or substantial portions of the Software.

  You may distribute a work that contains part or all of the source code
  generated by mwrap under the terms of your choice.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
  THE SOFTWARE.
*/

#include <mex.h>
#include <stdio.h>
#include <string.h>


#ifndef ulong
#  define ulong unsigned long
#endif
#ifndef uint
#  define uint  unsigned int
#endif
#ifndef uchar
#  define uchar unsigned char
#endif


/*
 * Support for 32-bit and 64-bit MEX files
 */
#ifndef mwSize
#  define mwSize int
#endif
#ifndef mwIndex
#  define mwIndex int
#endif
#ifndef mwSignedIndex
#  define mwSignedIndex int
#endif


/*
 * Records for call profile.
 */
int* mexprofrecord_= NULL;


/*
 * Support routines for copying data into and out of the MEX stubs
 */

void* mxWrapGetP(const mxArray* a, const char* fmt, const char** e)
{
    void* p = 0;
    mxArray* ap;
    if (mxGetClassID(a) == mxDOUBLE_CLASS && 
        mxGetM(a)*mxGetN(a) == 1 && *mxGetPr(a) == 0)
        return p;
    if (mxIsChar(a)) {
        char pbuf[128];
        mxGetString(a, pbuf, sizeof(pbuf));
        sscanf(pbuf, fmt, &p);
    } 
#ifdef R2008OO
    else if (ap = mxGetProperty(a, 0, "mwptr")) {
        return mxWrapGetP(ap, fmt, e);
    }
#endif
    if (p == 0)
        *e = "Invalid pointer";
    return p;
}

mxArray* mxWrapCreateP(void* p, const char* fmt)
{
    if (p == 0) {
        mxArray* z = mxCreateDoubleMatrix(1,1, mxREAL);
        *mxGetPr(z) = 0;
        return z;
    } else {
        char pbuf[128];
        sprintf(pbuf, fmt, p);
        return mxCreateString(pbuf);
    }
}

mxArray* mxWrapStrncpy(const char* s)
{
    if (s) {
        return mxCreateString(s);
    } else {
        mxArray* z = mxCreateDoubleMatrix(1,1, mxREAL);
        *mxGetPr(z) = 0;
        return z;
    }
}

double mxWrapGetScalar(const mxArray* a, const char** e)
{
    if (!a || mxGetClassID(a) != mxDOUBLE_CLASS || mxGetM(a)*mxGetN(a) != 1) {
        *e = "Invalid scalar argument";
        return 0;
    }
    return *mxGetPr(a);
}

char* mxWrapGetString(const mxArray* a, const char** e)
{
    char* s;
    mwSize slen;
    if (!a || (!mxIsChar(a) && mxGetM(a)*mxGetN(a) > 0)) {
        *e = "Invalid string argument";
        return NULL;
    }
    slen = mxGetM(a)*mxGetN(a) + 1;
    s = (char*) mxMalloc(slen);
    if (mxGetM(a)*mxGetN(a) == 0)
        *s = 0;
    else
        mxGetString(a, s, slen);
    return s;
}


#define mxWrapGetArrayDef(func, T) \
T* func(const mxArray* a, const char** e)     \
{ \
    T* array; \
    mwSize arraylen; \
    mwIndex i; \
    T* p; \
    double* q; \
    if (!a || mxGetClassID(a) != mxDOUBLE_CLASS) { \
        *e = "Invalid array argument"; \
        return 0; \
    } \
    arraylen = mxGetM(a)*mxGetN(a); \
    array = (T*) mxMalloc(mxGetM(a)*mxGetN(a) * sizeof(T)); \
    p = array; \
    q = mxGetPr(a); \
    for (i = 0; i < arraylen; ++i) \
        *p++ = (T) (*q++); \
    return array; \
}


#define mxWrapCopyDef(func, T) \
void func(mxArray* a, const T* q, mwSize n) \
{ \
    mwIndex i; \
    double* p = mxGetPr(a); \
    for (i = 0; i < n; ++i) \
        *p++ = *q++; \
}


#define mxWrapReturnDef(func, T) \
mxArray* func(const T* q, mwSize m, mwSize n) \
{ \
    mwIndex i; \
    double* p; \
    if (!q) { \
        return mxCreateDoubleMatrix(0,0, mxREAL); \
    } else { \
        mxArray* a = mxCreateDoubleMatrix(m,n, mxREAL); \
        p = mxGetPr(a); \
        for (i = 0; i < m*n; ++i) \
            *p++ = *q++; \
        return a; \
    } \
}


#define mxWrapGetScalarZDef(func, T, ZT, setz) \
void func(T* z, const mxArray* a) \
{ \
    double* pr = mxGetPr(a); \
    double* pi = mxGetPi(a); \
    setz(z, (ZT) *pr, (pi ? (ZT) *pi : (ZT) 0)); \
}


#define mxWrapGetArrayZDef(func, T, ZT, setz) \
T* func(const mxArray* a, const char** e) \
{ \
    T* array; \
    mwSize arraylen; \
    mwIndex i; \
    T* p; \
    double* qr; \
    double* qi; \
    if (!a || mxGetClassID(a) != mxDOUBLE_CLASS) { \
        *e = "Invalid array argument"; \
        return 0; \
    } \
    arraylen = mxGetM(a)*mxGetN(a); \
    array = (T*) mxMalloc(mxGetM(a)*mxGetN(a) * sizeof(T)); \
    p = array; \
    qr = mxGetPr(a); \
    qi = mxGetPi(a); \
    for (i = 0; i < arraylen; ++i) { \
        ZT val_qr = *qr++; \
        ZT val_qi = (qi ? (ZT) *qi++ : (ZT) 0); \
        setz(p, val_qr, val_qi); \
        ++p; \
    } \
    return array; \
}


#define mxWrapCopyZDef(func, T, real, imag) \
void func(mxArray* a, const T* q, mwSize n) \
{ \
    mwIndex i; \
    double* pr = mxGetPr(a); \
    double* pi = mxGetPi(a); \
    for (i = 0; i < n; ++i) { \
        *pr++ = real(*q); \
        *pi++ = imag(*q); \
        ++q; \
    } \
}


#define mxWrapReturnZDef(func, T, real, imag) \
mxArray* func(const T* q, mwSize m, mwSize n) \
{ \
    mwIndex i; \
    double* pr; \
    double* pi; \
    if (!q) { \
        return mxCreateDoubleMatrix(0,0, mxCOMPLEX); \
    } else { \
        mxArray* a = mxCreateDoubleMatrix(m,n, mxCOMPLEX); \
        pr = mxGetPr(a); \
        pi = mxGetPi(a); \
        for (i = 0; i < m*n; ++i) { \
            *pr++ = real(*q); \
            *pi++ = imag(*q); \
            ++q; \
        } \
        return a; \
    } \
}

#include <complex.h>

typedef _Complex double dcomplex;
#define real_dcomplex(z) creal(z)
#define imag_dcomplex(z) cimag(z)
#define setz_dcomplex(z,r,i)  *z = r + i*_Complex_I

typedef _Complex float fcomplex;
#define real_fcomplex(z) crealf(z)
#define imag_fcomplex(z) cimagf(z)
#define setz_fcomplex(z,r,i)  *z = r + i*_Complex_I

/* Array copier definitions */
mxWrapGetArrayDef(mxWrapGetArray_bool, bool)
mxWrapCopyDef    (mxWrapCopy_bool,     bool)
mxWrapReturnDef  (mxWrapReturn_bool,   bool)
mxWrapGetArrayDef(mxWrapGetArray_char, char)
mxWrapCopyDef    (mxWrapCopy_char,     char)
mxWrapReturnDef  (mxWrapReturn_char,   char)
mxWrapGetArrayDef(mxWrapGetArray_double, double)
mxWrapCopyDef    (mxWrapCopy_double,     double)
mxWrapReturnDef  (mxWrapReturn_double,   double)
mxWrapGetArrayDef(mxWrapGetArray_float, float)
mxWrapCopyDef    (mxWrapCopy_float,     float)
mxWrapReturnDef  (mxWrapReturn_float,   float)
mxWrapGetArrayDef(mxWrapGetArray_int, int)
mxWrapCopyDef    (mxWrapCopy_int,     int)
mxWrapReturnDef  (mxWrapReturn_int,   int)
mxWrapGetArrayDef(mxWrapGetArray_long, long)
mxWrapCopyDef    (mxWrapCopy_long,     long)
mxWrapReturnDef  (mxWrapReturn_long,   long)
mxWrapGetArrayDef(mxWrapGetArray_mwIndex, mwIndex)
mxWrapCopyDef    (mxWrapCopy_mwIndex,     mwIndex)
mxWrapReturnDef  (mxWrapReturn_mwIndex,   mwIndex)
mxWrapGetArrayDef(mxWrapGetArray_mwSignedIndex, mwSignedIndex)
mxWrapCopyDef    (mxWrapCopy_mwSignedIndex,     mwSignedIndex)
mxWrapReturnDef  (mxWrapReturn_mwSignedIndex,   mwSignedIndex)
mxWrapGetArrayDef(mxWrapGetArray_mwSize, mwSize)
mxWrapCopyDef    (mxWrapCopy_mwSize,     mwSize)
mxWrapReturnDef  (mxWrapReturn_mwSize,   mwSize)
mxWrapGetArrayDef(mxWrapGetArray_size_t, size_t)
mxWrapCopyDef    (mxWrapCopy_size_t,     size_t)
mxWrapReturnDef  (mxWrapReturn_size_t,   size_t)
mxWrapGetArrayDef(mxWrapGetArray_uchar, uchar)
mxWrapCopyDef    (mxWrapCopy_uchar,     uchar)
mxWrapReturnDef  (mxWrapReturn_uchar,   uchar)
mxWrapGetArrayDef(mxWrapGetArray_uint, uint)
mxWrapCopyDef    (mxWrapCopy_uint,     uint)
mxWrapReturnDef  (mxWrapReturn_uint,   uint)
mxWrapGetArrayDef(mxWrapGetArray_ulong, ulong)
mxWrapCopyDef    (mxWrapCopy_ulong,     ulong)
mxWrapReturnDef  (mxWrapReturn_ulong,   ulong)
mxWrapGetScalarZDef(mxWrapGetScalar_fcomplex, fcomplex,
                    float, setz_fcomplex)
mxWrapGetArrayZDef (mxWrapGetArray_fcomplex, fcomplex,
                    float, setz_fcomplex)
mxWrapCopyZDef     (mxWrapCopy_fcomplex, fcomplex,
                    real_fcomplex, imag_fcomplex)
mxWrapReturnZDef   (mxWrapReturn_fcomplex, fcomplex,
                    real_fcomplex, imag_fcomplex)
mxWrapGetScalarZDef(mxWrapGetScalar_dcomplex, dcomplex,
                    double, setz_dcomplex)
mxWrapGetArrayZDef (mxWrapGetArray_dcomplex, dcomplex,
                    double, setz_dcomplex)
mxWrapCopyZDef     (mxWrapCopy_dcomplex, dcomplex,
                    real_dcomplex, imag_dcomplex)
mxWrapReturnZDef   (mxWrapReturn_dcomplex, dcomplex,
                    real_dcomplex, imag_dcomplex)

#if defined(MWF77_CAPS)
#define MWF77_dirft2d1 DIRFT2D1
#define MWF77_dirft2d2 DIRFT2D2
#define MWF77_dirft2d3 DIRFT2D3
#define MWF77_nufft2d1f90 NUFFT2D1F90
#define MWF77_nufft2d2f90 NUFFT2D2F90
#define MWF77_nufft2d3f90 NUFFT2D3F90
#elif defined(MWF77_UNDERSCORE1)
#define MWF77_dirft2d1 dirft2d1_
#define MWF77_dirft2d2 dirft2d2_
#define MWF77_dirft2d3 dirft2d3_
#define MWF77_nufft2d1f90 nufft2d1f90_
#define MWF77_nufft2d2f90 nufft2d2f90_
#define MWF77_nufft2d3f90 nufft2d3f90_
#else /* f2c convention */
#define MWF77_dirft2d1 dirft2d1_
#define MWF77_dirft2d2 dirft2d2_
#define MWF77_dirft2d3 dirft2d3_
#define MWF77_nufft2d1f90 nufft2d1f90_
#define MWF77_nufft2d2f90 nufft2d2f90_
#define MWF77_nufft2d3f90 nufft2d3f90_
#endif

#ifdef __cplusplus
extern "C" { /* Prevent C++ name mangling */
#endif

#ifndef MWF77_RETURN
#define MWF77_RETURN int
#endif

MWF77_RETURN MWF77_dirft2d1(int*, double*, double*, dcomplex*, int*, int*, int*, dcomplex*);
MWF77_RETURN MWF77_dirft2d2(int*, double*, double*, dcomplex*, int*, int*, int*, dcomplex*);
MWF77_RETURN MWF77_dirft2d3(int*, double*, double*, dcomplex*, int*, int*, double*, double*, dcomplex*);
MWF77_RETURN MWF77_nufft2d1f90(int*, double*, double*, dcomplex*, int*, double*, int*, int*, dcomplex*, int*);
MWF77_RETURN MWF77_nufft2d2f90(int*, double*, double*, dcomplex*, int*, double*, int*, int*, dcomplex*, int*);
MWF77_RETURN MWF77_nufft2d3f90(int*, double*, double*, dcomplex*, int*, double*, int*, double*, double*, dcomplex*, int*);

#ifdef __cplusplus
} /* end extern C */
#endif

/* ---- nufft2d.mw: 52 ----
 * dirft2d1(int[1] nj, double[] xj, double[] yj, dcomplex[] cj, int[1] iflag, int[1] ms, int[1] mt, inout dcomplex[] fk);
 */
const char* stubids1_ = "dirft2d1(i int[x], i double[], i double[], i dcomplex[], i int[x], i int[x], i int[x], io dcomplex[])";

void mexStub1(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    int*        in0_ =0; /* nj         */
    double*     in1_ =0; /* xj         */
    double*     in2_ =0; /* yj         */
    dcomplex*   in3_ =0; /* cj         */
    int*        in4_ =0; /* iflag      */
    int*        in5_ =0; /* ms         */
    int*        in6_ =0; /* mt         */
    dcomplex*   in7_ =0; /* fk         */
    mwSize      dim8_;   /* 1          */
    mwSize      dim9_;   /* 1          */
    mwSize      dim10_;   /* 1          */
    mwSize      dim11_;   /* 1          */

    dim8_ = (mwSize) mxWrapGetScalar(prhs[8], &mw_err_txt_);
    dim9_ = (mwSize) mxWrapGetScalar(prhs[9], &mw_err_txt_);
    dim10_ = (mwSize) mxWrapGetScalar(prhs[10], &mw_err_txt_);
    dim11_ = (mwSize) mxWrapGetScalar(prhs[11], &mw_err_txt_);

    if (mxGetM(prhs[0])*mxGetN(prhs[0]) != dim8_) {
        mw_err_txt_ = "Bad argument size: nj";        goto mw_err_label;
    }

    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != dim9_) {
        mw_err_txt_ = "Bad argument size: iflag";        goto mw_err_label;
    }

    if (mxGetM(prhs[5])*mxGetN(prhs[5]) != dim10_) {
        mw_err_txt_ = "Bad argument size: ms";        goto mw_err_label;
    }

    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != dim11_) {
        mw_err_txt_ = "Bad argument size: mt";        goto mw_err_label;
    }

    if (mxGetM(prhs[0])*mxGetN(prhs[0]) != 0) {
        in0_ = mxWrapGetArray_int(prhs[0], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in0_ = NULL;
    if (mxGetM(prhs[1])*mxGetN(prhs[1]) != 0) {
        in1_ = mxGetPr(prhs[1]);
    } else
        in1_ = NULL;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxGetPr(prhs[2]);
    } else
        in2_ = NULL;
    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 0) {
        in3_ = mxWrapGetArray_dcomplex(prhs[3], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in3_ = NULL;
    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != 0) {
        in4_ = mxWrapGetArray_int(prhs[4], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in4_ = NULL;
    if (mxGetM(prhs[5])*mxGetN(prhs[5]) != 0) {
        in5_ = mxWrapGetArray_int(prhs[5], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in5_ = NULL;
    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != 0) {
        in6_ = mxWrapGetArray_int(prhs[6], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in6_ = NULL;
    if (mxGetM(prhs[7])*mxGetN(prhs[7]) != 0) {
        in7_ = mxWrapGetArray_dcomplex(prhs[7], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in7_ = NULL;
    if (mexprofrecord_)
        mexprofrecord_[1]++;
    MWF77_dirft2d1(in0_, in1_, in2_, in3_, in4_, in5_, in6_, in7_);
    plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[7]), mxGetN(prhs[7]), mxCOMPLEX);
    mxWrapCopy_dcomplex(plhs[0], in7_, mxGetM(prhs[7])*mxGetN(prhs[7]));

mw_err_label:
    if (in0_)  mxFree(in0_);
    if (in3_)  mxFree(in3_);
    if (in4_)  mxFree(in4_);
    if (in5_)  mxFree(in5_);
    if (in6_)  mxFree(in6_);
    if (in7_)  mxFree(in7_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- nufft2d.mw: 86 ----
 * dirft2d2(int[1] nj, double[] xj, double[] yj, inout dcomplex[] cj, int[1] iflag, int[1] ms, int[1] mt, dcomplex[] fk);
 */
const char* stubids2_ = "dirft2d2(i int[x], i double[], i double[], io dcomplex[], i int[x], i int[x], i int[x], i dcomplex[])";

void mexStub2(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    int*        in0_ =0; /* nj         */
    double*     in1_ =0; /* xj         */
    double*     in2_ =0; /* yj         */
    dcomplex*   in3_ =0; /* cj         */
    int*        in4_ =0; /* iflag      */
    int*        in5_ =0; /* ms         */
    int*        in6_ =0; /* mt         */
    dcomplex*   in7_ =0; /* fk         */
    mwSize      dim8_;   /* 1          */
    mwSize      dim9_;   /* 1          */
    mwSize      dim10_;   /* 1          */
    mwSize      dim11_;   /* 1          */

    dim8_ = (mwSize) mxWrapGetScalar(prhs[8], &mw_err_txt_);
    dim9_ = (mwSize) mxWrapGetScalar(prhs[9], &mw_err_txt_);
    dim10_ = (mwSize) mxWrapGetScalar(prhs[10], &mw_err_txt_);
    dim11_ = (mwSize) mxWrapGetScalar(prhs[11], &mw_err_txt_);

    if (mxGetM(prhs[0])*mxGetN(prhs[0]) != dim8_) {
        mw_err_txt_ = "Bad argument size: nj";        goto mw_err_label;
    }

    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != dim9_) {
        mw_err_txt_ = "Bad argument size: iflag";        goto mw_err_label;
    }

    if (mxGetM(prhs[5])*mxGetN(prhs[5]) != dim10_) {
        mw_err_txt_ = "Bad argument size: ms";        goto mw_err_label;
    }

    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != dim11_) {
        mw_err_txt_ = "Bad argument size: mt";        goto mw_err_label;
    }

    if (mxGetM(prhs[0])*mxGetN(prhs[0]) != 0) {
        in0_ = mxWrapGetArray_int(prhs[0], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in0_ = NULL;
    if (mxGetM(prhs[1])*mxGetN(prhs[1]) != 0) {
        in1_ = mxGetPr(prhs[1]);
    } else
        in1_ = NULL;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxGetPr(prhs[2]);
    } else
        in2_ = NULL;
    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 0) {
        in3_ = mxWrapGetArray_dcomplex(prhs[3], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in3_ = NULL;
    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != 0) {
        in4_ = mxWrapGetArray_int(prhs[4], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in4_ = NULL;
    if (mxGetM(prhs[5])*mxGetN(prhs[5]) != 0) {
        in5_ = mxWrapGetArray_int(prhs[5], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in5_ = NULL;
    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != 0) {
        in6_ = mxWrapGetArray_int(prhs[6], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in6_ = NULL;
    if (mxGetM(prhs[7])*mxGetN(prhs[7]) != 0) {
        in7_ = mxWrapGetArray_dcomplex(prhs[7], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in7_ = NULL;
    if (mexprofrecord_)
        mexprofrecord_[2]++;
    MWF77_dirft2d2(in0_, in1_, in2_, in3_, in4_, in5_, in6_, in7_);
    plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[3]), mxGetN(prhs[3]), mxCOMPLEX);
    mxWrapCopy_dcomplex(plhs[0], in3_, mxGetM(prhs[3])*mxGetN(prhs[3]));

mw_err_label:
    if (in0_)  mxFree(in0_);
    if (in3_)  mxFree(in3_);
    if (in4_)  mxFree(in4_);
    if (in5_)  mxFree(in5_);
    if (in6_)  mxFree(in6_);
    if (in7_)  mxFree(in7_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- nufft2d.mw: 120 ----
 * dirft2d3(int[1] nj, double[] xj, double[] yj, dcomplex[] cj, int[1] iflag, int[1] nk, double[] sk, double[] tk, inout dcomplex[] fk);
 */
const char* stubids3_ = "dirft2d3(i int[x], i double[], i double[], i dcomplex[], i int[x], i int[x], i double[], i double[], io dcomplex[])";

void mexStub3(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    int*        in0_ =0; /* nj         */
    double*     in1_ =0; /* xj         */
    double*     in2_ =0; /* yj         */
    dcomplex*   in3_ =0; /* cj         */
    int*        in4_ =0; /* iflag      */
    int*        in5_ =0; /* nk         */
    double*     in6_ =0; /* sk         */
    double*     in7_ =0; /* tk         */
    dcomplex*   in8_ =0; /* fk         */
    mwSize      dim9_;   /* 1          */
    mwSize      dim10_;   /* 1          */
    mwSize      dim11_;   /* 1          */

    dim9_ = (mwSize) mxWrapGetScalar(prhs[9], &mw_err_txt_);
    dim10_ = (mwSize) mxWrapGetScalar(prhs[10], &mw_err_txt_);
    dim11_ = (mwSize) mxWrapGetScalar(prhs[11], &mw_err_txt_);

    if (mxGetM(prhs[0])*mxGetN(prhs[0]) != dim9_) {
        mw_err_txt_ = "Bad argument size: nj";        goto mw_err_label;
    }

    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != dim10_) {
        mw_err_txt_ = "Bad argument size: iflag";        goto mw_err_label;
    }

    if (mxGetM(prhs[5])*mxGetN(prhs[5]) != dim11_) {
        mw_err_txt_ = "Bad argument size: nk";        goto mw_err_label;
    }

    if (mxGetM(prhs[0])*mxGetN(prhs[0]) != 0) {
        in0_ = mxWrapGetArray_int(prhs[0], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in0_ = NULL;
    if (mxGetM(prhs[1])*mxGetN(prhs[1]) != 0) {
        in1_ = mxGetPr(prhs[1]);
    } else
        in1_ = NULL;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxGetPr(prhs[2]);
    } else
        in2_ = NULL;
    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 0) {
        in3_ = mxWrapGetArray_dcomplex(prhs[3], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in3_ = NULL;
    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != 0) {
        in4_ = mxWrapGetArray_int(prhs[4], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in4_ = NULL;
    if (mxGetM(prhs[5])*mxGetN(prhs[5]) != 0) {
        in5_ = mxWrapGetArray_int(prhs[5], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in5_ = NULL;
    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != 0) {
        in6_ = mxGetPr(prhs[6]);
    } else
        in6_ = NULL;
    if (mxGetM(prhs[7])*mxGetN(prhs[7]) != 0) {
        in7_ = mxGetPr(prhs[7]);
    } else
        in7_ = NULL;
    if (mxGetM(prhs[8])*mxGetN(prhs[8]) != 0) {
        in8_ = mxWrapGetArray_dcomplex(prhs[8], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in8_ = NULL;
    if (mexprofrecord_)
        mexprofrecord_[3]++;
    MWF77_dirft2d3(in0_, in1_, in2_, in3_, in4_, in5_, in6_, in7_, in8_);
    plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[8]), mxGetN(prhs[8]), mxCOMPLEX);
    mxWrapCopy_dcomplex(plhs[0], in8_, mxGetM(prhs[8])*mxGetN(prhs[8]));

mw_err_label:
    if (in0_)  mxFree(in0_);
    if (in3_)  mxFree(in3_);
    if (in4_)  mxFree(in4_);
    if (in5_)  mxFree(in5_);
    if (in8_)  mxFree(in8_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- nufft2d.mw: 162 ----
 * nufft2d1f90(int[1] nj, double[] xj, double[] yj, dcomplex[] cj, int[1] iflag, double[1] eps, int[1] ms, int[1] mt, inout dcomplex[] fk, inout int[1] ier);
 */
const char* stubids4_ = "nufft2d1f90(i int[x], i double[], i double[], i dcomplex[], i int[x], i double[x], i int[x], i int[x], io dcomplex[], io int[x])";

void mexStub4(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    int*        in0_ =0; /* nj         */
    double*     in1_ =0; /* xj         */
    double*     in2_ =0; /* yj         */
    dcomplex*   in3_ =0; /* cj         */
    int*        in4_ =0; /* iflag      */
    double*     in5_ =0; /* eps        */
    int*        in6_ =0; /* ms         */
    int*        in7_ =0; /* mt         */
    dcomplex*   in8_ =0; /* fk         */
    int*        in9_ =0; /* ier        */
    mwSize      dim10_;   /* 1          */
    mwSize      dim11_;   /* 1          */
    mwSize      dim12_;   /* 1          */
    mwSize      dim13_;   /* 1          */
    mwSize      dim14_;   /* 1          */
    mwSize      dim15_;   /* 1          */

    dim10_ = (mwSize) mxWrapGetScalar(prhs[10], &mw_err_txt_);
    dim11_ = (mwSize) mxWrapGetScalar(prhs[11], &mw_err_txt_);
    dim12_ = (mwSize) mxWrapGetScalar(prhs[12], &mw_err_txt_);
    dim13_ = (mwSize) mxWrapGetScalar(prhs[13], &mw_err_txt_);
    dim14_ = (mwSize) mxWrapGetScalar(prhs[14], &mw_err_txt_);
    dim15_ = (mwSize) mxWrapGetScalar(prhs[15], &mw_err_txt_);

    if (mxGetM(prhs[0])*mxGetN(prhs[0]) != dim10_) {
        mw_err_txt_ = "Bad argument size: nj";        goto mw_err_label;
    }

    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != dim11_) {
        mw_err_txt_ = "Bad argument size: iflag";        goto mw_err_label;
    }

    if (mxGetM(prhs[5])*mxGetN(prhs[5]) != dim12_) {
        mw_err_txt_ = "Bad argument size: eps";        goto mw_err_label;
    }

    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != dim13_) {
        mw_err_txt_ = "Bad argument size: ms";        goto mw_err_label;
    }

    if (mxGetM(prhs[7])*mxGetN(prhs[7]) != dim14_) {
        mw_err_txt_ = "Bad argument size: mt";        goto mw_err_label;
    }

    if (mxGetM(prhs[9])*mxGetN(prhs[9]) != dim15_) {
        mw_err_txt_ = "Bad argument size: ier";        goto mw_err_label;
    }

    if (mxGetM(prhs[0])*mxGetN(prhs[0]) != 0) {
        in0_ = mxWrapGetArray_int(prhs[0], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in0_ = NULL;
    if (mxGetM(prhs[1])*mxGetN(prhs[1]) != 0) {
        in1_ = mxGetPr(prhs[1]);
    } else
        in1_ = NULL;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxGetPr(prhs[2]);
    } else
        in2_ = NULL;
    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 0) {
        in3_ = mxWrapGetArray_dcomplex(prhs[3], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in3_ = NULL;
    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != 0) {
        in4_ = mxWrapGetArray_int(prhs[4], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in4_ = NULL;
    if (mxGetM(prhs[5])*mxGetN(prhs[5]) != 0) {
        in5_ = mxGetPr(prhs[5]);
    } else
        in5_ = NULL;
    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != 0) {
        in6_ = mxWrapGetArray_int(prhs[6], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in6_ = NULL;
    if (mxGetM(prhs[7])*mxGetN(prhs[7]) != 0) {
        in7_ = mxWrapGetArray_int(prhs[7], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in7_ = NULL;
    if (mxGetM(prhs[8])*mxGetN(prhs[8]) != 0) {
        in8_ = mxWrapGetArray_dcomplex(prhs[8], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in8_ = NULL;
    if (mxGetM(prhs[9])*mxGetN(prhs[9]) != 0) {
        in9_ = mxWrapGetArray_int(prhs[9], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in9_ = NULL;
    if (mexprofrecord_)
        mexprofrecord_[4]++;
    MWF77_nufft2d1f90(in0_, in1_, in2_, in3_, in4_, in5_, in6_, in7_, in8_, in9_);
    plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[8]), mxGetN(prhs[8]), mxCOMPLEX);
    mxWrapCopy_dcomplex(plhs[0], in8_, mxGetM(prhs[8])*mxGetN(prhs[8]));
    plhs[1] = mxCreateDoubleMatrix(dim15_, 1, mxREAL);
    mxWrapCopy_int(plhs[1], in9_, dim15_);

mw_err_label:
    if (in0_)  mxFree(in0_);
    if (in3_)  mxFree(in3_);
    if (in4_)  mxFree(in4_);
    if (in6_)  mxFree(in6_);
    if (in7_)  mxFree(in7_);
    if (in8_)  mxFree(in8_);
    if (in9_)  mxFree(in9_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- nufft2d.mw: 201 ----
 * nufft2d2f90(int[1] nj, double[] xj, double[] yj, inout dcomplex[] cj, int[1] iflag, double[1] eps, int[1] ms, int[1] mt, dcomplex[] fk, inout int[1] ier);
 */
const char* stubids5_ = "nufft2d2f90(i int[x], i double[], i double[], io dcomplex[], i int[x], i double[x], i int[x], i int[x], i dcomplex[], io int[x])";

void mexStub5(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    int*        in0_ =0; /* nj         */
    double*     in1_ =0; /* xj         */
    double*     in2_ =0; /* yj         */
    dcomplex*   in3_ =0; /* cj         */
    int*        in4_ =0; /* iflag      */
    double*     in5_ =0; /* eps        */
    int*        in6_ =0; /* ms         */
    int*        in7_ =0; /* mt         */
    dcomplex*   in8_ =0; /* fk         */
    int*        in9_ =0; /* ier        */
    mwSize      dim10_;   /* 1          */
    mwSize      dim11_;   /* 1          */
    mwSize      dim12_;   /* 1          */
    mwSize      dim13_;   /* 1          */
    mwSize      dim14_;   /* 1          */
    mwSize      dim15_;   /* 1          */

    dim10_ = (mwSize) mxWrapGetScalar(prhs[10], &mw_err_txt_);
    dim11_ = (mwSize) mxWrapGetScalar(prhs[11], &mw_err_txt_);
    dim12_ = (mwSize) mxWrapGetScalar(prhs[12], &mw_err_txt_);
    dim13_ = (mwSize) mxWrapGetScalar(prhs[13], &mw_err_txt_);
    dim14_ = (mwSize) mxWrapGetScalar(prhs[14], &mw_err_txt_);
    dim15_ = (mwSize) mxWrapGetScalar(prhs[15], &mw_err_txt_);

    if (mxGetM(prhs[0])*mxGetN(prhs[0]) != dim10_) {
        mw_err_txt_ = "Bad argument size: nj";        goto mw_err_label;
    }

    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != dim11_) {
        mw_err_txt_ = "Bad argument size: iflag";        goto mw_err_label;
    }

    if (mxGetM(prhs[5])*mxGetN(prhs[5]) != dim12_) {
        mw_err_txt_ = "Bad argument size: eps";        goto mw_err_label;
    }

    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != dim13_) {
        mw_err_txt_ = "Bad argument size: ms";        goto mw_err_label;
    }

    if (mxGetM(prhs[7])*mxGetN(prhs[7]) != dim14_) {
        mw_err_txt_ = "Bad argument size: mt";        goto mw_err_label;
    }

    if (mxGetM(prhs[9])*mxGetN(prhs[9]) != dim15_) {
        mw_err_txt_ = "Bad argument size: ier";        goto mw_err_label;
    }

    if (mxGetM(prhs[0])*mxGetN(prhs[0]) != 0) {
        in0_ = mxWrapGetArray_int(prhs[0], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in0_ = NULL;
    if (mxGetM(prhs[1])*mxGetN(prhs[1]) != 0) {
        in1_ = mxGetPr(prhs[1]);
    } else
        in1_ = NULL;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxGetPr(prhs[2]);
    } else
        in2_ = NULL;
    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 0) {
        in3_ = mxWrapGetArray_dcomplex(prhs[3], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in3_ = NULL;
    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != 0) {
        in4_ = mxWrapGetArray_int(prhs[4], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in4_ = NULL;
    if (mxGetM(prhs[5])*mxGetN(prhs[5]) != 0) {
        in5_ = mxGetPr(prhs[5]);
    } else
        in5_ = NULL;
    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != 0) {
        in6_ = mxWrapGetArray_int(prhs[6], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in6_ = NULL;
    if (mxGetM(prhs[7])*mxGetN(prhs[7]) != 0) {
        in7_ = mxWrapGetArray_int(prhs[7], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in7_ = NULL;
    if (mxGetM(prhs[8])*mxGetN(prhs[8]) != 0) {
        in8_ = mxWrapGetArray_dcomplex(prhs[8], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in8_ = NULL;
    if (mxGetM(prhs[9])*mxGetN(prhs[9]) != 0) {
        in9_ = mxWrapGetArray_int(prhs[9], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in9_ = NULL;
    if (mexprofrecord_)
        mexprofrecord_[5]++;
    MWF77_nufft2d2f90(in0_, in1_, in2_, in3_, in4_, in5_, in6_, in7_, in8_, in9_);
    plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[3]), mxGetN(prhs[3]), mxCOMPLEX);
    mxWrapCopy_dcomplex(plhs[0], in3_, mxGetM(prhs[3])*mxGetN(prhs[3]));
    plhs[1] = mxCreateDoubleMatrix(dim15_, 1, mxREAL);
    mxWrapCopy_int(plhs[1], in9_, dim15_);

mw_err_label:
    if (in0_)  mxFree(in0_);
    if (in3_)  mxFree(in3_);
    if (in4_)  mxFree(in4_);
    if (in6_)  mxFree(in6_);
    if (in7_)  mxFree(in7_);
    if (in8_)  mxFree(in8_);
    if (in9_)  mxFree(in9_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ---- nufft2d.mw: 241 ----
 * nufft2d3f90(int[1] nj, double[] xj, double[] yj, dcomplex[] cj, int[1] iflag, double[1] eps, int[1] nk, double[] sk, double[] tk, inout dcomplex[] fk, inout int[1] ier);
 */
const char* stubids6_ = "nufft2d3f90(i int[x], i double[], i double[], i dcomplex[], i int[x], i double[x], i int[x], i double[], i double[], io dcomplex[], io int[x])";

void mexStub6(int nlhs, mxArray* plhs[],
              int nrhs, const mxArray* prhs[])
{
    const char* mw_err_txt_ = 0;
    int*        in0_ =0; /* nj         */
    double*     in1_ =0; /* xj         */
    double*     in2_ =0; /* yj         */
    dcomplex*   in3_ =0; /* cj         */
    int*        in4_ =0; /* iflag      */
    double*     in5_ =0; /* eps        */
    int*        in6_ =0; /* nk         */
    double*     in7_ =0; /* sk         */
    double*     in8_ =0; /* tk         */
    dcomplex*   in9_ =0; /* fk         */
    int*        in10_ =0; /* ier        */
    mwSize      dim11_;   /* 1          */
    mwSize      dim12_;   /* 1          */
    mwSize      dim13_;   /* 1          */
    mwSize      dim14_;   /* 1          */
    mwSize      dim15_;   /* 1          */

    dim11_ = (mwSize) mxWrapGetScalar(prhs[11], &mw_err_txt_);
    dim12_ = (mwSize) mxWrapGetScalar(prhs[12], &mw_err_txt_);
    dim13_ = (mwSize) mxWrapGetScalar(prhs[13], &mw_err_txt_);
    dim14_ = (mwSize) mxWrapGetScalar(prhs[14], &mw_err_txt_);
    dim15_ = (mwSize) mxWrapGetScalar(prhs[15], &mw_err_txt_);

    if (mxGetM(prhs[0])*mxGetN(prhs[0]) != dim11_) {
        mw_err_txt_ = "Bad argument size: nj";        goto mw_err_label;
    }

    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != dim12_) {
        mw_err_txt_ = "Bad argument size: iflag";        goto mw_err_label;
    }

    if (mxGetM(prhs[5])*mxGetN(prhs[5]) != dim13_) {
        mw_err_txt_ = "Bad argument size: eps";        goto mw_err_label;
    }

    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != dim14_) {
        mw_err_txt_ = "Bad argument size: nk";        goto mw_err_label;
    }

    if (mxGetM(prhs[10])*mxGetN(prhs[10]) != dim15_) {
        mw_err_txt_ = "Bad argument size: ier";        goto mw_err_label;
    }

    if (mxGetM(prhs[0])*mxGetN(prhs[0]) != 0) {
        in0_ = mxWrapGetArray_int(prhs[0], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in0_ = NULL;
    if (mxGetM(prhs[1])*mxGetN(prhs[1]) != 0) {
        in1_ = mxGetPr(prhs[1]);
    } else
        in1_ = NULL;
    if (mxGetM(prhs[2])*mxGetN(prhs[2]) != 0) {
        in2_ = mxGetPr(prhs[2]);
    } else
        in2_ = NULL;
    if (mxGetM(prhs[3])*mxGetN(prhs[3]) != 0) {
        in3_ = mxWrapGetArray_dcomplex(prhs[3], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in3_ = NULL;
    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != 0) {
        in4_ = mxWrapGetArray_int(prhs[4], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in4_ = NULL;
    if (mxGetM(prhs[5])*mxGetN(prhs[5]) != 0) {
        in5_ = mxGetPr(prhs[5]);
    } else
        in5_ = NULL;
    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != 0) {
        in6_ = mxWrapGetArray_int(prhs[6], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in6_ = NULL;
    if (mxGetM(prhs[7])*mxGetN(prhs[7]) != 0) {
        in7_ = mxGetPr(prhs[7]);
    } else
        in7_ = NULL;
    if (mxGetM(prhs[8])*mxGetN(prhs[8]) != 0) {
        in8_ = mxGetPr(prhs[8]);
    } else
        in8_ = NULL;
    if (mxGetM(prhs[9])*mxGetN(prhs[9]) != 0) {
        in9_ = mxWrapGetArray_dcomplex(prhs[9], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in9_ = NULL;
    if (mxGetM(prhs[10])*mxGetN(prhs[10]) != 0) {
        in10_ = mxWrapGetArray_int(prhs[10], &mw_err_txt_);
        if (mw_err_txt_)
            goto mw_err_label;
    } else
        in10_ = NULL;
    if (mexprofrecord_)
        mexprofrecord_[6]++;
    MWF77_nufft2d3f90(in0_, in1_, in2_, in3_, in4_, in5_, in6_, in7_, in8_, in9_, in10_);
    plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[9]), mxGetN(prhs[9]), mxCOMPLEX);
    mxWrapCopy_dcomplex(plhs[0], in9_, mxGetM(prhs[9])*mxGetN(prhs[9]));
    plhs[1] = mxCreateDoubleMatrix(dim15_, 1, mxREAL);
    mxWrapCopy_int(plhs[1], in10_, dim15_);

mw_err_label:
    if (in0_)  mxFree(in0_);
    if (in3_)  mxFree(in3_);
    if (in4_)  mxFree(in4_);
    if (in6_)  mxFree(in6_);
    if (in9_)  mxFree(in9_);
    if (in10_)  mxFree(in10_);
    if (mw_err_txt_)
        mexErrMsgTxt(mw_err_txt_);
}

/* ----
 */
void mexFunction(int nlhs, mxArray* plhs[],
                 int nrhs, const mxArray* prhs[])
{
    char id[512];
    if (nrhs == 0) {
        mexPrintf("Mex function installed\n");
        return;
    }

    if (mxGetString(prhs[0], id, sizeof(id)) != 0)
        mexErrMsgTxt("Identifier should be a string");
    else if (strcmp(id, stubids1_) == 0)
        mexStub1(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids2_) == 0)
        mexStub2(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids3_) == 0)
        mexStub3(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids4_) == 0)
        mexStub4(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids5_) == 0)
        mexStub5(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, stubids6_) == 0)
        mexStub6(nlhs,plhs, nrhs-1,prhs+1);
    else if (strcmp(id, "*profile on*") == 0) {
        if (!mexprofrecord_) {
            mexprofrecord_ = (int*) malloc(7 * sizeof(int));
            mexLock();
        }
        memset(mexprofrecord_, 0, 7 * sizeof(int));
    } else if (strcmp(id, "*profile off*") == 0) {
        if (mexprofrecord_) {
            free(mexprofrecord_);
            mexUnlock();
        }
        mexprofrecord_ = NULL;
    } else if (strcmp(id, "*profile report*") == 0) {
        if (!mexprofrecord_)
            mexPrintf("Profiler inactive\n");
        mexPrintf("%d calls to nufft2d.mw:52\n", mexprofrecord_[1]);
        mexPrintf("%d calls to nufft2d.mw:86\n", mexprofrecord_[2]);
        mexPrintf("%d calls to nufft2d.mw:120\n", mexprofrecord_[3]);
        mexPrintf("%d calls to nufft2d.mw:162\n", mexprofrecord_[4]);
        mexPrintf("%d calls to nufft2d.mw:201\n", mexprofrecord_[5]);
        mexPrintf("%d calls to nufft2d.mw:241\n", mexprofrecord_[6]);
    } else if (strcmp(id, "*profile log*") == 0) {
        FILE* logfp;
        if (nrhs != 2 || mxGetString(prhs[1], id, sizeof(id)) != 0)
            mexErrMsgTxt("Must have two string arguments");
        logfp = fopen(id, "w+");
        if (!logfp)
            mexErrMsgTxt("Cannot open log for output");
        if (!mexprofrecord_)
            fprintf(logfp, "Profiler inactive\n");
        fprintf(logfp, "%d calls to nufft2d.mw:52\n", mexprofrecord_[1]);
        fprintf(logfp, "%d calls to nufft2d.mw:86\n", mexprofrecord_[2]);
        fprintf(logfp, "%d calls to nufft2d.mw:120\n", mexprofrecord_[3]);
        fprintf(logfp, "%d calls to nufft2d.mw:162\n", mexprofrecord_[4]);
        fprintf(logfp, "%d calls to nufft2d.mw:201\n", mexprofrecord_[5]);
        fprintf(logfp, "%d calls to nufft2d.mw:241\n", mexprofrecord_[6]);
        fclose(logfp);
    } else
        mexErrMsgTxt("Unknown identifier");
}

