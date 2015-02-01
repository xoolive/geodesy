/*
 * Freely inspiring myself from Cephes library
 */

#ifndef __INTEL_COMPILER

#  ifndef TOOLS_SSE2_MATH
#  define TOOLS_SSE2_MATH

#  include <emmintrin.h>

// Microsoft compiler
#  if defined(_MSC_VER)
#    define ALIGN16_BEG __declspec(align(16))
#    define ALIGN16_END
#  elif defined(__GNUC__)
#    define ALIGN16_BEG
#    define ALIGN16_END __attribute__((aligned (16)))
#  endif

static const ALIGN16_BEG int _mm_cst_sign_mask_ps[4] ALIGN16_END =
{ 0x80000000, 0x80000000, 0x80000000, 0x80000000 };
static const ALIGN16_BEG int _mm_cst_inv_sign_mask_ps[4] ALIGN16_END =
{ ~0x80000000, ~0x80000000, ~0x80000000, ~0x80000000 };
static const ALIGN16_BEG int _mm_cst_sign_mask_pd[4] ALIGN16_END =
{ 0, 0x80000000, 0, 0x80000000 };
static const ALIGN16_BEG int _mm_cst_inv_sign_mask_pd[4] ALIGN16_END =
{ ~0, ~0x80000000, ~0, ~0x80000000 };

static const ALIGN16_BEG int _mm_cst_one[4] ALIGN16_END = { 1, 1, 1, 1 };
static const ALIGN16_BEG int _mm_cst_inv1[4] ALIGN16_END = { ~1, ~1, ~1, ~1 };
static const ALIGN16_BEG int _mm_cst_two[4] ALIGN16_END = { 2, 2, 2, 2 };
static const ALIGN16_BEG int _mm_cst_four[4] ALIGN16_END = { 4, 4, 4, 4 };

static const __m128 _mm_cst_ps_one   = _mm_set1_ps(1.0f);
static const __m128 _mm_cst_ps_mone  = _mm_set1_ps(-1.0f);
static const __m128 _mm_cst_ps_mtwo  = _mm_set1_ps(-2.0f);
static const __m128 _mm_cst_ps_0p5   = _mm_set1_ps(0.5f);

static const __m128d _mm_cst_pd_one   = _mm_set1_pd(1.0);
static const __m128d _mm_cst_pd_mone  = _mm_set1_pd(-1.0);
static const __m128d _mm_cst_pd_mtwo  = _mm_set1_pd(-2.0);
static const __m128d _mm_cst_pd_0p5   = _mm_set1_pd(0.5);

static const __m128 _mm_cst_ps_tan3pio8 = _mm_set1_ps(2.414213562373095f);
static const __m128 _mm_cst_ps_tanpio8  = _mm_set1_ps(0.4142135623730950f);

static const unsigned short T3P8[]       = { 0x9de6, 0x333f, 0x504f, 0x4003};
static const __m128d _mm_cst_pd_tan3pio8 = _mm_set1_pd(*(double*) T3P8);
static const __m128d _mm_cst_pd_tanpio8  = _mm_set1_pd(0.66); // see atan.c

static const __m128 _mm_cst_ps_fopi     = _mm_set1_ps(1.27323954473516f);

static const __m128 _mm_cst_ps_pi    = _mm_set1_ps(3.14159265358979f);
static const __m128 _mm_cst_ps_mpi   = _mm_set1_ps(-3.14159265358979f);
static const __m128 _mm_cst_ps_pio2  = _mm_set1_ps(1.5707963267948966f);
static const __m128 _mm_cst_ps_mpio2 = _mm_set1_ps(-1.5707963267948966f);
static const __m128 _mm_cst_ps_pio4  = _mm_set1_ps(0.7853981633974483f);

static const unsigned short PI[4]     = { 0x2d18, 0x5444, 0x21fb, 0x4009 };
static const unsigned short PIO2[4]   = { 0x2d18, 0x5444, 0x21fb, 0x3ff9 };
static const unsigned short PIO4[4]   = { 0x2d18, 0x5444, 0x21fb, 0x3fe9 };

static const __m128d _mm_cst_pd_pi    = _mm_set1_pd(*(double*) PI);
static const __m128d _mm_cst_pd_mpi   = _mm_set1_pd(0. - *(double*) PI);
static const __m128d _mm_cst_pd_pio2  = _mm_set1_pd(*(double*) PIO2);
static const __m128d _mm_cst_pd_mpio2 = _mm_set1_pd(0. - *(double*) PIO2);
static const __m128d _mm_cst_pd_pio4  = _mm_set1_pd(*(double*) PIO4);

/*
 * Extended precision modular arithmetic:
 *
 * within bc -l:
 *    obase = 16
 *    a(1)
 *    .C90FDAA22168C234A
 *
 * pi/4 is then decomposed:
 *    DP1 + DP2 + DP3 = pi / 4 (and sign)
 *
 * See all binary representation DP1: .C900 [00...]
 *                               DP2: .000FDA00 [00...]
 *                               DP3: .000000A22168 [etc.]
 *
 *  ->  Note the alignment
 */
static const __m128 _mm_cst_ps_DP1 = _mm_set1_ps(0.78515625f);
static const __m128 _mm_cst_ps_DP2 = _mm_set1_ps(2.4187564849853515625e-4f);
static const __m128 _mm_cst_ps_DP3 = _mm_set1_ps(3.77489497744594108e-8f);

static unsigned short P1[] = { 0x0000, 0x4000, 0x21fb, 0x3fe9 };
static unsigned short P2[] = { 0x0000, 0x0000, 0x442d, 0x3e64 };
static unsigned short P3[] = { 0x5170, 0x98cc, 0x4698, 0x3ce8 };

static const __m128d _mm_cst_pd_DP1 = _mm_set1_pd(*(double*) P1);
static const __m128d _mm_cst_pd_DP2 = _mm_set1_pd(*(double*) P2);
static const __m128d _mm_cst_pd_DP3 = _mm_set1_pd(*(double*) P3);

static const __m128 _mm_cst_ps_sincof_p0 = _mm_set1_ps(-1.9515295891e-4f);
static const __m128 _mm_cst_ps_sincof_p1 = _mm_set1_ps( 8.3321608736e-3f);
static const __m128 _mm_cst_ps_sincof_p2 = _mm_set1_ps(-1.6666654611e-1f);

static unsigned short sincof0[] = { 0x9ccd, 0x1fd1, 0xd8fd, 0x3de5 };
static unsigned short sincof1[] = { 0x1f5d, 0xa929, 0xe5e5, 0xbe5a };
static unsigned short sincof2[] = { 0x48a1, 0x567d, 0x1de3, 0x3ec7 };
static unsigned short sincof3[] = { 0xdf03, 0x19bf, 0x01a0, 0xbf2a };
static unsigned short sincof4[] = { 0xf7d0, 0x1110, 0x1111, 0x3f81 };
static unsigned short sincof5[] = { 0x5548, 0x5555, 0x5555, 0xbfc5 };

static const __m128d _mm_cst_pd_sincof_p0 = _mm_set1_pd(*(double*) sincof0);
static const __m128d _mm_cst_pd_sincof_p1 = _mm_set1_pd(*(double*) sincof1);
static const __m128d _mm_cst_pd_sincof_p2 = _mm_set1_pd(*(double*) sincof2);
static const __m128d _mm_cst_pd_sincof_p3 = _mm_set1_pd(*(double*) sincof3);
static const __m128d _mm_cst_pd_sincof_p4 = _mm_set1_pd(*(double*) sincof4);
static const __m128d _mm_cst_pd_sincof_p5 = _mm_set1_pd(*(double*) sincof5);

static const __m128 _mm_cst_ps_coscof_p0 = _mm_set1_ps( 2.443315711809948e-5f);
static const __m128 _mm_cst_ps_coscof_p1 = _mm_set1_ps(-1.388731625493765e-3f);
static const __m128 _mm_cst_ps_coscof_p2 = _mm_set1_ps( 4.166664568298827e-2f);


static unsigned short coscof0[] = { 0x1a9b, 0xa086, 0xfa49, 0xbda8 };
static unsigned short coscof1[] = { 0x3f05, 0x7b4e, 0xee9d, 0x3e21 };
static unsigned short coscof2[] = { 0x4bc6, 0x7eac, 0x7e4f, 0xbe92 };
static unsigned short coscof3[] = { 0x44f5, 0x19c8, 0x01a0, 0x3efa };
static unsigned short coscof4[] = { 0x4f91, 0x16c1, 0xc16c, 0xbf56 };
static unsigned short coscof5[] = { 0x554b, 0x5555, 0x5555, 0x3fa5 };

static const __m128d _mm_cst_pd_coscof_p0 = _mm_set1_pd(*(double*) coscof0);
static const __m128d _mm_cst_pd_coscof_p1 = _mm_set1_pd(*(double*) coscof1);
static const __m128d _mm_cst_pd_coscof_p2 = _mm_set1_pd(*(double*) coscof2);
static const __m128d _mm_cst_pd_coscof_p3 = _mm_set1_pd(*(double*) coscof3);
static const __m128d _mm_cst_pd_coscof_p4 = _mm_set1_pd(*(double*) coscof4);
static const __m128d _mm_cst_pd_coscof_p5 = _mm_set1_pd(*(double*) coscof5);

static const __m128 _mm_cst_ps_asincof_p0 = _mm_set1_ps(4.2163199048e-2f);
static const __m128 _mm_cst_ps_asincof_p1 = _mm_set1_ps(2.4181311049e-2f);
static const __m128 _mm_cst_ps_asincof_p2 = _mm_set1_ps(4.5470025998e-2f);
static const __m128 _mm_cst_ps_asincof_p3 = _mm_set1_ps(7.4953002686e-2f);
static const __m128 _mm_cst_ps_asincof_p4 = _mm_set1_ps(1.6666752422e-1f);

static const short asincof_p0[] = { 0x8ad3, 0x0bd4, 0x6b9b, 0x3f71 };
static const short asincof_p1[] = { 0x5c16, 0x333e, 0x4341, 0xbfe3 };
static const short asincof_p2[] = { 0x2dd9, 0x178a, 0xc74b, 0x4015 };
static const short asincof_p3[] = { 0x907b, 0xde27, 0x4331, 0xc030 };
static const short asincof_p4[] = { 0x9259, 0xda77, 0x9007, 0x4033 };
static const short asincof_p5[] = { 0xafd5, 0x06ce, 0x656c, 0xc020 };

static const __m128d _mm_cst_pd_asincof_p0 = _mm_set1_pd(*(double*) asincof_p0);
static const __m128d _mm_cst_pd_asincof_p1 = _mm_set1_pd(*(double*) asincof_p1);
static const __m128d _mm_cst_pd_asincof_p2 = _mm_set1_pd(*(double*) asincof_p2);
static const __m128d _mm_cst_pd_asincof_p3 = _mm_set1_pd(*(double*) asincof_p3);
static const __m128d _mm_cst_pd_asincof_p4 = _mm_set1_pd(*(double*) asincof_p4);
static const __m128d _mm_cst_pd_asincof_p5 = _mm_set1_pd(*(double*) asincof_p5);

static const short asincof_q0[] = { 0x0eab, 0x0b5e, 0x7b59, 0xc02d };
static const short asincof_q1[] = { 0x9054, 0x25fe, 0x9fc0, 0x4051 };
static const short asincof_q2[] = { 0x76d7, 0x6d35, 0x65bb, 0xc062 };
static const short asincof_q3[] = { 0xbf9d, 0x84ff, 0x7056, 0x4061 };
static const short asincof_q4[] = { 0x07ac, 0x0a36, 0x9822, 0xc048 };

static const __m128d _mm_cst_pd_asincof_q0 = _mm_set1_pd(*(double*) asincof_q0);
static const __m128d _mm_cst_pd_asincof_q1 = _mm_set1_pd(*(double*) asincof_q1);
static const __m128d _mm_cst_pd_asincof_q2 = _mm_set1_pd(*(double*) asincof_q2);
static const __m128d _mm_cst_pd_asincof_q3 = _mm_set1_pd(*(double*) asincof_q3);
static const __m128d _mm_cst_pd_asincof_q4 = _mm_set1_pd(*(double*) asincof_q4);


static const short asincof_r0[] = { 0x9f08, 0x988e, 0x4fc3, 0x3f68 };
static const short asincof_r1[] = { 0x290f, 0x59f9, 0x0792, 0xbfe2 };
static const short asincof_r2[] = { 0x3e6a, 0xbaf3, 0xdff5, 0x401b };
static const short asincof_r3[] = { 0xab68, 0xac01, 0x91aa, 0xc039 };
static const short asincof_r4[] = { 0x081d, 0x40f3, 0x8962, 0x403c };

static const __m128d _mm_cst_pd_asincof_r0 = _mm_set1_pd(*(double*) asincof_r0);
static const __m128d _mm_cst_pd_asincof_r1 = _mm_set1_pd(*(double*) asincof_r1);
static const __m128d _mm_cst_pd_asincof_r2 = _mm_set1_pd(*(double*) asincof_r2);
static const __m128d _mm_cst_pd_asincof_r3 = _mm_set1_pd(*(double*) asincof_r3);
static const __m128d _mm_cst_pd_asincof_r4 = _mm_set1_pd(*(double*) asincof_r4);


static const short asincof_s0[] = { 0x5d8c, 0xb6bf, 0xf2a2, 0xc035 };
static const short asincof_s1[] = { 0x7f42, 0xaf6a, 0x6219, 0x4062 };
static const short asincof_s2[] = { 0x63ee, 0x9590, 0xfe08, 0xc077 };
static const short asincof_s3[] = { 0x44be, 0xb0b6, 0x6709, 0x4075 };

static const __m128d _mm_cst_pd_asincof_s0 = _mm_set1_pd(*(double*) asincof_s0);
static const __m128d _mm_cst_pd_asincof_s1 = _mm_set1_pd(*(double*) asincof_s1);
static const __m128d _mm_cst_pd_asincof_s2 = _mm_set1_pd(*(double*) asincof_s2);
static const __m128d _mm_cst_pd_asincof_s3 = _mm_set1_pd(*(double*) asincof_s3);


static const __m128 _mm_cst_ps_atancof_p0 = _mm_set1_ps( 8.05374449538e-2f);
static const __m128 _mm_cst_ps_atancof_p1 = _mm_set1_ps(-1.38776856032e-1f);
static const __m128 _mm_cst_ps_atancof_p2 = _mm_set1_ps( 1.99777106478e-1f);
static const __m128 _mm_cst_ps_atancof_p3 = _mm_set1_ps(-3.33329491539e-1f);

static const short atancof_p0[] = { 0x2594, 0xa1f7, 0x007f, 0xbfec };
static const short atancof_p1[] = { 0x807a, 0x5b6b, 0x2854, 0xc030 };
static const short atancof_p2[] = { 0x0273, 0x3688, 0xc08c, 0xc052 };
static const short atancof_p3[] = { 0xba25, 0x2d05, 0xb8bf, 0xc05e };
static const short atancof_p4[] = { 0xec8e, 0xfd28, 0x3669, 0xc050 };

static const __m128d _mm_cst_pd_atancof_p0 = _mm_set1_pd(*(double*) atancof_p0);
static const __m128d _mm_cst_pd_atancof_p1 = _mm_set1_pd(*(double*) atancof_p1);
static const __m128d _mm_cst_pd_atancof_p2 = _mm_set1_pd(*(double*) atancof_p2);
static const __m128d _mm_cst_pd_atancof_p3 = _mm_set1_pd(*(double*) atancof_p3);
static const __m128d _mm_cst_pd_atancof_p4 = _mm_set1_pd(*(double*) atancof_p4);

static const short atancof_q0[] = { 0x603c, 0x5b14, 0xdbc4, 0x4038 };
static const short atancof_q1[] = { 0xfa25, 0x43b8, 0xa0dd, 0x4064 };
static const short atancof_q2[] = { 0xbe3b, 0xd2e2, 0x0e18, 0x407b };
static const short atancof_q3[] = { 0x49ea, 0x13b0, 0x563f, 0x407e };
static const short atancof_q4[] = { 0x62ec, 0xfbbd, 0x519e, 0x4068 };


static const __m128d _mm_cst_pd_atancof_q0 = _mm_set1_pd(*(double*) atancof_q0);
static const __m128d _mm_cst_pd_atancof_q1 = _mm_set1_pd(*(double*) atancof_q1);
static const __m128d _mm_cst_pd_atancof_q2 = _mm_set1_pd(*(double*) atancof_q2);
static const __m128d _mm_cst_pd_atancof_q3 = _mm_set1_pd(*(double*) atancof_q3);
static const __m128d _mm_cst_pd_atancof_q4 = _mm_set1_pd(*(double*) atancof_q4);

/*
 * Range reduction is into intervals of pi/4.  The reduction error is nearly
 * eliminated by contriving an extended precision modular arithmetic.
 *
 * Two polynomial approximating functions are employed.
 * Between 0 and pi/4 the sine is approximated by
 *     x  +  x**3 P(x**2).
 * Between pi/4 and pi/2 the cosine is represented as
 *     1  -  x**2 Q(x**2).
 *
 * Partial loss of accuracy begins to occur at x = 2**30 = 1.074e9. The loss
 * is not gradual, but jumps suddenly to about 1 part in 10e7. Results may be
 * meaningless for x > 2**49 = 5.6e14.  The routine as implemented flags a
 * TLOSS error for x > 2**30 and returns 0.0.
 *
 */
__m128 _mm_sin_ps(__m128 x) {
  // Extract the sign bit and work with absolute values
  __m128 signbit = _mm_and_ps(x, *(__m128*) _mm_cst_sign_mask_ps);
  x = _mm_and_ps(x, *(__m128*) _mm_cst_inv_sign_mask_ps);

  // Cephes method for scaling between 0 et pi/4
  __m128 y = _mm_mul_ps(x, _mm_cst_ps_fopi);
  __m128i yf = _mm_cvttps_epi32(y); // floor
  // see j = (j+1) & (~1) in Cephes
  yf = _mm_add_epi32(yf, *(__m128i*) _mm_cst_one);
  yf = _mm_and_si128(yf, *(__m128i*) _mm_cst_inv1);
  y = _mm_cvtepi32_ps(yf);

  // Swap sign flag
  __m128i flag = _mm_and_si128(yf, *(__m128i*) _mm_cst_four);
  flag = _mm_slli_epi32(flag, 29); // flag << 29

  // Polynom selection mask
  yf = _mm_and_si128(yf, *(__m128i*) _mm_cst_two);
  yf = _mm_cmpeq_epi32(yf, _mm_setzero_si128());

  __m128 swapsign = _mm_castsi128_ps(flag);
  __m128 polymask = _mm_castsi128_ps(yf);
  signbit = _mm_xor_ps(signbit, swapsign);

  // Extended precision modular arithmetic (Cody and Waite)
  // z = ((x - y * DP1) - y * DP2) - y * DP3
  __m128 xmm1 = _mm_mul_ps(y, _mm_cst_ps_DP1);
  __m128 xmm2 = _mm_mul_ps(y, _mm_cst_ps_DP2);
  __m128 xmm3 = _mm_mul_ps(y, _mm_cst_ps_DP3);
  x = _mm_sub_ps(x, xmm1);
  x = _mm_sub_ps(x, xmm2);
  x = _mm_sub_ps(x, xmm3);

  // First polynom x \in [0, pi/4]
  __m128 z = _mm_mul_ps(x, x);

  y = _mm_mul_ps(_mm_cst_ps_coscof_p0, z);
  y = _mm_add_ps(y, _mm_cst_ps_coscof_p1);
  y = _mm_mul_ps(y, z);
  y = _mm_add_ps(y, _mm_cst_ps_coscof_p2);
  y = _mm_mul_ps(y, z);
  y = _mm_mul_ps(y, z);
  __m128 tmp = _mm_mul_ps(z, _mm_cst_ps_0p5);
  y = _mm_sub_ps(y, tmp);
  y = _mm_add_ps(y, _mm_cst_ps_one);

  // Second polynom
  __m128 y2 = _mm_mul_ps(_mm_cst_ps_sincof_p0, z);
  y2 = _mm_add_ps(y2, _mm_cst_ps_sincof_p1);
  y2 = _mm_mul_ps(y2, z);
  y2 = _mm_add_ps(y2, _mm_cst_ps_sincof_p2);
  y2 = _mm_mul_ps(y2, z);
  y2 = _mm_mul_ps(y2, x);
  y2 = _mm_add_ps(y2, x);

  // Select the correct result
  y2 = _mm_and_ps(polymask, y2);
  y = _mm_andnot_ps(polymask, y);
  y = _mm_add_ps(y, y2);

  // Update the sign
  y = _mm_xor_ps(y, signbit);

  return y;
}

__m128d _mm_sin_pd(__m128d x) {
  // Extract the sign bit and work with absolute values
  __m128d signbit = _mm_and_pd(x, *(__m128d*) _mm_cst_sign_mask_pd);
  x = _mm_and_pd(x, *(__m128d*) _mm_cst_inv_sign_mask_pd);

  // Cephes method for scaling between 0 et pi/4
  __m128d y = _mm_div_pd(x, _mm_cst_pd_pio4);
  __m64 yf = _mm_cvttpd_pi32(y); // floor
  // see j = (j+1) & (~1) in Cephes
  yf = _mm_add_pi32(yf, *(__m64*) _mm_cst_one);
  yf = _mm_and_si64(yf, *(__m64*) _mm_cst_inv1);
  y = _mm_cvtpi32_pd(yf);

  // Swap sign flag
  __m64 flag = _mm_and_si64(yf, *(__m64*) _mm_cst_four);
  flag = _mm_slli_pi32(flag, 29); // flag << 29

  // Polynom selection mask
  yf = _mm_and_si64(yf, *(__m64*) _mm_cst_two);
  yf = _mm_cmpeq_pi32(yf, _mm_setzero_si64());

  static ALIGN16_BEG int signs[4] ALIGN16_END;
  signs[0] = 0;
  signs[1] = ((int*)(&flag))[0];
  signs[2] = 0;
  signs[3] = ((int*)(&flag))[1];

  static ALIGN16_BEG int poly[4] ALIGN16_END;
  poly[0] = ((int*)(&yf))[0];
  poly[1] = ((int*)(&yf))[0];
  poly[2] = ((int*)(&yf))[1];
  poly[3] = ((int*)(&yf))[1];

  __m128d swapsign = _mm_load_pd((double*) signs);
  __m128d polymask = _mm_load_pd((double*) poly);
  signbit = _mm_xor_pd(signbit, swapsign);

  // Extended precision modular arithmetic (Cody and Waite)
  // z = ((x - y * DP1) - y * DP2) - y * DP3
  __m128d xmm1 = _mm_mul_pd(y, _mm_cst_pd_DP1);
  __m128d xmm2 = _mm_mul_pd(y, _mm_cst_pd_DP2);
  __m128d xmm3 = _mm_mul_pd(y, _mm_cst_pd_DP3);
  x = _mm_sub_pd(x, xmm1);
  x = _mm_sub_pd(x, xmm2);
  x = _mm_sub_pd(x, xmm3);

  // First polynom x \in [0, pi/4]
  __m128d z = _mm_mul_pd(x, x);

  y = _mm_mul_pd(_mm_cst_pd_coscof_p0, z);
  y = _mm_add_pd(y, _mm_cst_pd_coscof_p1);
  y = _mm_mul_pd(y, z);
  y = _mm_add_pd(y, _mm_cst_pd_coscof_p2);
  y = _mm_mul_pd(y, z);
  y = _mm_add_pd(y, _mm_cst_pd_coscof_p3);
  y = _mm_mul_pd(y, z);
  y = _mm_add_pd(y, _mm_cst_pd_coscof_p4);
  y = _mm_mul_pd(y, z);
  y = _mm_add_pd(y, _mm_cst_pd_coscof_p5);
  y = _mm_mul_pd(y, z);
  y = _mm_mul_pd(y, z);
  __m128d tmp = _mm_mul_pd(z, _mm_cst_pd_0p5);
  y = _mm_sub_pd(y, tmp);
  y = _mm_add_pd(y, _mm_cst_pd_one);

  // Second polynom
  __m128d y2 = _mm_mul_pd(_mm_cst_pd_sincof_p0, z);
  y2 = _mm_add_pd(y2, _mm_cst_pd_sincof_p1);
  y2 = _mm_mul_pd(y2, z);
  y2 = _mm_add_pd(y2, _mm_cst_pd_sincof_p2);
  y2 = _mm_mul_pd(y2, z);
  y2 = _mm_add_pd(y2, _mm_cst_pd_sincof_p3);
  y2 = _mm_mul_pd(y2, z);
  y2 = _mm_add_pd(y2, _mm_cst_pd_sincof_p4);
  y2 = _mm_mul_pd(y2, z);
  y2 = _mm_add_pd(y2, _mm_cst_pd_sincof_p5);
  y2 = _mm_mul_pd(y2, z);
  y2 = _mm_mul_pd(y2, x);
  y2 = _mm_add_pd(y2, x);

  // Select the correct result
  y2 = _mm_and_pd(polymask, y2);
  y = _mm_andnot_pd(polymask, y);
  y = _mm_add_pd(y, y2);

  // Update the sign
  y = _mm_xor_pd(y, signbit);

  return y;
}

/*
 * Range reduction is into intervals of pi/4.  The reduction
 * error is nearly eliminated by contriving an extended precision
 * modular arithmetic.
 *
 * Two polynomial approximating functions are employed.
 * Between 0 and pi/4 the cosine is approximated by
 *      1  -  x**2 Q(x**2).
 * Between pi/4 and pi/2 the sine is represented as
 *      x  +  x**3 P(x**2).
 */
__m128 _mm_cos_ps(__m128 x)
{
  // Extract the absolute values
  x = _mm_and_ps(x, *(__m128*) _mm_cst_inv_sign_mask_ps);

  // Cephes method for scaling between 0 et pi/4
  __m128 y = _mm_mul_ps(x, _mm_cst_ps_fopi);
  __m128i yf = _mm_cvttps_epi32(y); // floor
  // see j = (j+1) & (~1) in Cephes
  yf = _mm_add_epi32(yf, *(__m128i*) _mm_cst_one);
  yf = _mm_and_si128(yf, *(__m128i*) _mm_cst_inv1);
  y = _mm_cvtepi32_ps(yf);

  yf = _mm_sub_epi32(yf, *(__m128i*) _mm_cst_two);

  // Swap sign flag
  __m128i flag = _mm_andnot_si128(yf, *(__m128i*) _mm_cst_four);
  flag = _mm_slli_epi32(flag, 29); // flag << 29

  // Polynom selection mask
  yf = _mm_and_si128(yf, *(__m128i*) _mm_cst_two);
  yf = _mm_cmpeq_epi32(yf, _mm_setzero_si128());

  __m128 signbit = _mm_castsi128_ps(flag);
  __m128 polymask = _mm_castsi128_ps(yf);

  // Extended precision modular arithmetic (Cody and Waite)
  // z = ((x - y * DP1) - y * DP2) - y * DP3
  __m128 xmm1 = _mm_mul_ps(y, _mm_cst_ps_DP1);
  __m128 xmm2 = _mm_mul_ps(y, _mm_cst_ps_DP2);
  __m128 xmm3 = _mm_mul_ps(y, _mm_cst_ps_DP3);
  x = _mm_sub_ps(x, xmm1);
  x = _mm_sub_ps(x, xmm2);
  x = _mm_sub_ps(x, xmm3);

  // First polynom x \in [0, pi/4]
  __m128 z = _mm_mul_ps(x, x);

  y = _mm_mul_ps(_mm_cst_ps_coscof_p0, z);
  y = _mm_add_ps(y, _mm_cst_ps_coscof_p1);
  y = _mm_mul_ps(y, z);
  y = _mm_add_ps(y, _mm_cst_ps_coscof_p2);
  y = _mm_mul_ps(y, z);
  y = _mm_mul_ps(y, z);
  __m128 tmp = _mm_mul_ps(z, _mm_cst_ps_0p5);
  y = _mm_sub_ps(y, tmp);
  y = _mm_add_ps(y, _mm_cst_ps_one);

  // Second polynom
  __m128 y2 = _mm_mul_ps(_mm_cst_ps_sincof_p0, z);
  y2 = _mm_add_ps(y2, _mm_cst_ps_sincof_p1);
  y2 = _mm_mul_ps(y2, z);
  y2 = _mm_add_ps(y2, _mm_cst_ps_sincof_p2);
  y2 = _mm_mul_ps(y2, z);
  y2 = _mm_mul_ps(y2, x);
  y2 = _mm_add_ps(y2, x);

  // Select the correct result
  y2 = _mm_and_ps(polymask, y2);
  y = _mm_andnot_ps(polymask, y);
  y = _mm_add_ps(y, y2);

  // Update the sign
  y = _mm_xor_ps(y, signbit);

  return y;
}

__m128d _mm_cos_pd(__m128d x)
{
  // Extract the absolute values
  x = _mm_and_pd(x, *(__m128d*) _mm_cst_inv_sign_mask_pd);

  // Cephes method for scaling between 0 et pi/4
  __m128d y = _mm_div_pd(x, _mm_cst_pd_pio4);
  __m64 yf = _mm_cvttpd_pi32(y); // floor
  // see j = (j+1) & (~1) in Cephes
  yf = _mm_add_pi32(yf, *(__m64*) _mm_cst_one);
  yf = _mm_and_si64(yf, *(__m64*) _mm_cst_inv1);
  y = _mm_cvtpi32_pd(yf);

  yf = _mm_sub_pi32(yf, *(__m64*) _mm_cst_two);

  // Swap sign flag
  __m64 flag = _mm_andnot_si64(yf, *(__m64*) _mm_cst_four);
  flag = _mm_slli_pi32(flag, 29); // flag << 29

  // Polynom selection mask
  yf = _mm_and_si64(yf, *(__m64*) _mm_cst_two);
  yf = _mm_cmpeq_pi32(yf, _mm_setzero_si64());

  static ALIGN16_BEG int signs[4] ALIGN16_END;
  signs[0] = 0;
  signs[1] = ((int*)(&flag))[0];
  signs[2] = 0;
  signs[3] = ((int*)(&flag))[1];

  static ALIGN16_BEG int poly[4] ALIGN16_END;
  poly[0] = ((int*)(&yf))[0];
  poly[1] = ((int*)(&yf))[0];
  poly[2] = ((int*)(&yf))[1];
  poly[3] = ((int*)(&yf))[1];

  __m128d signbit = _mm_load_pd((double*) signs);
  __m128d polymask = _mm_load_pd((double*) poly);

  // Extended precision modular arithmetic (Cody and Waite)
  // z = ((x - y * DP1) - y * DP2) - y * DP3
  __m128d xmm1 = _mm_mul_pd(y, _mm_cst_pd_DP1);
  __m128d xmm2 = _mm_mul_pd(y, _mm_cst_pd_DP2);
  __m128d xmm3 = _mm_mul_pd(y, _mm_cst_pd_DP3);
  x = _mm_sub_pd(x, xmm1);
  x = _mm_sub_pd(x, xmm2);
  x = _mm_sub_pd(x, xmm3);

  // First polynom x \in [0, pi/4]
  __m128d z = _mm_mul_pd(x, x);

  y = _mm_mul_pd(_mm_cst_pd_coscof_p0, z);
  y = _mm_add_pd(y, _mm_cst_pd_coscof_p1);
  y = _mm_mul_pd(y, z);
  y = _mm_add_pd(y, _mm_cst_pd_coscof_p2);
  y = _mm_mul_pd(y, z);
  y = _mm_add_pd(y, _mm_cst_pd_coscof_p3);
  y = _mm_mul_pd(y, z);
  y = _mm_add_pd(y, _mm_cst_pd_coscof_p4);
  y = _mm_mul_pd(y, z);
  y = _mm_add_pd(y, _mm_cst_pd_coscof_p5);
  y = _mm_mul_pd(y, z);
  y = _mm_mul_pd(y, z);
  __m128d tmp = _mm_mul_pd(z, _mm_cst_pd_0p5);
  y = _mm_sub_pd(y, tmp);
  y = _mm_add_pd(y, _mm_cst_pd_one);

  // Second polynom
  __m128d y2 = _mm_mul_pd(_mm_cst_pd_sincof_p0, z);
  y2 = _mm_add_pd(y2, _mm_cst_pd_sincof_p1);
  y2 = _mm_mul_pd(y2, z);
  y2 = _mm_add_pd(y2, _mm_cst_pd_sincof_p2);
  y2 = _mm_mul_pd(y2, z);
  y2 = _mm_add_pd(y2, _mm_cst_pd_sincof_p3);
  y2 = _mm_mul_pd(y2, z);
  y2 = _mm_add_pd(y2, _mm_cst_pd_sincof_p4);
  y2 = _mm_mul_pd(y2, z);
  y2 = _mm_add_pd(y2, _mm_cst_pd_sincof_p5);
  y2 = _mm_mul_pd(y2, z);
  y2 = _mm_mul_pd(y2, x);
  y2 = _mm_add_pd(y2, x);

  // Select the correct result
  y2 = _mm_and_pd(polymask, y2);
  y = _mm_andnot_pd(polymask, y);
  y = _mm_add_pd(y, y2);

  // Update the sign
  y = _mm_xor_pd(y, signbit);

  return y;
}

/*
 * Sine and cosine computations are very similar.
 * Both functions can be merged for computing sine and cosine at minimum
 * computation cost
 */
__m128 _mm_sincos_ps(__m128*c, __m128 x)
{
  // Extract the sign bit and work with absolute values
  __m128 sin_signbit = _mm_and_ps(x, *(__m128*) _mm_cst_sign_mask_ps);
  x = _mm_and_ps(x, *(__m128*) _mm_cst_inv_sign_mask_ps);

  // Cephes method for scaling between 0 et pi/4
  __m128 y = _mm_mul_ps(x, _mm_cst_ps_fopi);
  __m128i yf = _mm_cvttps_epi32(y); // floor
  // see j = (j+1) & (~1) in Cephes
  yf = _mm_add_epi32(yf, *(__m128i*) _mm_cst_one);
  yf = _mm_and_si128(yf, *(__m128i*) _mm_cst_inv1);
  y = _mm_cvtepi32_ps(yf);

  // Fork yf to yf_cos for cosine
  __m128i yf_cos = _mm_sub_epi32(yf, *(__m128i*) _mm_cst_two);
  yf_cos = _mm_andnot_si128(yf_cos, *(__m128i*) _mm_cst_four);
  yf_cos = _mm_slli_epi32(yf_cos, 29); // flag << 29
  __m128 cos_signbit = _mm_castsi128_ps(yf_cos);

  // Back to sine computation: swap sign flag
  __m128i flag = _mm_and_si128(yf, *(__m128i*) _mm_cst_four);
  flag = _mm_slli_epi32(flag, 29); // flag << 29

  // Polynom selection mask
  yf = _mm_and_si128(yf, *(__m128i*) _mm_cst_two);
  yf = _mm_cmpeq_epi32(yf, _mm_setzero_si128());

  __m128 swapsign = _mm_castsi128_ps(flag);
  __m128 polymask = _mm_castsi128_ps(yf);
  sin_signbit = _mm_xor_ps(sin_signbit, swapsign);

  // Extended precision modular arithmetic (Cody and Waite)
  // z = ((x - y * DP1) - y * DP2) - y * DP3
  __m128 xmm1 = _mm_mul_ps(y, _mm_cst_ps_DP1);
  __m128 xmm2 = _mm_mul_ps(y, _mm_cst_ps_DP2);
  __m128 xmm3 = _mm_mul_ps(y, _mm_cst_ps_DP3);
  x = _mm_sub_ps(x, xmm1);
  x = _mm_sub_ps(x, xmm2);
  x = _mm_sub_ps(x, xmm3);

  // First polynom x \in [0, pi/4]
  __m128 z = _mm_mul_ps(x, x);

  y = _mm_mul_ps(_mm_cst_ps_coscof_p0, z);
  y = _mm_add_ps(y, _mm_cst_ps_coscof_p1);
  y = _mm_mul_ps(y, z);
  y = _mm_add_ps(y, _mm_cst_ps_coscof_p2);
  y = _mm_mul_ps(y, z);
  y = _mm_mul_ps(y, z);
  __m128 tmp = _mm_mul_ps(z, _mm_cst_ps_0p5);
  y = _mm_sub_ps(y, tmp);
  y = _mm_add_ps(y, _mm_cst_ps_one);

  // Second polynom
  __m128 y2 = _mm_mul_ps(_mm_cst_ps_sincof_p0, z);
  y2 = _mm_add_ps(y2, _mm_cst_ps_sincof_p1);
  y2 = _mm_mul_ps(y2, z);
  y2 = _mm_add_ps(y2, _mm_cst_ps_sincof_p2);
  y2 = _mm_mul_ps(y2, z);
  y2 = _mm_mul_ps(y2, x);
  y2 = _mm_add_ps(y2, x);

  // Select the correct result
  __m128 ysin2 = _mm_and_ps(polymask, y2);
  __m128 ysin1 = _mm_andnot_ps(polymask, y);

  __m128 ycos2 = _mm_sub_ps(y2, ysin2);
  __m128 ycos1 = _mm_sub_ps(y, ysin1);

  __m128 ys = _mm_add_ps(ysin1, ysin2);
  __m128 yc = _mm_add_ps(ycos1, ycos2);

  // Update the sign
  *c = _mm_xor_ps(yc, cos_signbit);
  return _mm_xor_ps(ys, sin_signbit);

}

__m128d _mm_sincos_pd(__m128d*c, __m128d x)
{
  // Extract the sign bit and work with absolute values
  __m128d sin_signbit = _mm_and_pd(x, *(__m128d*) _mm_cst_sign_mask_pd);
  x = _mm_and_pd(x, *(__m128d*) _mm_cst_inv_sign_mask_pd);

  // Cephes method for scaling between 0 et pi/4
  __m128d y = _mm_div_pd(x, _mm_cst_pd_pio4);
  __m64 yf = _mm_cvttpd_pi32(y); // floor
  // see j = (j+1) & (~1) in Cephes
  yf = _mm_add_pi32(yf, *(__m64*) _mm_cst_one);
  yf = _mm_and_si64(yf, *(__m64*) _mm_cst_inv1);
  y = _mm_cvtpi32_pd(yf);

  // Fork yf to yf_cos for cosine
  __m64 yf_cos = _mm_sub_pi32(yf, *(__m64*) _mm_cst_two);
  yf_cos = _mm_andnot_si64(yf_cos, *(__m64*) _mm_cst_four);
  yf_cos = _mm_slli_pi32(yf_cos, 29); // flag << 29

  static ALIGN16_BEG int cosigns[4] ALIGN16_END;
  cosigns[0] = 0;
  cosigns[1] = ((int*)(&yf_cos))[0];
  cosigns[2] = 0;
  cosigns[3] = ((int*)(&yf_cos))[1];
  __m128d cos_signbit = _mm_load_pd((double*) cosigns);

  // Back to sine computation: swap sign flag
  __m64 flag = _mm_and_si64(yf, *(__m64*) _mm_cst_four);
  flag = _mm_slli_pi32(flag, 29); // flag << 29

  // Polynom selection mask
  yf = _mm_and_si64(yf, *(__m64*) _mm_cst_two);
  yf = _mm_cmpeq_pi32(yf, _mm_setzero_si64());

  static ALIGN16_BEG int signs[4] ALIGN16_END;
  signs[0] = 0;
  signs[1] = ((int*)(&flag))[0];
  signs[2] = 0;
  signs[3] = ((int*)(&flag))[1];

  static ALIGN16_BEG int poly[4] ALIGN16_END;
  poly[0] = ((int*)(&yf))[0];
  poly[1] = ((int*)(&yf))[0];
  poly[2] = ((int*)(&yf))[1];
  poly[3] = ((int*)(&yf))[1];

  __m128d swapsign = _mm_load_pd((double*) signs);
  __m128d polymask = _mm_load_pd((double*) poly);
  sin_signbit = _mm_xor_pd(sin_signbit, swapsign);

  // Extended precision modular arithmetic (Cody and Waite)
  // z = ((x - y * DP1) - y * DP2) - y * DP3
  __m128d xmm1 = _mm_mul_pd(y, _mm_cst_pd_DP1);
  __m128d xmm2 = _mm_mul_pd(y, _mm_cst_pd_DP2);
  __m128d xmm3 = _mm_mul_pd(y, _mm_cst_pd_DP3);
  x = _mm_sub_pd(x, xmm1);
  x = _mm_sub_pd(x, xmm2);
  x = _mm_sub_pd(x, xmm3);

  // First polynom x \in [0, pi/4]
  __m128d z = _mm_mul_pd(x, x);

  y = _mm_mul_pd(_mm_cst_pd_coscof_p0, z);
  y = _mm_add_pd(y, _mm_cst_pd_coscof_p1);
  y = _mm_mul_pd(y, z);
  y = _mm_add_pd(y, _mm_cst_pd_coscof_p2);
  y = _mm_mul_pd(y, z);
  y = _mm_add_pd(y, _mm_cst_pd_coscof_p3);
  y = _mm_mul_pd(y, z);
  y = _mm_add_pd(y, _mm_cst_pd_coscof_p4);
  y = _mm_mul_pd(y, z);
  y = _mm_add_pd(y, _mm_cst_pd_coscof_p5);
  y = _mm_mul_pd(y, z);
  y = _mm_mul_pd(y, z);
  __m128d tmp = _mm_mul_pd(z, _mm_cst_pd_0p5);
  y = _mm_sub_pd(y, tmp);
  y = _mm_add_pd(y, _mm_cst_pd_one);

  // Second polynom
  __m128d y2 = _mm_mul_pd(_mm_cst_pd_sincof_p0, z);
  y2 = _mm_add_pd(y2, _mm_cst_pd_sincof_p1);
  y2 = _mm_mul_pd(y2, z);
  y2 = _mm_add_pd(y2, _mm_cst_pd_sincof_p2);
  y2 = _mm_mul_pd(y2, z);
  y2 = _mm_add_pd(y2, _mm_cst_pd_sincof_p3);
  y2 = _mm_mul_pd(y2, z);
  y2 = _mm_add_pd(y2, _mm_cst_pd_sincof_p4);
  y2 = _mm_mul_pd(y2, z);
  y2 = _mm_add_pd(y2, _mm_cst_pd_sincof_p5);
  y2 = _mm_mul_pd(y2, z);
  y2 = _mm_mul_pd(y2, x);
  y2 = _mm_add_pd(y2, x);

  // Select the correct result
  __m128d ysin2 = _mm_and_pd(polymask, y2);
  __m128d ysin1 = _mm_andnot_pd(polymask, y);

  __m128d ycos2 = _mm_sub_pd(y2, ysin2);
  __m128d ycos1 = _mm_sub_pd(y, ysin1);

  __m128d ys = _mm_add_pd(ysin1, ysin2);
  __m128d yc = _mm_add_pd(ycos1, ycos2);

  // Update the sign
  *c = _mm_xor_pd(yc, cos_signbit);
  return _mm_xor_pd(ys, sin_signbit);

}


/* Returns radian angle between -pi/2 and +pi/2 whose cosine is x.
 *
 * Analytically, acos(x) = pi/2 - asin(x).
 *
 * However if |x| is near 1, there is cancellation error in subtracting
 * asin(x) from pi/2.  Hence if x < -0.5,
 *    acos(x) = pi - 2.0 * asin( sqrt((1+x)/2) );
 *
 * or if x > +0.5,
 *    acos(x) = 2.0 * asin(  sqrt((1-x)/2)
 */

__m128 _mm_asin_ps(__m128 x)
{
  __m128 signbit = _mm_and_ps(x, *(__m128*) _mm_cst_sign_mask_ps);
  __m128 nanmask = _mm_cmplt_ps(x, _mm_cst_ps_one);
  x = _mm_and_ps(x, *(__m128*) _mm_cst_inv_sign_mask_ps);


  // if a > 0.5
  __m128 z1 = _mm_mul_ps(_mm_cst_ps_0p5,
                         _mm_sub_ps(_mm_cst_ps_one, x));
  __m128 x1 = _mm_sqrt_ps (z1);

  // else
  __m128 z = _mm_mul_ps(x, x);

  __m128 mask = _mm_cmpgt_ps(x, _mm_cst_ps_0p5);
  z = _mm_andnot_ps(mask, z);
  z = _mm_add_ps(z, _mm_and_ps(mask, z1 ));

  x = _mm_andnot_ps(mask, x);
  x = _mm_add_ps(x, _mm_and_ps(mask, x1));

  __m128 tmp = _mm_mul_ps(z, _mm_cst_ps_asincof_p0);

  tmp = _mm_add_ps(tmp, _mm_cst_ps_asincof_p1);
  tmp = _mm_mul_ps(tmp, z);
  tmp = _mm_add_ps(tmp, _mm_cst_ps_asincof_p2);
  tmp = _mm_mul_ps(tmp, z);
  tmp = _mm_add_ps(tmp, _mm_cst_ps_asincof_p3);
  tmp = _mm_mul_ps(tmp, z);
  tmp = _mm_add_ps(tmp, _mm_cst_ps_asincof_p4);
  tmp = _mm_mul_ps(tmp, z);

  tmp = _mm_mul_ps(tmp, x);
  tmp = _mm_add_ps(tmp, x);

  __m128 tmp2 = _mm_add_ps(tmp, tmp);
  tmp2 = _mm_sub_ps(_mm_cst_ps_pio2, tmp2);

  tmp = _mm_andnot_ps(mask, tmp);
  tmp = _mm_add_ps(tmp, _mm_and_ps(mask, tmp2));

  tmp = _mm_xor_ps(tmp, signbit);

  // if x > 1 return 0
  tmp = _mm_and_ps(tmp, nanmask);

  return tmp;
}


__m128d _mm_asin_pd(__m128d x)
{
  // Extract the sign bit and work with absolute values
  __m128d signbit = _mm_and_pd(x, *(__m128d*) _mm_cst_sign_mask_pd);
  __m128d nanmask = _mm_cmplt_pd(x, _mm_cst_pd_one);

  x = _mm_and_pd(x, *(__m128d*) _mm_cst_inv_sign_mask_pd);

  __m128d mask = _mm_cmpgt_pd(x, _mm_cst_pd_0p5);

  // if x > 0.5
  /* arcsin(1-x) = pi/2 - sqrt(2x)(1+R(x))  */
  __m128d z1   = _mm_sub_pd(_mm_cst_pd_one, x);
  __m128d z   = _mm_mul_pd(x, x);
  z = _mm_add_pd(_mm_and_pd(mask, z1),
                 _mm_andnot_pd(mask, z));

  __m128d tmp2 = _mm_mul_pd(z, _mm_cst_pd_asincof_r0);

  tmp2 = _mm_add_pd(tmp2, _mm_cst_pd_asincof_r1);
  tmp2 = _mm_mul_pd(tmp2, z);
  tmp2 = _mm_add_pd(tmp2, _mm_cst_pd_asincof_r2);
  tmp2 = _mm_mul_pd(tmp2, z);
  tmp2 = _mm_add_pd(tmp2, _mm_cst_pd_asincof_r3);
  tmp2 = _mm_mul_pd(tmp2, z);
  tmp2 = _mm_add_pd(tmp2, _mm_cst_pd_asincof_r4);
  tmp2 = _mm_mul_pd(tmp2, z);

  __m128d den2 = _mm_add_pd(z, _mm_cst_pd_asincof_s0);
  den2 = _mm_mul_pd(den2, z);
  den2 = _mm_add_pd(den2, _mm_cst_pd_asincof_s1);
  den2 = _mm_mul_pd(den2, z);
  den2 = _mm_add_pd(den2, _mm_cst_pd_asincof_s2);
  den2 = _mm_mul_pd(den2, z);
  den2 = _mm_add_pd(den2, _mm_cst_pd_asincof_s3);

  tmp2 = _mm_div_pd(tmp2, den2);

  __m128d zz = _mm_sqrt_pd(_mm_add_pd(z, z));
  zz = _mm_sub_pd(_mm_sub_pd(_mm_cst_pd_pio2, zz),
                  _mm_mul_pd(zz, tmp2));
  tmp2 = _mm_and_pd(mask, zz);


  // else
  __m128d tmp = _mm_mul_pd(z, _mm_cst_pd_asincof_p0);

  tmp = _mm_add_pd(tmp, _mm_cst_pd_asincof_p1);
  tmp = _mm_mul_pd(tmp, z);
  tmp = _mm_add_pd(tmp, _mm_cst_pd_asincof_p2);
  tmp = _mm_mul_pd(tmp, z);
  tmp = _mm_add_pd(tmp, _mm_cst_pd_asincof_p3);
  tmp = _mm_mul_pd(tmp, z);
  tmp = _mm_add_pd(tmp, _mm_cst_pd_asincof_p4);
  tmp = _mm_mul_pd(tmp, z);
  tmp = _mm_add_pd(tmp, _mm_cst_pd_asincof_p5);
  tmp = _mm_mul_pd(tmp, z);

  __m128d den = _mm_add_pd(z, _mm_cst_pd_asincof_q0);
  den = _mm_mul_pd(den, z);
  den = _mm_add_pd(den, _mm_cst_pd_asincof_q1);
  den = _mm_mul_pd(den, z);
  den = _mm_add_pd(den, _mm_cst_pd_asincof_q2);
  den = _mm_mul_pd(den, z);
  den = _mm_add_pd(den, _mm_cst_pd_asincof_q3);
  den = _mm_mul_pd(den, z);
  den = _mm_add_pd(den, _mm_cst_pd_asincof_q4);

  tmp = _mm_div_pd(tmp, den);

  tmp = _mm_mul_pd(tmp, x);
  tmp = _mm_add_pd(tmp, x);
  tmp = _mm_andnot_pd(mask, tmp);

  tmp = _mm_add_pd(tmp2, tmp);

  // Update the sign
  tmp = _mm_xor_pd(tmp, signbit);

  // if x > 1 return 0
  tmp = _mm_and_pd(tmp, nanmask);

  return tmp;
}

/*
 * Returns radian angle between -pi/2 and +pi/2 whose tangent is x.
 *
 * Range reduction is from four intervals into the interval from zero to
 * tan(pi/8). A polynomial approximates the function in this basic interval.
 */
__m128 _mm_atan_ps(__m128 x)
{
  // Extract the sign bit and work with absolute values
  __m128 signbit = _mm_and_ps(x, *(__m128*) _mm_cst_sign_mask_ps);
  __m128 y       = _mm_setzero_ps();

  x = _mm_and_ps(x, *(__m128*) _mm_cst_inv_sign_mask_ps);

  // Range reduction
  __m128 x2 = _mm_div_ps(_mm_cst_ps_mone, x);
  __m128 x3 = _mm_div_ps(_mm_sub_ps (x, _mm_cst_ps_one),
                         _mm_add_ps (x, _mm_cst_ps_one));

  // if x > tan(3 pi/8) (= 2.4142...)
  __m128 mask = _mm_cmpgt_ps(x, _mm_cst_ps_tan3pio8);
  __m128 y2 = _mm_and_ps(mask, _mm_cst_ps_pio2);
  x2 = _mm_and_ps(mask, x2);
  y = _mm_andnot_ps(mask, y);
  x = _mm_andnot_ps(mask, x);
  y = _mm_add_ps(y, y2);
  x = _mm_add_ps(x, x2);

  // if x > tan(pi/8) (= 0.4142...)
  mask = _mm_cmpgt_ps(x, _mm_cst_ps_tanpio8);
  y2 = _mm_and_ps(mask, _mm_cst_ps_pio4);
  x3 = _mm_and_ps(mask, x3);
  y = _mm_andnot_ps(mask, y);
  x = _mm_andnot_ps(mask, x);
  y = _mm_add_ps(y, y2);
  x = _mm_add_ps(x, x3);


  // Polynom computation
  __m128 z   = _mm_mul_ps(x, x);
  __m128 tmp = _mm_mul_ps(z, _mm_cst_ps_atancof_p0);

  tmp = _mm_add_ps(tmp, _mm_cst_ps_atancof_p1);
  tmp = _mm_mul_ps(tmp, z);
  tmp = _mm_add_ps(tmp, _mm_cst_ps_atancof_p2);
  tmp = _mm_mul_ps(tmp, z);
  tmp = _mm_add_ps(tmp, _mm_cst_ps_atancof_p3);
  tmp = _mm_mul_ps(tmp, z);
  tmp = _mm_mul_ps(tmp, x);
  tmp = _mm_add_ps(tmp, x);

  y = _mm_add_ps(y, tmp);

  // Update the sign
  y = _mm_xor_ps(y, signbit);

  return y;
}

__m128d _mm_atan_pd(__m128d x)
{
  // Extract the sign bit and work with absolute values
  __m128d signbit = _mm_and_pd(x, *(__m128d*) _mm_cst_sign_mask_pd);
  __m128d y       = _mm_setzero_pd();

  x = _mm_and_pd(x, *(__m128d*) _mm_cst_inv_sign_mask_pd);

  // Range reduction
  __m128d x2 = _mm_div_pd(_mm_cst_pd_mone, x);
  __m128d x3 = _mm_div_pd(_mm_sub_pd (x, _mm_cst_pd_one),
                          _mm_add_pd (x, _mm_cst_pd_one));

  // if x > tan(3 pi/8) (= 2.4142...)
  __m128d mask = _mm_cmpgt_pd(x, _mm_cst_pd_tan3pio8);
  __m128d y2 = _mm_and_pd(mask, _mm_cst_pd_pio2);
  x2 = _mm_and_pd(mask, x2);
  y = _mm_andnot_pd(mask, y);
  x = _mm_andnot_pd(mask, x);
  y = _mm_add_pd(y, y2);
  x = _mm_add_pd(x, x2);

  // if x > tan(pi/8) (= 0.4142...)
  mask = _mm_cmpgt_pd(x, _mm_cst_pd_tanpio8);
  y2 = _mm_and_pd(mask, _mm_cst_pd_pio4);
  x3 = _mm_and_pd(mask, x3);
  y = _mm_andnot_pd(mask, y);
  x = _mm_andnot_pd(mask, x);
  y = _mm_add_pd(y, y2);
  x = _mm_add_pd(x, x3);


  // Polynom computation
  __m128d z   = _mm_mul_pd(x, x);
  __m128d tmp = _mm_mul_pd(z, _mm_cst_pd_atancof_p0);

  tmp = _mm_add_pd(tmp, _mm_cst_pd_atancof_p1);
  tmp = _mm_mul_pd(tmp, z);
  tmp = _mm_add_pd(tmp, _mm_cst_pd_atancof_p2);
  tmp = _mm_mul_pd(tmp, z);
  tmp = _mm_add_pd(tmp, _mm_cst_pd_atancof_p3);
  tmp = _mm_mul_pd(tmp, z);
  tmp = _mm_add_pd(tmp, _mm_cst_pd_atancof_p4);
  tmp = _mm_mul_pd(tmp, z);

  __m128d den = _mm_add_pd(z, _mm_cst_pd_atancof_q0);
  den = _mm_mul_pd(den, z);
  den = _mm_add_pd(den, _mm_cst_pd_atancof_q1);
  den = _mm_mul_pd(den, z);
  den = _mm_add_pd(den, _mm_cst_pd_atancof_q2);
  den = _mm_mul_pd(den, z);
  den = _mm_add_pd(den, _mm_cst_pd_atancof_q3);
  den = _mm_mul_pd(den, z);
  den = _mm_add_pd(den, _mm_cst_pd_atancof_q4);

  tmp = _mm_div_pd(tmp, den);

  tmp = _mm_mul_pd(tmp, x);
  tmp = _mm_add_pd(tmp, x);

  y = _mm_add_pd(y, tmp);

  // Update the sign
  y = _mm_xor_pd(y, signbit);

  return y;
}

/*
 * Returns radian angle whose tangent is y/x between -PI and PI
 */
__m128 _mm_atan2_ps(__m128 y, __m128 x)
{

  static __m128 zero = _mm_setzero_ps();

  __m128 mask = _mm_cmplt_ps(x, zero);
  __m128 w    = _mm_and_ps(mask, _mm_cst_ps_pi);
  mask        = _mm_cmplt_ps(y, zero);
  __m128 mone = _mm_and_ps(mask, _mm_cst_ps_mtwo);
  // {-2, 0} + 1 -> {-1, 1} -> -pi, 0, pi
  w           = _mm_mul_ps(_mm_add_ps(mone, _mm_cst_ps_one), w);

  __m128 q = _mm_div_ps(y, x);
  q = _mm_add_ps(w, _mm_atan_ps(q));

  __m128 mask2;

  // atan2(yneg, 0) -> -pi/2
  mask  = _mm_cmpeq_ps(x, zero);
  mask2 = _mm_and_ps(mask, _mm_cmplt_ps(y, zero));
  q     = _mm_andnot_ps(mask2, q);
  w     = _mm_and_ps(mask2, _mm_cst_ps_mpio2);
  q     = _mm_or_ps(w, q);

  // atan2(ypos, 0) -> pi/2
  mask2 = _mm_and_ps(mask, _mm_cmpgt_ps(y, zero));
  q     = _mm_andnot_ps(mask2, q);
  w     = _mm_and_ps(mask2, _mm_cst_ps_pio2);
  q     = _mm_or_ps(w, q);

  // atan2(0, 0) -> 0
  mask2 = _mm_and_ps(mask, _mm_cmpeq_ps(y, zero));
  q     = _mm_andnot_ps(mask2, q);
  w     = _mm_and_ps(mask2, zero);
  q     = _mm_or_ps(w, q);

  // atan2(0, xneg) -> pi
  mask  = _mm_cmplt_ps(x, zero);
  mask2 = _mm_and_ps(mask, _mm_cmpeq_ps(y, zero));
  q     = _mm_andnot_ps(mask2, q);
  w     = _mm_and_ps(mask2, _mm_cst_ps_pi);
  q     = _mm_or_ps(w, q);

  return q;
}

__m128d _mm_atan2_pd(__m128d y, __m128d x)
{

  static __m128d zero = _mm_setzero_pd();

  __m128d mask = _mm_cmplt_pd(x, zero);
  __m128d w    = _mm_and_pd(mask, _mm_cst_pd_pi);
  mask        = _mm_cmplt_pd(y, zero);
  __m128d mone = _mm_and_pd(mask, _mm_cst_pd_mtwo);
  // {-2, 0} + 1 -> {-1, 1} -> -pi, 0, pi
  w           = _mm_mul_pd(_mm_add_pd(mone, _mm_cst_pd_one), w);

  __m128d q = _mm_div_pd(y, x);
  q = _mm_add_pd(w, _mm_atan_pd(q));

  __m128d mask2;

  // atan2(yneg, 0) -> -pi/2
  mask  = _mm_cmpeq_pd(x, zero);
  mask2 = _mm_and_pd(mask, _mm_cmplt_pd(y, zero));
  q     = _mm_andnot_pd(mask2, q);
  w     = _mm_and_pd(mask2, _mm_cst_pd_mpio2);
  q     = _mm_or_pd(w, q);

  // atan2(ypos, 0) -> pi/2
  mask2 = _mm_and_pd(mask, _mm_cmpgt_pd(y, zero));
  q     = _mm_andnot_pd(mask2, q);
  w     = _mm_and_pd(mask2, _mm_cst_pd_pio2);
  q     = _mm_or_pd(w, q);

  // atan2(0, 0) -> 0
  mask2 = _mm_and_pd(mask, _mm_cmpeq_pd(y, zero));
  q     = _mm_andnot_pd(mask2, q);
  w     = _mm_and_pd(mask2, zero);
  q     = _mm_or_pd(w, q);

  // atan2(0, xneg) -> pi
  mask  = _mm_cmplt_pd(x, zero);
  mask2 = _mm_and_pd(mask, _mm_cmpeq_pd(y, zero));
  q     = _mm_andnot_pd(mask2, q);
  w     = _mm_and_pd(mask2, _mm_cst_pd_pi);
  q     = _mm_or_pd(w, q);

  return q;
}

#  endif
#endif
