/*
 * Freely inspired (understand plagiarism) from
 * http://gruntthepeon.free.fr/ssemath/
 *
 * Filtered SSE out, renamed to please my eyes and added more explicit comments
 * from various other sources.
 *
 * atan2 and atan made by myself, inspiring myself from Cephes library
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

/*
 * TODO Find how to make this cleaner
 */

#  define M128I_CONST(N, V)\
  static const ALIGN16_BEG int _mm_cst_pi32_##N[4] ALIGN16_END = { V, V, V, V }

#  define M128T_CONST(N, T, V)\
  static const ALIGN16_BEG T _mm_cst_##N[4] ALIGN16_END = { V, V, V, V }

  M128T_CONST(sign_mask, int, 0x80000000);
  M128T_CONST(inv_sign_mask, int, ~0x80000000);

  M128I_CONST(one,   1);
  M128I_CONST(inv1, ~1);
  M128I_CONST(two,   2);
  M128I_CONST(four,  4);

  static const __m128 _mm_cst_ps_one   = _mm_set1_ps(1.0f);
  static const __m128 _mm_cst_ps_mone  = _mm_set1_ps(-1.0f);
  static const __m128 _mm_cst_ps_mtwo  = _mm_set1_ps(-2.0f);
  static const __m128 _mm_cst_ps_0p5   = _mm_set1_ps(0.5f);

  static const __m128 _mm_cst_ps_tan3pio8 = _mm_set1_ps(2.414213562373095f);
  static const __m128 _mm_cst_ps_tanpio8  = _mm_set1_ps(0.4142135623730950f);

  static const __m128 _mm_cst_ps_fopi     = _mm_set1_ps(1.27323954473516f);

  static const __m128 _mm_cst_ps_pi    = _mm_set1_ps(3.14159265358979f);
  static const __m128 _mm_cst_ps_mpi   = _mm_set1_ps(-3.14159265358979f);
  static const __m128 _mm_cst_ps_pio2  = _mm_set1_ps(1.5707963267948966f);
  static const __m128 _mm_cst_ps_mpio2 = _mm_set1_ps(-1.5707963267948966f);
  static const __m128 _mm_cst_ps_pio4  = _mm_set1_ps(0.7853981633974483f);

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
  static const __m128 _mm_cst_ps_DP1 = _mm_set1_ps((float) -0.78515625);
  static const __m128 _mm_cst_ps_DP2 = _mm_set1_ps((float) -2.4187564849853515625e-4);
  static const __m128 _mm_cst_ps_DP3 = _mm_set1_ps((float) -3.77489497744594108e-8);

  static const __m128 _mm_cst_ps_sincof_p0 = _mm_set1_ps((float) -1.9515295891e-4);
  static const __m128 _mm_cst_ps_sincof_p1 = _mm_set1_ps((float)  8.3321608736e-3);
  static const __m128 _mm_cst_ps_sincof_p2 = _mm_set1_ps((float) -1.6666654611e-1);

  static const __m128 _mm_cst_ps_coscof_p0 = _mm_set1_ps((float)  2.443315711809948e-5);
  static const __m128 _mm_cst_ps_coscof_p1 = _mm_set1_ps((float) -1.388731625493765e-3);
  static const __m128 _mm_cst_ps_coscof_p2 = _mm_set1_ps((float)  4.166664568298827e-2);

  static const __m128 _mm_cst_ps_atancof_p0 = _mm_set1_ps((float)  8.05374449538e-2);
  static const __m128 _mm_cst_ps_atancof_p1 = _mm_set1_ps((float) -1.38776856032e-1);
  static const __m128 _mm_cst_ps_atancof_p2 = _mm_set1_ps((float)  1.99777106478e-1);
  static const __m128 _mm_cst_ps_atancof_p3 = _mm_set1_ps((float) -3.33329491539e-1);

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
    __m128 signbit = _mm_and_ps(x, *(__m128*) _mm_cst_sign_mask);
    x = _mm_and_ps(x, *(__m128*) _mm_cst_inv_sign_mask);

    // Cephes method for scaling between 0 et pi/4
    __m128 y = _mm_mul_ps(x, _mm_cst_ps_fopi);
    __m128i yf = _mm_cvttps_epi32(y); // floor
    // see j = (j+1) & (~1) in Cephes
    yf = _mm_add_epi32(yf, *(__m128i*) _mm_cst_pi32_one);
    yf = _mm_and_si128(yf, *(__m128i*) _mm_cst_pi32_inv1);
    y = _mm_cvtepi32_ps(yf);

    // Swap sign flag
    __m128i flag = _mm_and_si128(yf, *(__m128i*) _mm_cst_pi32_four);
    flag = _mm_slli_epi32(flag, 29); // flag << 29

    // Polynom selection mask
    yf = _mm_and_si128(yf, *(__m128i*) _mm_cst_pi32_two);
    yf = _mm_cmpeq_epi32(yf, _mm_setzero_si128());

    __m128 swapsign = _mm_castsi128_ps(flag);
    __m128 polymask = _mm_castsi128_ps(yf);
    signbit = _mm_xor_ps(signbit, swapsign);

    // Extended precision modular arithmetic (Cody and Waite)
    // z = ((x - y * DP1) - y * DP2) - y * DP3
    __m128 xmm1 = _mm_mul_ps(y, _mm_cst_ps_DP1);
    __m128 xmm2 = _mm_mul_ps(y, _mm_cst_ps_DP2);
    __m128 xmm3 = _mm_mul_ps(y, _mm_cst_ps_DP3);
    x = _mm_add_ps(x, xmm1);
    x = _mm_add_ps(x, xmm2);
    x = _mm_add_ps(x, xmm3);

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
    x = _mm_and_ps(x, *(__m128*) _mm_cst_inv_sign_mask);

    // Cephes method for scaling between 0 et pi/4
    __m128 y = _mm_mul_ps(x, _mm_cst_ps_fopi);
    __m128i yf = _mm_cvttps_epi32(y); // floor
    // see j = (j+1) & (~1) in Cephes
    yf = _mm_add_epi32(yf, *(__m128i*) _mm_cst_pi32_one);
    yf = _mm_and_si128(yf, *(__m128i*) _mm_cst_pi32_inv1);
    y = _mm_cvtepi32_ps(yf);

    yf = _mm_sub_epi32(yf, *(__m128i*) _mm_cst_pi32_two);

    // Swap sign flag
    __m128i flag = _mm_andnot_si128(yf, *(__m128i*) _mm_cst_pi32_four);
    flag = _mm_slli_epi32(flag, 29); // flag << 29

    // Polynom selection mask
    yf = _mm_and_si128(yf, *(__m128i*) _mm_cst_pi32_two);
    yf = _mm_cmpeq_epi32(yf, _mm_setzero_si128());

    __m128 signbit = _mm_castsi128_ps(flag);
    __m128 polymask = _mm_castsi128_ps(yf);

    // Extended precision modular arithmetic (Cody and Waite)
    // z = ((x - y * DP1) - y * DP2) - y * DP3
    __m128 xmm1 = _mm_mul_ps(y, _mm_cst_ps_DP1);
    __m128 xmm2 = _mm_mul_ps(y, _mm_cst_ps_DP2);
    __m128 xmm3 = _mm_mul_ps(y, _mm_cst_ps_DP3);
    x = _mm_add_ps(x, xmm1);
    x = _mm_add_ps(x, xmm2);
    x = _mm_add_ps(x, xmm3);

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

  /*
   * Sine and cosine computations are very similar.
   * Both functions can be merged for computing sine and cosine at minimum
   * computation cost
   */
  __m128 _mm_sincos_ps(__m128*c, __m128 x)
  {
    // Extract the sign bit and work with absolute values
    __m128 sin_signbit = _mm_and_ps(x, *(__m128*) _mm_cst_sign_mask);
    x = _mm_and_ps(x, *(__m128*) _mm_cst_inv_sign_mask);

    // Cephes method for scaling between 0 et pi/4
    __m128 y = _mm_mul_ps(x, _mm_cst_ps_fopi);
    __m128i yf = _mm_cvttps_epi32(y); // floor
    // see j = (j+1) & (~1) in Cephes
    yf = _mm_add_epi32(yf, *(__m128i*) _mm_cst_pi32_one);
    yf = _mm_and_si128(yf, *(__m128i*) _mm_cst_pi32_inv1);
    y = _mm_cvtepi32_ps(yf);

    // Fork yf to yf_cos for cosine
    __m128i yf_cos = _mm_sub_epi32(yf, *(__m128i*) _mm_cst_pi32_two);
    yf_cos = _mm_andnot_si128(yf_cos, *(__m128i*) _mm_cst_pi32_four);
    yf_cos = _mm_slli_epi32(yf_cos, 29); // flag << 29
    __m128 cos_signbit = _mm_castsi128_ps(yf_cos);

    // Back to sine computation: swap sign flag
    __m128i flag = _mm_and_si128(yf, *(__m128i*) _mm_cst_pi32_four);
    flag = _mm_slli_epi32(flag, 29); // flag << 29

    // Polynom selection mask
    yf = _mm_and_si128(yf, *(__m128i*) _mm_cst_pi32_two);
    yf = _mm_cmpeq_epi32(yf, _mm_setzero_si128());

    __m128 swapsign = _mm_castsi128_ps(flag);
    __m128 polymask = _mm_castsi128_ps(yf);
    sin_signbit = _mm_xor_ps(sin_signbit, swapsign);

    // Extended precision modular arithmetic (Cody and Waite)
    // z = ((x - y * DP1) - y * DP2) - y * DP3
    __m128 xmm1 = _mm_mul_ps(y, _mm_cst_ps_DP1);
    __m128 xmm2 = _mm_mul_ps(y, _mm_cst_ps_DP2);
    __m128 xmm3 = _mm_mul_ps(y, _mm_cst_ps_DP3);
    x = _mm_add_ps(x, xmm1);
    x = _mm_add_ps(x, xmm2);
    x = _mm_add_ps(x, xmm3);

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

  /*
   * Returns radian angle between -pi/2 and +pi/2 whose tangent is x.
   *
   * Range reduction is from four intervals into the interval from zero to
   * tan(pi/8). A polynomial approximates the function in this basic interval.
   */
  __m128 _mm_atan_ps(__m128 x)
  {
    // Extract the sign bit and work with absolute values
    __m128 signbit = _mm_and_ps(x, *(__m128*) _mm_cst_sign_mask);
    __m128 y       = _mm_setzero_ps();

    x = _mm_and_ps(x, *(__m128*) _mm_cst_inv_sign_mask);

    // Range reduction
    __m128 x2 = _mm_div_ps( _mm_cst_ps_mone, x);
    __m128 x3 = _mm_div_ps( _mm_sub_ps (x, _mm_cst_ps_one),
                            _mm_add_ps (x, _mm_cst_ps_one) );

    // if x > tan(3 pi/8) (= 2.4142...)
    __m128 mask = _mm_cmpgt_ps( x, _mm_cst_ps_tan3pio8);
    __m128 y2 = _mm_and_ps(mask, _mm_cst_ps_pio2);
    x2 = _mm_and_ps(mask, x2);
    y = _mm_andnot_ps(mask, y);
    x = _mm_andnot_ps(mask, x);
    y = _mm_add_ps(y, y2);
    x = _mm_add_ps(x, x2);

    // if x > tan(pi/8) (= 0.4142...)
    mask = _mm_cmpgt_ps( x, _mm_cst_ps_tanpio8);
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

// #include <stdio.h>
// TODO TO BE TESTED
__m128 _mm_asin_ps(__m128 x)
{
  __m128 cosval = _mm_sub_ps(_mm_cst_ps_one, _mm_mul_ps(x,x));
//   cosval = _mm_mul_ps(cosval, _mm_rsqrt_ps(cosval));
  cosval = _mm_sqrt_ps(cosval);
//   ALIGN16_BEG float m[4] ALIGN16_END;
//   _mm_storeu_ps(m,_mm_atan2_ps(x, cosval) );
//   printf("asin %f %f %f %f\n", m[0], m[1], m[2], m[3]);
  return _mm_atan2_ps(x, cosval);
}

/*
inline __m128 _mm_acos_ps(__m128 x)
{
  __m128 tmp = _mm_atan_ps(
      _mm_rsqrt_ps(_mm_mul_ps(_mm_add_ps(_mm_cst_ps_one,x),
                              _mm_rcp_ps(_mm_sub_ps(_mm_cst_ps_one,x)))));
  return _mm_add_ps(tmp,tmp);
}
*/
#  endif
#endif
