#include "sse2_math.h"
#include "logger.h"

#include <math.h>

#include <stdio.h>
#include <stdlib.h>

void test_sincos_f(Logger& log)
{
  static float PI = acosf(-1);

  float ALIGN16_BEG angles[256] ALIGN16_END;
  float ALIGN16_BEG sines[256] ALIGN16_END;
  float ALIGN16_BEG cosines[256] ALIGN16_END;
  float ALIGN16_BEG sincos_sin[256] ALIGN16_END;
  float ALIGN16_BEG sincos_cos[256] ALIGN16_END;

  __m128 mmsin, mmcos;

  for (size_t i = 0; i<256; ++i)
    angles[i] = 2.0f * PI * i / 256;

  for (size_t i = 0; i<64; ++i)
  {
    _mm_storeu_ps(sines + 4*i, _mm_sin_ps(_mm_loadu_ps(angles + 4*i)));
    _mm_storeu_ps(cosines + 4*i, _mm_cos_ps(_mm_loadu_ps(angles + 4*i)));
    mmsin = _mm_sincos_ps(&mmcos, _mm_loadu_ps(angles + 4*i));
    _mm_storeu_ps(sincos_sin + 4*i, mmsin);
    _mm_storeu_ps(sincos_cos + 4*i, mmcos);
  }

  for (size_t i = 0; i<256; ++i)
  {
    log.testfloat(__LINE__, cosines[i], cos(angles[i]), "cos_f");
    log.testfloat(__LINE__, sincos_cos[i], cos(angles[i]), "sincos_cos_f");
    log.testfloat(__LINE__, sines[i], sin(angles[i]), "sin_f");
    log.testfloat(__LINE__, sincos_sin[i], sin(angles[i]), "sincos_sin_f");
  }
}

void test_sincos_d(Logger& log)
{
  static double PI = acos(-1);

  double ALIGN16_BEG angles[256] ALIGN16_END;
  double ALIGN16_BEG sines[256] ALIGN16_END;
  double ALIGN16_BEG cosines[256] ALIGN16_END;
  double ALIGN16_BEG sincos_sin[256] ALIGN16_END;
  double ALIGN16_BEG sincos_cos[256] ALIGN16_END;

  __m128d mmsin, mmcos;

  for (size_t i = 0; i<256; ++i)
    angles[i] = 2.0 * PI * i / 256;

  for (size_t i = 0; i<128; ++i)
  {
    _mm_storeu_pd(sines + 2*i, _mm_sin_pd(_mm_loadu_pd(angles + 2*i)));
    _mm_storeu_pd(cosines + 2*i, _mm_cos_pd(_mm_loadu_pd(angles + 2*i)));
    mmsin = _mm_sincos_pd(&mmcos, _mm_loadu_pd(angles + 2*i));
    _mm_storeu_pd(sincos_sin + 2*i, mmsin);
    _mm_storeu_pd(sincos_cos + 2*i, mmcos);
  }

  for (size_t i = 0; i<256; ++i)
  {
    log.testdouble(__LINE__, cosines[i], cos(angles[i]), "cos_d");
    log.testdouble(__LINE__, sincos_cos[i], cos(angles[i]), "sincos_cos_d");
    log.testdouble(__LINE__, sines[i], sin(angles[i]), "sin_d");
    log.testdouble(__LINE__, sincos_sin[i], sin(angles[i]), "sincos_sin_d");
  }
}

void test_asin_f(Logger& log)
{
  float ALIGN16_BEG vals[256] ALIGN16_END;
  float ALIGN16_BEG asins[256] ALIGN16_END;

  for (size_t i = 0; i<256; ++i)
    vals[i] = (i - 128.f) / 128.f;
  for (size_t i = 0; i<64; ++i)
    _mm_storeu_ps(asins + 4*i, _mm_asin_ps(_mm_loadu_ps(vals + 4*i)));
  for (size_t i = 0; i<256; ++i)
    log.testfloat(__LINE__, asins[i], asin(vals[i]), "asin_f");
}

void test_asin_d(Logger& log)
{
  double ALIGN16_BEG vals[256] ALIGN16_END;
  double ALIGN16_BEG asins[256] ALIGN16_END;

  for (size_t i = 0; i<256; ++i)
    vals[i] = (i - 128.) / 128.;
  for (size_t i = 0; i<128; ++i)
    _mm_storeu_pd(asins + 2*i, _mm_asin_pd(_mm_loadu_pd(vals + 2*i)));
  for (size_t i = 0; i<256; ++i)
    log.testdouble(__LINE__, asins[i], asin(vals[i]), "asin_d");
}


void test_atan_f(Logger& log)
{
  float ALIGN16_BEG vals[256] ALIGN16_END;
  float ALIGN16_BEG atans[256] ALIGN16_END;

  for (size_t i = 0; i<256; ++i)
    vals[i] = (i - 128.f) / 10.f;
  for (size_t i = 0; i<64; ++i)
    _mm_storeu_ps(atans + 4*i, _mm_atan_ps(_mm_loadu_ps(vals + 4*i)));
  for (size_t i = 0; i<256; ++i)
    log.testfloat(__LINE__, atans[i], atan(vals[i]), "atan_f");
}


void test_atan_d(Logger& log)
{
  double ALIGN16_BEG vals[256] ALIGN16_END;
  double ALIGN16_BEG atans[256] ALIGN16_END;

  for (size_t i = 0; i<256; ++i)
    vals[i] = (i - 128.) / 10.;
  for (size_t i = 0; i<128; ++i)
    _mm_storeu_pd(atans + 2*i, _mm_atan_pd(_mm_loadu_pd(vals + 2*i)));
  for (size_t i = 0; i<256; ++i)
    log.testdouble(__LINE__, atans[i], atan(vals[i]), "atan_d");
}


void test_atan2_f(Logger& log)
{
  float ALIGN16_BEG valt[16] ALIGN16_END;
  float ALIGN16_BEG valu[16] ALIGN16_END;
  float ALIGN16_BEG vals[256] ALIGN16_END;
  float ALIGN16_BEG atan2s[256] ALIGN16_END;

  for (size_t i = 0; i<16; ++i)
  {
    valt[i] = (i - 8.f) / 20.f;
    valu[i] = (i - 8.f) / 20.f;
  }

  for (size_t i = 0; i<16; ++i)
    for (size_t j = 0; j<16; ++j)
      vals[i+15*j] = valt[i] / valu[j];

  for (size_t i = 0; i<4; ++i)
    for (size_t j = 0; j<16; ++j)
    {
      __m128 _mm_valu = _mm_set_ps(valu[j], valu[j], valu[j], valu[j]);
      _mm_storeu_ps(atan2s + 4*i+16*j,
                    _mm_atan2_ps(_mm_loadu_ps(valt + 4*i), _mm_valu));
    }

  for (size_t i = 0; i<16; ++i)
    for (size_t j = 0; j<16; ++j)
      log.testfloat(__LINE__, atan2s[i+16*j], atan2(valt[i], valu[j]),
                    "atan2_f", 2 * FLT_EPSILON);

}


void test_atan2_d(Logger& log)
{
  double ALIGN16_BEG valt[16] ALIGN16_END;
  double ALIGN16_BEG valu[16] ALIGN16_END;
  double ALIGN16_BEG vals[256] ALIGN16_END;
  double ALIGN16_BEG atan2s[256] ALIGN16_END;

  for (size_t i = 0; i<16; ++i)
  {
    valt[i] = (i - 8.) / 20.;
    valu[i] = (i - 8.) / 20.;
  }

  for (size_t i = 0; i<16; ++i)
    for (size_t j = 0; j<16; ++j)
      vals[i+15*j] = valt[i] / valu[j];

  for (size_t i = 0; i<8; ++i)
    for (size_t j = 0; j<16; ++j)
    {
      __m128d _mm_valu = _mm_set_pd(valu[j], valu[j]);
      _mm_storeu_pd(atan2s + 2*i+16*j,
                    _mm_atan2_pd(_mm_loadu_pd(valt + 2*i), _mm_valu));
    }

  for (size_t i = 0; i<16; ++i)
    for (size_t j = 0; j<16; ++j)
      log.testdouble(__LINE__, atan2s[i+16*j], atan2(valt[i], valu[j]),
                    "atan2_d", 2 * DBL_EPSILON);

}

int main() {
  Logger log(__FILE__);

  test_sincos_f(log);
  test_sincos_d(log);
  test_asin_f(log);
  test_asin_d(log);
  test_atan_f(log);
  test_atan_d(log);
  test_atan2_f(log);
  test_atan2_d(log);

  return log.reportexit();
}
