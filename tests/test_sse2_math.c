#include "sse2_math.h"
#include "logger.h"

#include <math.h>

#include <stdio.h>
#include <stdlib.h>

void test_sincos(Logger& log)
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
    log.testfloat(__LINE__, cosines[i], cos(angles[i]), "cos");
    log.testfloat(__LINE__, sincos_cos[i], cos(angles[i]), "sincos_cos");
    log.testfloat(__LINE__, sines[i], sin(angles[i]), "sin");
    log.testfloat(__LINE__, sincos_sin[i], sin(angles[i]), "sincos_sin");
  }
}

void test_atan(Logger& log)
{
  float ALIGN16_BEG vals[256] ALIGN16_END;
  float ALIGN16_BEG atans[256] ALIGN16_END;

  for (size_t i = 0; i<256; ++i)
    vals[i] = (i - 128.f) / 10.f;
  for (size_t i = 0; i<64; ++i)
    _mm_storeu_ps(atans + 4*i, _mm_atan_ps(_mm_loadu_ps(vals + 4*i)));
  for (size_t i = 0; i<256; ++i)
    log.testfloat(__LINE__, atans[i], atan(vals[i]), "atan");
}


void test_atan2(Logger& log)
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
                    "atan2", 2 * FLT_EPSILON);

}

int main() {
  Logger log(__FILE__);

  test_sincos(log);
  test_atan(log);
  test_atan2(log);

  return log.reportexit();
}
