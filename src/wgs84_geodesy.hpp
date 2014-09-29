#include <cmath>
#include <cfloat>
#include <iostream>
#include <limits>

#include "spherical_geodesy.h"

template<typename T>
void WGS84Geodesy<T>::distanceAndBearing(
    const T& fromlat, const T& fromlon, const T& tolat, const T& tolon,
    T& distance, T& bearing1, T& bearing2)
{

  const T beta1 = atan((1.0 - wgs84_f) * tan(fromlat));
  const T beta2 = atan((1.0 - wgs84_f) * tan(tolat));
  const T cosB1 = cos(beta1);
  const T cosB2 = cos(beta2);
  const T sinB1 = sin(beta1);
  const T sinB2 = sin(beta2);

  T sigma, cosSigma, sinSigma;
  T cosBn, cosSigmaM, cos2SigmaM, v, K3;
  T omega, cosOm, sinOm, olddeltaOm, deltaOm = 0.0;
  int loopMax = 10;

  // Vincenty method
  do
  {
    omega = tolon - fromlon + deltaOm;

    cosOm = cos(omega);
    sinOm = sin(omega);

    sigma = atan2(sqrt(pow(cosB2*sinOm, 2) +
                       pow((cosB1*sinB2 - sinB1*cosB2*cosOm), 2)),
                  sinB1*sinB2 + cosB1*cosB2*cosOm );

    cosSigma = cos(sigma);
    sinSigma = sin(sigma);

    // Undetermined
    if (std::abs(sinSigma) < std::numeric_limits<T>::epsilon())
      return SphericalGeodesy<T>::distanceAndBearing(
          fromlat, fromlon, tolat, tolon, distance, bearing1, bearing2);

    cosBn = cosB1 * cosB2 * sinOm / sinSigma;

    // Undetermined
    if ((1.0-cosBn*cosBn) < std::numeric_limits<double>::epsilon())
      cosSigmaM = 0;
    else
      cosSigmaM = cosSigma - (2.0*sinB1*sinB2)/(1.0-cosBn*cosBn);

    v  = 0.25 * wgs84_f * (1.0-cosBn*cosBn);

    K3 = v * (1 +  wgs84_f + wgs84_f*wgs84_f
              - v * (3.0 + 7.0*wgs84_f - 13.0*v));

    cos2SigmaM = 2.0 * cosSigmaM * cosSigmaM - 1.0;

    olddeltaOm = deltaOm;

    deltaOm = (1.0 - K3) * wgs84_f * cosBn *
      (sigma + K3 * sinSigma * (cosSigmaM + K3*cosSigma*cos2SigmaM));

    loopMax--;

  } while ( (fabs(deltaOm-olddeltaOm) > 1e-12) && (loopMax > 0) );

  T secondEccentricSq =
    (wgs84_a*wgs84_a - wgs84_b*wgs84_b )/ (wgs84_b * wgs84_b);

  T t = 0.25 * secondEccentricSq * (1. - cosBn*cosBn);

  T K1 = 1.0 + t * (1. - 0.25 * t * (3. - t * (5. - 11. * t)));
  T K2 = t * (1 - t * (2. - 0.125 * t * (37. - 94. * t)));

  T cos3SigmaM = cosSigmaM * (4.0*cosSigmaM*cosSigmaM - 3.);

  T dSigma = K2 * sinSigma *
    (cosSigmaM + 0.25 * K2 * (cosSigma * (2.*cosSigmaM*cosSigmaM - 1.)
                              + K2/6. * (1. + 2. * (2.*cosSigma*cosSigma-1.))
                              * cos3SigmaM ) );

  bearing1 = atan2(cosB2 * sinOm, cosB1 * sinB2 - sinB1 * cosB2 * cosOm);
  bearing2 = atan2(cosB1 * sinOm, cosB1 * sinB2 * cosOm - sinB1 * cosB2  ) ;

  distance =  K1 * wgs84_b * (sigma - dSigma);
}


template<typename T>
void WGS84Geodesy<T>::destination(
    const T& fromlat, const T& fromlon, const T& bearing, const T& dist,
    T& tolat, T& tolon, T& bearing2)
{
  T sin_alpha1 = sin(bearing);
  T cos_alpha1 = cos(bearing);

  T tan_u1 = (1 - wgs84_f) * tan(fromlat);
  T cos_u1 = 1 / sqrt(1 + tan_u1*tan_u1);
  T sin_u1 = tan_u1 * cos_u1;

  T sigma1 = atan2(tan_u1, cos_alpha1);
  T sin_alpha = cos_u1 * sin_alpha1;
  T cos_sqalpha = 1 - sin_alpha * sin_alpha;
  T usq = cos_sqalpha *
    (wgs84_a * wgs84_a - wgs84_b * wgs84_b) / (wgs84_b * wgs84_b);
  T a = 1 + usq/16384.*(4096. + usq*(-768. + usq*(320 - 175*usq)));
  T b = usq/1024. * (256. + usq*(-128. + usq*(74. - 47.*usq)));

  T sigma = dist / (wgs84_b - a);
  T sigma_p = 0; //2 * PI;
  T cos2_sigmam, sin_sigma, cos_sigma, delta_sigma;

  do {
    cos2_sigmam = cos(2. * sigma1 + sigma);
    sin_sigma = sin(sigma);
    cos_sigma = cos(sigma);
    delta_sigma =
        b*sin_sigma*(cos2_sigmam +
                     b/4. * (cos_sigma*(-1. + 2. * cos2_sigmam * cos2_sigmam) -
                           b/6. * cos2_sigmam * (-3. + 4.*sin_sigma*sin_sigma)*
                           (-3. + 4 * cos2_sigmam * cos2_sigmam)));
    sigma_p = sigma;
    sigma = dist / (wgs84_b*a) + delta_sigma;
  } while (fabs(sigma - sigma_p) > 1e-12);

  T tmp = sin_u1*sin_sigma - cos_u1*cos_sigma*cos_alpha1;

  tolat = atan2(sin_u1*cos_sigma + cos_u1*sin_sigma*cos_alpha1,
                (1 - wgs84_f)*sqrt(sin_alpha*sin_alpha + tmp*tmp));

  T lambda = atan2(sin_sigma*sin_alpha1,
                   cos_u1*cos_sigma - sin_u1*sin_sigma*cos_alpha1);

  T c = wgs84_f/16.*cos_sqalpha*(4. + wgs84_f*(4 - 3*cos_sqalpha));
  T l = lambda - (1 - c) * wgs84_f * sin_alpha *
      (sigma + c*sin_sigma*(cos2_sigmam+c*cos_sigma*
                            (-1+2*cos2_sigmam*cos2_sigmam)));

  tolon = (fromlon + l);

  bearing2 = atan2(sin_alpha, -tmp);

}

