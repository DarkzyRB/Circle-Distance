/*
 * @File: circle_distance.dart
 * @Project: CircleDistance
 * @Author: Paul Salcedo F.
 * @Date: Wednesday, 15th December 2021 09:41 pm
 * @Email:  paulsalcedo.dev@gmail.com
 * @Last modified by:   Paul Salcedo F.
 * @Last modified time: Wednesday, 15th December 2021 09:37 pm
 * @License: Privativo
 * @Copyright: Copyright (c) 2021 Paul Salcedo F.
 */

import 'dart:math';

// -----------------------------------------------------------------------
//
// ! Main class
// -----------------------------------------------------------------------
class CircleDistance {
  double latitude1;
  double longitude1;

  double latitude2;
  double longitude2;

  CircleDistance.fromRadians({
    required this.latitude1,
    required this.longitude1,
    required this.latitude2,
    required this.longitude2,
  }) {
    latitude1 = latitude1;
    longitude1 = longitude1;

    latitude2 = latitude2;
    longitude2 = longitude2;

    _throwExceptionOnInvalidCoordinates();
  }

  CircleDistance.fromDegrees({
    required this.latitude1,
    required this.longitude1,
    required this.latitude2,
    required this.longitude2,
  }) {
    latitude1 = _radiansFromDegrees(latitude1);
    longitude1 = _radiansFromDegrees(longitude1);

    latitude2 = _radiansFromDegrees(latitude2);
    longitude2 = _radiansFromDegrees(longitude2);

    _throwExceptionOnInvalidCoordinates();
  }

  double haversineDistance() {
    return Haversine.distance(latitude1, longitude1, latitude2, longitude2);
  }

  double sphericalLawOfCosinesDistance() {
    return SphericalLawOfCosines.distance(
      latitude1,
      longitude1,
      latitude2,
      longitude2,
    );
  }

  double vincentyDistance() {
    return Vincenty.distance(latitude1, longitude1, latitude2, longitude2);
  }

  double _radiansFromDegrees(final double degrees) => degrees * (pi / 180.0);

  bool _isValidCoordinate(double latitude, longitude) =>
      _isValidLatitude(latitude) && _isValidLongitude(longitude);

  bool _isValidLatitude(double latitudeInRadians) =>
      !(latitudeInRadians < _radiansFromDegrees(-90.0) ||
          latitudeInRadians > _radiansFromDegrees(90.0));

  bool _isValidLongitude(double longitudeInRadians) =>
      !(longitudeInRadians < _radiansFromDegrees(-180.0) ||
          longitudeInRadians > _radiansFromDegrees(180.0));

  void _throwExceptionOnInvalidCoordinates() {
    String invalidDescription = """
            A coordinate is considered invalid if it meets at least one of the following criteria:
            - Its latitude is greater than 90 degrees or less than -90 degrees.
            - Its longitude is greater than 180 degrees or less than -180 degrees.
            
            see https://en.wikipedia.org/wiki/Decimal_degrees
        """;

    if (!_isValidCoordinate(latitude1, longitude1)) {
      throw FormatException(
          "Invalid coordinates at latitude1|longitude1\n$invalidDescription");
    }

    if (!_isValidCoordinate(latitude2, longitude2)) {
      throw FormatException(
          "Invalid coordinates at latitude2|longitude2\n$invalidDescription");
    }
  }
}

// -----------------------------------------------------------------------
//
// ! Formulas
// -----------------------------------------------------------------------

class Haversine {
  static double distance(double lat1, lon1, lat2, lon2) {
    double earthRadius = 6378137.0; // WGS84 major axis
    double distance = 2 *
        earthRadius *
        asin(sqrt(pow(sin(lat2 - lat1) / 2, 2) +
            cos(lat1) * cos(lat2) * pow(sin(lon2 - lon1) / 2, 2)));
    return distance;
  }
}

class SphericalLawOfCosines {
  static double distance(double lat1, lon1, lat2, lon2) {
    double distance =
        acos(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon2 - lon1));
    if (distance < 0) distance = distance + pi;

    double earthRadius = 6378137.0; // WGS84 major axis
    return earthRadius * distance;
  }
}

class Vincenty {
  static double distance(double lat1, lon1, lat2, lon2) {
    int MAXITERS = 20;

    double a = 6378137.0; // WGS84 major axis
    double b = 6356752.3142; // WGS84 semi-major axis
    double f = (a - b) / a;
    double aSqMinusBSqOverBSq = (a * a - b * b) / (b * b);
    double L = lon2 - lon1;
    double A = 0.0;
    double U1 = atan((1.0 - f) * tan(lat1));
    double U2 = atan((1.0 - f) * tan(lat2));
    double cosU1 = cos(U1);
    double cosU2 = cos(U2);
    double sinU1 = sin(U1);
    double sinU2 = sin(U2);
    double cosU1cosU2 = cosU1 * cosU2;
    double sinU1sinU2 = sinU1 * sinU2;
    double sigma = 0.0;
    double deltaSigma = 0.0;
    double cosSqAlpha = 0.0;
    double cos2SM = 0.0;
    double cosSigma = 0.0;
    double sinSigma = 0.0;
    double cosLambda = 0.0;
    double sinLambda = 0.0;
    double lambda = L; // initial guess
    for (int iter = 0; iter < MAXITERS; iter++) {
      double lambdaOrig = lambda;
      cosLambda = cos(lambda);
      sinLambda = sin(lambda);
      double t1 = cosU2 * sinLambda;
      double t2 = cosU1 * sinU2 - sinU1 * cosU2 * cosLambda;
      double sinSqSigma = t1 * t1 + t2 * t2; // (14)
      sinSigma = sqrt(sinSqSigma);
      cosSigma = sinU1sinU2 + cosU1cosU2 * cosLambda; // (15)
      sigma = atan2(sinSigma, cosSigma); // (16)
      double sinAlpha =
          (sinSigma == 0) ? 0.0 : cosU1cosU2 * sinLambda / sinSigma; // (17)
      cosSqAlpha = 1.0 - sinAlpha * sinAlpha;
      cos2SM = (cosSqAlpha == 0)
          ? 0.0
          : cosSigma - 2.0 * sinU1sinU2 / cosSqAlpha; // (18)
      double uSquared = cosSqAlpha * aSqMinusBSqOverBSq; // defn
      A = 1 +
          (uSquared / 16384.0) * // (3)
              (4096.0 +
                  uSquared * (-768 + uSquared * (320.0 - 175.0 * uSquared)));
      double B = (uSquared / 1024.0) * // (4)
          (256.0 + uSquared * (-128.0 + uSquared * (74.0 - 47.0 * uSquared)));
      double C = (f / 16.0) *
          cosSqAlpha *
          (4.0 + f * (4.0 - 3.0 * cosSqAlpha)); // (10)
      double cos2SMSq = cos2SM * cos2SM;
      deltaSigma = B *
          sinSigma * // (6)
          (cos2SM +
              (B / 4.0) *
                  (cosSigma * (-1.0 + 2.0 * cos2SMSq) -
                      (B / 6.0) *
                          cos2SM *
                          (-3.0 + 4.0 * sinSigma * sinSigma) *
                          (-3.0 + 4.0 * cos2SMSq)));
      lambda = L +
          (1.0 - C) *
              f *
              sinAlpha *
              (sigma +
                  C *
                      sinSigma *
                      (cos2SM +
                          C *
                              cosSigma *
                              (-1.0 + 2.0 * cos2SM * cos2SM))); // (11)
      double delta = (lambda - lambdaOrig) / lambda;
      if (delta.abs() < 1.0e-12) {
        break;
      }
    }
    double distance = (b * A * (sigma - deltaSigma));
    return distance;
  }
}
