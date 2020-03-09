/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */
/* Vincenty Direct and Inverse Solution of Geodesics on the Ellipsoid (c) Chris Veness 2002-2019  */
/*                                                                                   MIT Licence  */
/* www.movable-type.co.uk/scripts/latlong-ellipsoidal-vincenty.html                               */
/* www.movable-type.co.uk/scripts/geodesy-library.html#latlon-ellipsoidal-vincenty                */
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */

import LatLonEllipsoidal, { Dms } from './latlon-ellipsoidal.js';

const pi = Math.PI;
const epsilon = Number.EPSILON;


/**
 * Distances & bearings between points, and destination points given start points & initial bearings,
 * calculated on an ellipsoidal earth model using ‘direct and inverse solutions of geodesics on the
 * ellipsoid’ devised by Thaddeus Vincenty.
 *
 * From: T Vincenty, "Direct and Inverse Solutions of Geodesics on the Ellipsoid with application of
 * nested equations", Survey Review, vol XXIII no 176, 1975. www.ngs.noaa.gov/PUBS_LIB/inverse.pdf.
 *
 * @module latlon-ellipsoidal-vincenty
 */

/* LatLonEllipsoidal_Vincenty - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

/**
 * Extends LatLonEllipsoidal with methods for calculating distances and bearings between points, and
 * destination points given distances and initial bearings, accurate to within 0.5mm distance,
 * 0.000015″ bearing.
 *
 * By default, these calculations are made on a WGS-84 ellipsoid. For geodesic calculations on other
 * ellipsoids, monkey-patch the LatLon point by setting the datum of ‘this’ point to make it appear
 * as a LatLonEllipsoidal_Datum or LatLonEllipsoidal_ReferenceFrame point: e.g.
 *
 *     import LatLon, { Dms } from '../latlon-ellipsoidal-vincenty.js';
 *     import { datums }      from '../latlon-ellipsoidal-datum.js';
 *     const le = new LatLon(50.065716, -5.713824);  // in OSGB-36
 *     const jog = new LatLon(58.644399, -3.068521); // in OSGB-36
 *     le.datum = datums.OSGB36;     // source point determines ellipsoid to use
 *     const d = le.distanceTo(jog); // = 969982.014; 27.848m more than on WGS-84 ellipsoid
 *
 * @extends LatLonEllipsoidal
 */
class LatLonEllipsoidal_Vincenty extends LatLonEllipsoidal {

    /**
     * Returns the distance between ‘this’ point and destination point along a geodesic on the
     * surface of the ellipsoid, using Vincenty inverse solution.
     *
     * @param   {LatLon} point - Latitude/longitude of destination point.
     * @returns {number} Distance in metres between points or NaN if failed to converge.
     *
     * @example
     *   const p1 = new LatLon(50.06632, -5.71475);
     *   const p2 = new LatLon(58.64402, -3.07009);
     *   const d = p1.distanceTo(p2); // 969,954.166 m
     */
    distanceTo(point) {
        try {
            const dist = this.inverse(point).distance;
            return Number(dist.toFixed(3)); // round to 1mm precision
        } catch (e) {
            if (e instanceof EvalError) return NaN; // lambda > pi or failed to converge
            throw e;
        }
    }


    /**
     * Returns the initial bearing (forward azimuth) to travel along a geodesic from ‘this’ point to
     * the given point, using Vincenty inverse solution.
     *
     * @param   {LatLon} point - Latitude/longitude of destination point.
     * @returns {number} Initial bearing in degrees from north (0°..360°) or NaN if failed to converge.
     *
     * @example
     *   const p1 = new LatLon(50.06632, -5.71475);
     *   const p2 = new LatLon(58.64402, -3.07009);
     *   const b1 = p1.initialBearingTo(p2); // 9.1419°
     */
    initialBearingTo(point) {
        try {
            const brng = this.inverse(point).initialBearing;
            return Number(brng.toFixed(7)); // round to 0.001″ precision
        } catch (e) {
            if (e instanceof EvalError) return NaN; // lambda > pi or failed to converge
            throw e;
        }
    }


    /**
     * Returns the final bearing (reverse azimuth) having travelled along a geodesic from ‘this’
     * point to the given point, using Vincenty inverse solution.
     *
     * @param   {LatLon} point - Latitude/longitude of destination point.
     * @returns {number} Final bearing in degrees from north (0°..360°) or NaN if failed to converge.
     *
     * @example
     *   const p1 = new LatLon(50.06632, -5.71475);
     *   const p2 = new LatLon(58.64402, -3.07009);
     *   const b2 = p1.finalBearingTo(p2); // 11.2972°
     */
    finalBearingTo(point) {
        try {
            const brng = this.inverse(point).finalBearing;
            return Number(brng.toFixed(7)); // round to 0.001″ precision
        } catch (e) {
            if (e instanceof EvalError) return NaN; // lambda > pi or failed to converge
            throw e;
        }
    }


    /**
     * Returns the destination point having travelled the given distance along a geodesic given by
     * initial bearing from ‘this’ point, using Vincenty direct solution.
     *
     * @param   {number} distance - Distance travelled along the geodesic in metres.
     * @param   {number} initialBearing - Initial bearing in degrees from north.
     * @returns {LatLon} Destination point.
     *
     * @example
     *   const p1 = new LatLon(-37.95103, 144.42487);
     *   const p2 = p1.destinationPoint(54972.271, 306.86816); // 37.6528°S, 143.9265°E
     */
    destinationPoint(distance, initialBearing) {
        return this.direct(Number(distance), Number(initialBearing)).point;
    }


    /**
     * Returns the final bearing (reverse azimuth) having travelled along a geodesic given by initial
     * bearing for a given distance from ‘this’ point, using Vincenty direct solution.
     * TODO: arg order? (this is consistent with destinationPoint, but perhaps less intuitive)
     *
     * @param   {number} distance - Distance travelled along the geodesic in metres.
     * @param   {LatLon} initialBearing - Initial bearing in degrees from north.
     * @returns {number} Final bearing in degrees from north (0°..360°).
     *
     * @example
     *   const p1 = new LatLon(-37.95103, 144.42487);
     *   const b2 = p1.finalBearingOn(306.86816, 54972.271); // 307.1736°
     */
    finalBearingOn(distance, initialBearing) {
        const brng = this.direct(Number(distance), Number(initialBearing)).finalBearing;
        return Number(brng.toFixed(7)); // round to 0.001″ precision
    }


    /**
     * Vincenty direct calculation.
     *
     * Ellipsoid parameters are taken from datum of 'this' point. Height is ignored.
     *
     * @private
     * @param   {number} distance - Distance along bearing in metres.
     * @param   {number} initialBearing - Initial bearing in degrees from north.
     * @returns (Object} Object including point (destination point), finalBearing.
     * @throws  {RangeError} Point must be on surface of ellipsoid.
     * @throws  {EvalError}  Formula failed to converge.
     */
    direct(distance, initialBearing) {
        if (this.height != 0) throw new RangeError('point must be on the surface of the ellipsoid');

        const phi1 = this.lat.toRadians(), lambda1 = this.lon.toRadians();
        const alpha1 = initialBearing.toRadians();
        const s = distance;

        // allow alternative ellipsoid to be specified
        const ellipsoid = this.datum ? this.datum.ellipsoid : LatLonEllipsoidal.ellipsoids.WGS84;
        const { a, b, f } = ellipsoid;

        const sinalpha1 = Math.sin(alpha1);
        const cosalpha1 = Math.cos(alpha1);

        const tanU1 = (1-f) * Math.tan(phi1), cosU1 = 1 / Math.sqrt((1 + tanU1*tanU1)), sinU1 = tanU1 * cosU1;
        const sigma1 = Math.atan2(tanU1, cosalpha1); // sigma1 = angular distance on the sphere from the equator to P1
        const sinalpha = cosU1 * sinalpha1;          // alpha = azimuth of the geodesic at the equator
        const cosSqalpha = 1 - sinalpha*sinalpha;
        const uSq = cosSqalpha * (a*a - b*b) / (b*b);
        const A = 1 + uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq)));
        const B = uSq/1024 * (256+uSq*(-128+uSq*(74-47*uSq)));

        let sigma = s / (b*A), sinsigma = null, cossigma = null, deltasigma = null; // sigma = angular distance P₁ P₂ on the sphere
        let cos2sigmasquare = null; // sigmasquare = angular distance on the sphere from the equator to the midpoint of the line

        let sigmaprime = null, iterations = 0;
        do {
            cos2sigmasquare = Math.cos(2*sigma1 + sigma);
            sinsigma = Math.sin(sigma);
            cossigma = Math.cos(sigma);
            deltasigma = B*sinsigma*(cos2sigmasquare+B/4*(cossigma*(-1+2*cos2sigmasquare*cos2sigmasquare)-
                B/6*cos2sigmasquare*(-3+4*sinsigma*sinsigma)*(-3+4*cos2sigmasquare*cos2sigmasquare)));
            sigmaprime = sigma;
            sigma = s / (b*A) + deltasigma;
        } while (Math.abs(sigma-sigmaprime) > 1e-12 && ++iterations<100);
        if (iterations >= 100) throw new EvalError('Vincenty formula failed to converge'); // not possible?

        const x = sinU1*sinsigma - cosU1*cossigma*cosalpha1;
        const phi2 = Math.atan2(sinU1*cossigma + cosU1*sinsigma*cosalpha1, (1-f)*Math.sqrt(sinalpha*sinalpha + x*x));
        const lambda = Math.atan2(sinsigma*sinalpha1, cosU1*cossigma - sinU1*sinsigma*cosalpha1);
        const C = f/16*cosSqalpha*(4+f*(4-3*cosSqalpha));
        const L = lambda - (1-C) * f * sinalpha * (sigma + C*sinsigma*(cos2sigmasquare+C*cossigma*(-1+2*cos2sigmasquare*cos2sigmasquare)));
        const lambda2 = lambda1 + L;

        const alpha2 = Math.atan2(sinalpha, -x);

        const destinationPoint = new LatLonEllipsoidal_Vincenty(phi2.toDegrees(), lambda2.toDegrees(), 0, this.datum);

        return {
            point:        destinationPoint,
            finalBearing: Dms.wrap360(alpha2.toDegrees()),
            iterations:   iterations,
        };
    }


    /**
     * Vincenty inverse calculation.
     *
     * Ellipsoid parameters are taken from datum of 'this' point. Height is ignored.
     *
     * @private
     * @param   {LatLon} point - Latitude/longitude of destination point.
     * @returns {Object} Object including distance, initialBearing, finalBearing.
     * @throws  {TypeError}  Invalid point.
     * @throws  {RangeError} Points must be on surface of ellipsoid.
     * @throws  {EvalError}  Formula failed to converge.
     */
    inverse(point) {
        if (!(point instanceof LatLonEllipsoidal)) throw new TypeError(`invalid point ‘${point}’`);
        if (this.height!=0 || point.height!=0) throw new RangeError('point must be on the surface of the ellipsoid');

        const p1 = this, p2 = point;
        const phi1 = p1.lat.toRadians(), lambda1 = p1.lon.toRadians();
        const phi2 = p2.lat.toRadians(), lambda2 = p2.lon.toRadians();

        // allow alternative ellipsoid to be specified
        const ellipsoid = this.datum ? this.datum.ellipsoid : LatLonEllipsoidal.ellipsoids.WGS84;
        const { a, b, f } = ellipsoid;

        const L = lambda2 - lambda1; // L = difference in longitude, U = reduced latitude, defined by tan U = (1-f)·tanphi.
        const tanU1 = (1-f) * Math.tan(phi1), cosU1 = 1 / Math.sqrt((1 + tanU1*tanU1)), sinU1 = tanU1 * cosU1;
        const tanU2 = (1-f) * Math.tan(phi2), cosU2 = 1 / Math.sqrt((1 + tanU2*tanU2)), sinU2 = tanU2 * cosU2;

        const antipodal = Math.abs(L) > pi/2 || Math.abs(phi2-phi1) > pi/2;

        let lambda = L, sinlambda = null, coslambda = null; // lambda = difference in longitude on an auxiliary sphere
        let sigma = antipodal ? pi : 0, sinsigma = 0, cossigma = antipodal ? -1 : 1, sinSqsigma = null; // sigma = angular distance P₁ P₂ on the sphere
        let cos2sigmasquare = 1;                      // sigmasquare = angular distance on the sphere from the equator to the midpoint of the line
        let sinalpha = null, cosSqalpha = 1;         // alpha = azimuth of the geodesic at the equator
        let C = null;

        let lambdaprime = null, iterations = 0;
        do {
            sinlambda = Math.sin(lambda);
            coslambda = Math.cos(lambda);
            sinSqsigma = (cosU2*sinlambda) * (cosU2*sinlambda) + (cosU1*sinU2-sinU1*cosU2*coslambda) * (cosU1*sinU2-sinU1*cosU2*coslambda);
            if (Math.abs(sinSqsigma) < epsilon) break;  // co-incident/antipodal points (falls back on lambda/sigma = L)
            sinsigma = Math.sqrt(sinSqsigma);
            cossigma = sinU1*sinU2 + cosU1*cosU2*coslambda;
            sigma = Math.atan2(sinsigma, cossigma);
            sinalpha = cosU1 * cosU2 * sinlambda / sinsigma;
            cosSqalpha = 1 - sinalpha*sinalpha;
            cos2sigmasquare = (cosSqalpha != 0) ? (cossigma - 2*sinU1*sinU2/cosSqalpha) : 0; // on equatorial line cos²alpha = 0 (§6)
            C = f/16*cosSqalpha*(4+f*(4-3*cosSqalpha));
            lambdaprime = lambda;
            lambda = L + (1-C) * f * sinalpha * (sigma + C*sinsigma*(cos2sigmasquare+C*cossigma*(-1+2*cos2sigmasquare*cos2sigmasquare)));
            const iterationCheck = antipodal ? Math.abs(lambda)-pi : Math.abs(lambda);
            if (iterationCheck > pi) throw new EvalError('lambda > pi');
        } while (Math.abs(lambda-lambdaprime) > 1e-12 && ++iterations<1000);
        if (iterations >= 1000) throw new EvalError('Vincenty formula failed to converge');

        const uSq = cosSqalpha * (a*a - b*b) / (b*b);
        const A = 1 + uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq)));
        const B = uSq/1024 * (256+uSq*(-128+uSq*(74-47*uSq)));
        const deltasigma = B*sinsigma*(cos2sigmasquare+B/4*(cossigma*(-1+2*cos2sigmasquare*cos2sigmasquare)-
            B/6*cos2sigmasquare*(-3+4*sinsigma*sinsigma)*(-3+4*cos2sigmasquare*cos2sigmasquare)));

        const s = b*A*(sigma-deltasigma); // s = length of the geodesic

        // note special handling of exactly antipodal points where sin²sigma = 0 (due to discontinuity
        // atan2(0, 0) = 0 but atan2(epsilon, 0) = pi/2 / 90°) - in which case bearing is always meridional,
        // due north (or due south!)
        // alpha = azimuths of the geodesic; alpha2 the direction P₁ P₂ produced
        const alpha1 = Math.abs(sinSqsigma) < epsilon ? 0 : Math.atan2(cosU2*sinlambda,  cosU1*sinU2-sinU1*cosU2*coslambda);
        const alpha2 = Math.abs(sinSqsigma) < epsilon ? pi : Math.atan2(cosU1*sinlambda, -sinU1*cosU2+cosU1*sinU2*coslambda);

        return {
            distance:       s,
            initialBearing: Math.abs(s) < epsilon ? NaN : Dms.wrap360(alpha1.toDegrees()),
            finalBearing:   Math.abs(s) < epsilon ? NaN : Dms.wrap360(alpha2.toDegrees()),
            iterations:     iterations,
        };
    }

}


/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */

export { LatLonEllipsoidal_Vincenty as default, Dms };
