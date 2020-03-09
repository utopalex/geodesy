/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */
/* Latitude/longitude spherical geodesy tools                         (c) Chris Veness 2002-2019  */
/*                                                                                   MIT Licence  */
/* www.movable-type.co.uk/scripts/latlong.html                                                    */
/* www.movable-type.co.uk/scripts/geodesy-library.html#latlon-spherical                           */
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */

import Dms from './dms.js';

const pi = Math.PI;


/**
 * Library of geodesy functions for operations on a spherical earth model.
 *
 * Includes distances, bearings, destinations, etc, for both great circle paths and rhumb lines,
 * and other related functions.
 *
 * All calculations are done using simple spherical trigonometric formulae.
 *
 * @module latlon-spherical
 */

// note greek letters (e.g. phi, lambda, theta) are used for angles in radians to distinguish from angles in
// degrees (e.g. lat, lon, brng)


/* LatLonSpherical - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */


/**
 * Latitude/longitude points on a spherical model earth, and methods for calculating distances,
 * bearings, destinations, etc on (orthodromic) great-circle paths and (loxodromic) rhumb lines.
 */
class LatLonSpherical {

    /**
     * Creates a latitude/longitude point on the earth’s surface, using a spherical model earth.
     *
     * @param  {number} lat - Latitude (in degrees).
     * @param  {number} lon - Longitude (in degrees).
     * @throws {TypeError} Invalid lat/lon.
     *
     * @example
     *   import LatLon from '/js/geodesy/latlon-spherical.js';
     *   const p = new LatLon(52.205, 0.119);
     */
    constructor(lat, lon) {
        if (isNaN(lat)) throw new TypeError(`invalid lat ‘${lat}’`);
        if (isNaN(lon)) throw new TypeError(`invalid lon ‘${lon}’`);

        this._lat = Dms.wrap90(Number(lat));
        this._lon = Dms.wrap180(Number(lon));
    }


    /**
     * Latitude in degrees north from equator (including aliases lat, latitude): can be set as
     * numeric or hexagesimal (deg-min-sec); returned as numeric.
     */
    get lat()       { return this._lat; }
    get latitude()  { return this._lat; }
    set lat(lat) {
        this._lat = isNaN(lat) ? Dms.wrap90(Dms.parse(lat)) : Dms.wrap90(Number(lat));
        if (isNaN(this._lat)) throw new TypeError(`invalid lat ‘${lat}’`);
    }
    set latitude(lat) {
        this._lat = isNaN(lat) ? Dms.wrap90(Dms.parse(lat)) : Dms.wrap90(Number(lat));
        if (isNaN(this._lat)) throw new TypeError(`invalid latitude ‘${lat}’`);
    }

    /**
     * Longitude in degrees east from international reference meridian (including aliases lon, lng,
     * longitude): can be set as numeric or hexagesimal (deg-min-sec); returned as numeric.
     */
    get lon()       { return this._lon; }
    get lng()       { return this._lon; }
    get longitude() { return this._lon; }
    set lon(lon) {
        this._lon = isNaN(lon) ? Dms.wrap180(Dms.parse(lon)) : Dms.wrap180(Number(lon));
        if (isNaN(this._lon)) throw new TypeError(`invalid lon ‘${lon}’`);
    }
    set lng(lon) {
        this._lon = isNaN(lon) ? Dms.wrap180(Dms.parse(lon)) : Dms.wrap180(Number(lon));
        if (isNaN(this._lon)) throw new TypeError(`invalid lng ‘${lon}’`);
    }
    set longitude(lon) {
        this._lon = isNaN(lon) ? Dms.wrap180(Dms.parse(lon)) : Dms.wrap180(Number(lon));
        if (isNaN(this._lon)) throw new TypeError(`invalid longitude ‘${lon}’`);
    }


    /** Conversion factors; 1000 * LatLon.metresToKm gives 1. */
    static get metresToKm()            { return 1/1000; }
    /** Conversion factors; 1000 * LatLon.metresToMiles gives 0.621371192237334. */
    static get metresToMiles()         { return 1/1609.344; }
    /** Conversion factors; 1000 * LatLon.metresToMiles gives 0.5399568034557236. */
    static get metresToNauticalMiles() { return 1/1852; }


    /**
     * Parses a latitude/longitude point from a variety of formats.
     *
     * Latitude & longitude (in degrees) can be supplied as two separate parameters, as a single
     * comma-separated lat/lon string, or as a single object with { lat, lon } or GeoJSON properties.
     *
     * The latitude/longitude values may be numeric or strings; they may be signed decimal or
     * deg-min-sec (hexagesimal) suffixed by compass direction (NSEW); a variety of separators are
     * accepted. Examples -3.62, '3 37 12W', '3°37′12″W'.
     *
     * Thousands/decimal separators must be comma/dot; use Dms.fromLocale to convert locale-specific
     * thousands/decimal separators.
     *
     * @param   {number|string|Object} lat|latlon - Latitude (in degrees) or comma-separated lat/lon or lat/lon object.
     * @param   {number|string}        [lon]      - Longitude (in degrees).
     * @returns {LatLon} Latitude/longitude point.
     * @throws  {TypeError} Invalid point.
     *
     * @example
     *   const p1 = LatLon.parse(52.205, 0.119);                                    // numeric pair (≡ new LatLon)
     *   const p2 = LatLon.parse('52.205', '0.119');                                // numeric string pair (≡ new LatLon)
     *   const p3 = LatLon.parse('52.205, 0.119');                                  // single string numerics
     *   const p4 = LatLon.parse('52°12′18.0″N', '000°07′08.4″E');                  // DMS pair
     *   const p5 = LatLon.parse('52°12′18.0″N, 000°07′08.4″E');                    // single string DMS
     *   const p6 = LatLon.parse({ lat: 52.205, lon: 0.119 });                      // { lat, lon } object numeric
     *   const p7 = LatLon.parse({ lat: '52°12′18.0″N', lng: '000°07′08.4″E' });    // { lat, lng } object DMS
     *   const p8 = LatLon.parse({ type: 'Point', coordinates: [ 0.119, 52.205] }); // GeoJSON
     */
    static parse(...args) {
        if (args.length == 0) throw new TypeError('invalid (empty) point');
        if (args[0]===null || args[1]===null) throw new TypeError('invalid (null) point');

        let lat=undefined, lon=undefined;

        if (args.length == 2) { // regular (lat, lon) arguments
            [ lat, lon ] = args;
            lat = Dms.wrap90(Dms.parse(lat));
            lon = Dms.wrap180(Dms.parse(lon));
            if (isNaN(lat) || isNaN(lon)) throw new TypeError(`invalid point ‘${args.toString()}’`);
        }

        if (args.length == 1 && typeof args[0] == 'string') { // single comma-separated lat,lon string
            [ lat, lon ] = args[0].split(',');
            lat = Dms.wrap90(Dms.parse(lat));
            lon = Dms.wrap180(Dms.parse(lon));
            if (isNaN(lat) || isNaN(lon)) throw new TypeError(`invalid point ‘${args[0]}’`);
        }

        if (args.length == 1 && typeof args[0] == 'object') { // single { lat, lon } object
            const ll = args[0];
            if (ll.type == 'Point' && Array.isArray(ll.coordinates)) { // GeoJSON
                [ lon, lat ] = ll.coordinates;
            } else { // regular { lat, lon } object
                if (ll.latitude  != undefined) lat = ll.latitude;
                if (ll.lat       != undefined) lat = ll.lat;
                if (ll.longitude != undefined) lon = ll.longitude;
                if (ll.lng       != undefined) lon = ll.lng;
                if (ll.lon       != undefined) lon = ll.lon;
                lat = Dms.wrap90(Dms.parse(lat));
                lon = Dms.wrap180(Dms.parse(lon));
            }
            if (isNaN(lat) || isNaN(lon)) throw new TypeError(`invalid point ‘${JSON.stringify(args[0])}’`);
        }

        if (isNaN(lat) || isNaN(lon)) throw new TypeError(`invalid point ‘${args.toString()}’`);

        return new LatLonSpherical(lat, lon);
    }


    /**
     * Returns the distance along the surface of the earth from ‘this’ point to destination point.
     *
     * Uses haversine formula: a = sin²(deltaphi/2) + cosphi1·cosphi2 · sin²(deltalambda/2); d = 2 · atan2(√a, √(a-1)).
     *
     * @param   {LatLon} point - Latitude/longitude of destination point.
     * @param   {number} [radius=6371e3] - Radius of earth (defaults to mean radius in metres).
     * @returns {number} Distance between this point and destination point, in same units as radius.
     * @throws  {TypeError} Invalid radius.
     *
     * @example
     *   const p1 = new LatLon(52.205, 0.119);
     *   const p2 = new LatLon(48.857, 2.351);
     *   const d = p1.distanceTo(p2);       // 404.3×10³ m
     *   const m = p1.distanceTo(p2, 3959); // 251.2 miles
     */
    distanceTo(point, radius=6371e3) {
        if (!(point instanceof LatLonSpherical)) point = LatLonSpherical.parse(point); // allow literal forms
        if (isNaN(radius)) throw new TypeError(`invalid radius ‘${radius}’`);

        // a = sin²(deltaphi/2) + cos(phi1)⋅cos(phi2)⋅sin²(deltalambda/2)
        // delta = 2·atan2(√(a), √(1−a))
        // see mathforum.org/library/drmath/view/51879.html for derivation

        const R = radius;
        const phi1 = this.lat.toRadians(),  lambda1 = this.lon.toRadians();
        const phi2 = point.lat.toRadians(), lambda2 = point.lon.toRadians();
        const deltaphi = phi2 - phi1;
        const deltalambda = lambda2 - lambda1;

        const a = Math.sin(deltaphi/2)*Math.sin(deltaphi/2) + Math.cos(phi1)*Math.cos(phi2) * Math.sin(deltalambda/2)*Math.sin(deltalambda/2);
        const c = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1-a));
        const d = R * c;

        return d;
    }


    /**
     * Returns the initial bearing from ‘this’ point to destination point.
     *
     * @param   {LatLon} point - Latitude/longitude of destination point.
     * @returns {number} Initial bearing in degrees from north (0°..360°).
     *
     * @example
     *   const p1 = new LatLon(52.205, 0.119);
     *   const p2 = new LatLon(48.857, 2.351);
     *   const b1 = p1.initialBearingTo(p2); // 156.2°
     */
    initialBearingTo(point) {
        if (!(point instanceof LatLonSpherical)) point = LatLonSpherical.parse(point); // allow literal forms
        if (this.equals(point)) return NaN; // coincident points

        // tantheta = sindeltalambda⋅cosphi2 / cosphi1⋅sinphi2 − sinphi1⋅cosphi2⋅cosdeltalambda
        // see mathforum.org/library/drmath/view/55417.html for derivation

        const phi1 = this.lat.toRadians();
        const phi2 = point.lat.toRadians();
        const deltalambda = (point.lon - this.lon).toRadians();

        const x = Math.cos(phi1) * Math.sin(phi2) - Math.sin(phi1) * Math.cos(phi2) * Math.cos(deltalambda);
        const y = Math.sin(deltalambda) * Math.cos(phi2);
        const theta = Math.atan2(y, x);

        const bearing = theta.toDegrees();

        return Dms.wrap360(bearing);
    }


    /**
     * Returns final bearing arriving at destination point from ‘this’ point; the final bearing will
     * differ from the initial bearing by varying degrees according to distance and latitude.
     *
     * @param   {LatLon} point - Latitude/longitude of destination point.
     * @returns {number} Final bearing in degrees from north (0°..360°).
     *
     * @example
     *   const p1 = new LatLon(52.205, 0.119);
     *   const p2 = new LatLon(48.857, 2.351);
     *   const b2 = p1.finalBearingTo(p2); // 157.9°
     */
    finalBearingTo(point) {
        if (!(point instanceof LatLonSpherical)) point = LatLonSpherical.parse(point); // allow literal forms

        // get initial bearing from destination point to this point & reverse it by adding 180°

        const bearing = point.initialBearingTo(this) + 180;

        return Dms.wrap360(bearing);
    }


    /**
     * Returns the midpoint between ‘this’ point and destination point.
     *
     * @param   {LatLon} point - Latitude/longitude of destination point.
     * @returns {LatLon} Midpoint between this point and destination point.
     *
     * @example
     *   const p1 = new LatLon(52.205, 0.119);
     *   const p2 = new LatLon(48.857, 2.351);
     *   const pMid = p1.midpointTo(p2); // 50.5363°N, 001.2746°E
     */
    midpointTo(point) {
        if (!(point instanceof LatLonSpherical)) point = LatLonSpherical.parse(point); // allow literal forms

        // phim = atan2( sinphi1 + sinphi2, √( (cosphi1 + cosphi2⋅cosdeltalambda)² + cos²phi2⋅sin²deltalambda ) )
        // lambdam = lambda1 + atan2(cosphi2⋅sindeltalambda, cosphi1 + cosphi2⋅cosdeltalambda)
        // midpoint is sum of vectors to two points: mathforum.org/library/drmath/view/51822.html

        const phi1 = this.lat.toRadians();
        const lambda1 = this.lon.toRadians();
        const phi2 = point.lat.toRadians();
        const deltalambda = (point.lon - this.lon).toRadians();

        // get cartesian coordinates for the two points
        const A = { x: Math.cos(phi1), y: 0, z: Math.sin(phi1) }; // place point A on prime meridian y=0
        const B = { x: Math.cos(phi2)*Math.cos(deltalambda), y: Math.cos(phi2)*Math.sin(deltalambda), z: Math.sin(phi2) };

        // vector to midpoint is sum of vectors to two points (no need to normalise)
        const C = { x: A.x + B.x, y: A.y + B.y, z: A.z + B.z };

        const phim = Math.atan2(C.z, Math.sqrt(C.x*C.x + C.y*C.y));
        const lambdam = lambda1 + Math.atan2(C.y, C.x);

        const lat = phim.toDegrees();
        const lon = lambdam.toDegrees();

        return new LatLonSpherical(lat, lon);
    }


    /**
     * Returns the point at given fraction between ‘this’ point and given point.
     *
     * @param   {LatLon} point - Latitude/longitude of destination point.
     * @param   {number} fraction - Fraction between the two points (0 = this point, 1 = specified point).
     * @returns {LatLon} Intermediate point between this point and destination point.
     *
     * @example
     *   const p1 = new LatLon(52.205, 0.119);
     *   const p2 = new LatLon(48.857, 2.351);
     *   const pInt = p1.intermediatePointTo(p2, 0.25); // 51.3721°N, 000.7073°E
     */
    intermediatePointTo(point, fraction) {
        if (!(point instanceof LatLonSpherical)) point = LatLonSpherical.parse(point); // allow literal forms
        if (this.equals(point)) return new LatLonSpherical(this.lat, this.lon); // coincident points

        const phi1 = this.lat.toRadians(), lambda1 = this.lon.toRadians();
        const phi2 = point.lat.toRadians(), lambda2 = point.lon.toRadians();

        // distance between points
        const deltaphi = phi2 - phi1;
        const deltalambda = lambda2 - lambda1;
        const a = Math.sin(deltaphi/2) * Math.sin(deltaphi/2)
            + Math.cos(phi1) * Math.cos(phi2) * Math.sin(deltalambda/2) * Math.sin(deltalambda/2);
        const delta = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1-a));

        const A = Math.sin((1-fraction)*delta) / Math.sin(delta);
        const B = Math.sin(fraction*delta) / Math.sin(delta);

        const x = A * Math.cos(phi1) * Math.cos(lambda1) + B * Math.cos(phi2) * Math.cos(lambda2);
        const y = A * Math.cos(phi1) * Math.sin(lambda1) + B * Math.cos(phi2) * Math.sin(lambda2);
        const z = A * Math.sin(phi1) + B * Math.sin(phi2);

        const phi3 = Math.atan2(z, Math.sqrt(x*x + y*y));
        const lambda3 = Math.atan2(y, x);

        const lat = phi3.toDegrees();
        const lon = lambda3.toDegrees();

        return new LatLonSpherical(lat, lon);
    }


    /**
     * Returns the destination point from ‘this’ point having travelled the given distance on the
     * given initial bearing (bearing normally varies around path followed).
     *
     * @param   {number} distance - Distance travelled, in same units as earth radius (default: metres).
     * @param   {number} bearing - Initial bearing in degrees from north.
     * @param   {number} [radius=6371e3] - (Mean) radius of earth (defaults to radius in metres).
     * @returns {LatLon} Destination point.
     *
     * @example
     *   const p1 = new LatLon(51.47788, -0.00147);
     *   const p2 = p1.destinationPoint(7794, 300.7); // 51.5136°N, 000.0983°W
     */
    destinationPoint(distance, bearing, radius=6371e3) {
        // sinphi2 = sinphi1⋅cosdelta + cosphi1⋅sindelta⋅costheta
        // tandeltalambda = sintheta⋅sindelta⋅cosphi1 / cosdelta−sinphi1⋅sinphi2
        // see mathforum.org/library/drmath/view/52049.html for derivation

        const delta = distance / radius; // angular distance in radians
        const theta = Number(bearing).toRadians();

        const phi1 = this.lat.toRadians(), lambda1 = this.lon.toRadians();

        const sinphi2 = Math.sin(phi1) * Math.cos(delta) + Math.cos(phi1) * Math.sin(delta) * Math.cos(theta);
        const phi2 = Math.asin(sinphi2);
        const y = Math.sin(theta) * Math.sin(delta) * Math.cos(phi1);
        const x = Math.cos(delta) - Math.sin(phi1) * sinphi2;
        const lambda2 = lambda1 + Math.atan2(y, x);

        const lat = phi2.toDegrees();
        const lon = lambda2.toDegrees();

        return new LatLonSpherical(lat, lon);
    }


    /**
     * Returns the point of intersection of two paths defined by point and bearing.
     *
     * @param   {LatLon}      p1 - prime point.
     * @param   {number}      brng1 - Initial bearing from prime point.
     * @param   {LatLon}      p2 - Second point.
     * @param   {number}      brng2 - Initial bearing from second point.
     * @returns {LatLon|null} Destination point (null if no unique intersection defined).
     *
     * @example
     *   const p1 = new LatLon(51.8853, 0.2545), brng1 = 108.547;
     *   const p2 = new LatLon(49.0034, 2.5735), brng2 =  32.435;
     *   const pInt = LatLon.intersection(p1, brng1, p2, brng2); // 50.9078°N, 004.5084°E
     */
    static intersection(p1, brng1, p2, brng2) {
        if (!(p1 instanceof LatLonSpherical)) p1 = LatLonSpherical.parse(p1); // allow literal forms
        if (!(p2 instanceof LatLonSpherical)) p2 = LatLonSpherical.parse(p2); // allow literal forms
        if (isNaN(brng1)) throw new TypeError(`invalid brng1 ‘${brng1}’`);
        if (isNaN(brng2)) throw new TypeError(`invalid brng2 ‘${brng2}’`);

        // see www.edwilliams.org/avform.htm#Intersection

        const phi1 = p1.lat.toRadians(), lambda1 = p1.lon.toRadians();
        const phi2 = p2.lat.toRadians(), lambda2 = p2.lon.toRadians();
        const theta13 = Number(brng1).toRadians(), theta23 = Number(brng2).toRadians();
        const deltaphi = phi2 - phi1, deltalambda = lambda2 - lambda1;

        // angular distance p1-p2
        const delta12 = 2 * Math.asin(Math.sqrt(Math.sin(deltaphi/2) * Math.sin(deltaphi/2)
            + Math.cos(phi1) * Math.cos(phi2) * Math.sin(deltalambda/2) * Math.sin(deltalambda/2)));
        if (Math.abs(delta12) < Number.EPSILON) return new LatLonSpherical(p1.lat, p1.lon); // coincident points

        // initial/final bearings between points
        const costhetaa = (Math.sin(phi2) - Math.sin(phi1)*Math.cos(delta12)) / (Math.sin(delta12)*Math.cos(phi1));
        const costhetab = (Math.sin(phi1) - Math.sin(phi2)*Math.cos(delta12)) / (Math.sin(delta12)*Math.cos(phi2));
        const thetaa = Math.acos(Math.min(Math.max(costhetaa, -1), 1)); // protect against rounding errors
        const thetab = Math.acos(Math.min(Math.max(costhetab, -1), 1)); // protect against rounding errors

        const theta12 = Math.sin(lambda2-lambda1)>0 ? thetaa : 2*pi-thetaa;
        const theta21 = Math.sin(lambda2-lambda1)>0 ? 2*pi-thetab : thetab;

        const alpha1 = theta13 - theta12; // angle 2-1-3
        const alpha2 = theta21 - theta23; // angle 1-2-3

        if (Math.sin(alpha1) == 0 && Math.sin(alpha2) == 0) return null; // infinite intersections
        if (Math.sin(alpha1) * Math.sin(alpha2) < 0) return null;        // ambiguous intersection (antipodal?)

        const cosalpha3 = -Math.cos(alpha1)*Math.cos(alpha2) + Math.sin(alpha1)*Math.sin(alpha2)*Math.cos(delta12);

        const delta13 = Math.atan2(Math.sin(delta12)*Math.sin(alpha1)*Math.sin(alpha2), Math.cos(alpha2) + Math.cos(alpha1)*cosalpha3);

        const phi3 = Math.asin(Math.min(Math.max(Math.sin(phi1)*Math.cos(delta13) + Math.cos(phi1)*Math.sin(delta13)*Math.cos(theta13), -1), 1));

        const deltalambda13 = Math.atan2(Math.sin(theta13)*Math.sin(delta13)*Math.cos(phi1), Math.cos(delta13) - Math.sin(phi1)*Math.sin(phi3));
        const lambda3 = lambda1 + deltalambda13;

        const lat = phi3.toDegrees();
        const lon = lambda3.toDegrees();

        return new LatLonSpherical(lat, lon);
    }


    /**
     * Returns (signed) distance from ‘this’ point to great circle defined by start-point and
     * end-point.
     *
     * @param   {LatLon} pathStart - Start point of great circle path.
     * @param   {LatLon} pathEnd - End point of great circle path.
     * @param   {number} [radius=6371e3] - (Mean) radius of earth (defaults to radius in metres).
     * @returns {number} Distance to great circle (-ve if to left, +ve if to right of path).
     *
     * @example
     *   const pCurrent = new LatLon(53.2611, -0.7972);
     *   const p1 = new LatLon(53.3206, -1.7297);
     *   const p2 = new LatLon(53.1887, 0.1334);
     *   const d = pCurrent.crossTrackDistanceTo(p1, p2);  // -307.5 m
     */
    crossTrackDistanceTo(pathStart, pathEnd, radius=6371e3) {
        if (!(pathStart instanceof LatLonSpherical)) pathStart = LatLonSpherical.parse(pathStart); // allow literal forms
        if (!(pathEnd instanceof LatLonSpherical)) pathEnd = LatLonSpherical.parse(pathEnd);       // allow literal forms
        const R = radius;

        if (this.equals(pathStart)) return 0;

        const delta13 = pathStart.distanceTo(this, R) / R;
        const theta13 = pathStart.initialBearingTo(this).toRadians();
        const theta12 = pathStart.initialBearingTo(pathEnd).toRadians();

        const deltaxt = Math.asin(Math.sin(delta13) * Math.sin(theta13 - theta12));

        return deltaxt * R;
    }


    /**
     * Returns how far ‘this’ point is along a path from from start-point, heading towards end-point.
     * That is, if a perpendicular is drawn from ‘this’ point to the (great circle) path, the
     * along-track distance is the distance from the start point to where the perpendicular crosses
     * the path.
     *
     * @param   {LatLon} pathStart - Start point of great circle path.
     * @param   {LatLon} pathEnd - End point of great circle path.
     * @param   {number} [radius=6371e3] - (Mean) radius of earth (defaults to radius in metres).
     * @returns {number} Distance along great circle to point nearest ‘this’ point.
     *
     * @example
     *   const pCurrent = new LatLon(53.2611, -0.7972);
     *   const p1 = new LatLon(53.3206, -1.7297);
     *   const p2 = new LatLon(53.1887,  0.1334);
     *   const d = pCurrent.alongTrackDistanceTo(p1, p2);  // 62.331 km
     */
    alongTrackDistanceTo(pathStart, pathEnd, radius=6371e3) {
        if (!(pathStart instanceof LatLonSpherical)) pathStart = LatLonSpherical.parse(pathStart); // allow literal forms
        if (!(pathEnd instanceof LatLonSpherical)) pathEnd = LatLonSpherical.parse(pathEnd);       // allow literal forms
        const R = radius;

        if (this.equals(pathStart)) return 0;

        const delta13 = pathStart.distanceTo(this, R) / R;
        const theta13 = pathStart.initialBearingTo(this).toRadians();
        const theta12 = pathStart.initialBearingTo(pathEnd).toRadians();

        const deltaxt = Math.asin(Math.sin(delta13) * Math.sin(theta13-theta12));

        const deltaat = Math.acos(Math.cos(delta13) / Math.abs(Math.cos(deltaxt)));

        return deltaat*Math.sign(Math.cos(theta12-theta13)) * R;
    }


    /**
     * Returns maximum latitude reached when travelling on a great circle on given bearing from
     * ‘this’ point (‘Clairaut’s formula’). Negate the result for the minimum latitude (in the
     * southern hemisphere).
     *
     * The maximum latitude is independent of longitude; it will be the same for all points on a
     * given latitude.
     *
     * @param   {number} bearing - Initial bearing.
     * @returns {number} Maximum latitude reached.
     */
    maxLatitude(bearing) {
        const theta = Number(bearing).toRadians();

        const phi = this.lat.toRadians();

        const phiMax = Math.acos(Math.abs(Math.sin(theta) * Math.cos(phi)));

        return phiMax.toDegrees();
    }


    /**
     * Returns the pair of meridians at which a great circle defined by two points crosses the given
     * latitude. If the great circle doesn't reach the given latitude, null is returned.
     *
     * @param   {LatLon}      point1 - prime point defining great circle.
     * @param   {LatLon}      point2 - Second point defining great circle.
     * @param   {number}      latitude - Latitude crossings are to be determined for.
     * @returns {Object|null} Object containing { lon1, lon2 } or null if given latitude not reached.
     */
    static crossingParallels(point1, point2, latitude) {
        if (point1.equals(point2)) return null; // coincident points

        const phi = Number(latitude).toRadians();

        const phi1 = point1.lat.toRadians();
        const lambda1 = point1.lon.toRadians();
        const phi2 = point2.lat.toRadians();
        const lambda2 = point2.lon.toRadians();

        const deltalambda = lambda2 - lambda1;

        const x = Math.sin(phi1) * Math.cos(phi2) * Math.cos(phi) * Math.sin(deltalambda);
        const y = Math.sin(phi1) * Math.cos(phi2) * Math.cos(phi) * Math.cos(deltalambda) - Math.cos(phi1) * Math.sin(phi2) * Math.cos(phi);
        const z = Math.cos(phi1) * Math.cos(phi2) * Math.sin(phi) * Math.sin(deltalambda);

        if (z * z > x * x + y * y) return null; // great circle doesn't reach latitude

        const lambdam = Math.atan2(-y, x);               // longitude at max latitude
        const deltalambdai = Math.acos(z / Math.sqrt(x*x + y*y)); // deltalambda from lambdam to intersection points

        const lambdai1 = lambda1 + lambdam - deltalambdai;
        const lambdai2 = lambda1 + lambdam + deltalambdai;

        const lon1 = lambdai1.toDegrees();
        const lon2 = lambdai2.toDegrees();

        return {
            lon1: Dms.wrap180(lon1),
            lon2: Dms.wrap180(lon2),
        };
    }


    /* Rhumb - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */


    /**
     * Returns the distance travelling from ‘this’ point to destination point along a rhumb line.
     *
     * @param   {LatLon} point - Latitude/longitude of destination point.
     * @param   {number} [radius=6371e3] - (Mean) radius of earth (defaults to radius in metres).
     * @returns {number} Distance in km between this point and destination point (same units as radius).
     *
     * @example
     *   const p1 = new LatLon(51.127, 1.338);
     *   const p2 = new LatLon(50.964, 1.853);
     *   const d = p1.distanceTo(p2); //  40.31 km
     */
    rhumbDistanceTo(point, radius=6371e3) {
        if (!(point instanceof LatLonSpherical)) point = LatLonSpherical.parse(point); // allow literal forms

        // see www.edwilliams.org/avform.htm#Rhumb

        const R = radius;
        const phi1 = this.lat.toRadians();
        const phi2 = point.lat.toRadians();
        const deltaphi = phi2 - phi1;
        let deltalambda = Math.abs(point.lon - this.lon).toRadians();
        // if dLon over 180° take shorter rhumb line across the anti-meridian:
        if (Math.abs(deltalambda) > pi) deltalambda = deltalambda > 0 ? -(2 * pi - deltalambda) : (2 * pi + deltalambda);

        // on Mercator projection, longitude distances shrink by latitude; q is the 'stretch factor'
        // q becomes ill-conditioned along E-W line (0/0); use empirical tolerance to avoid it
        const deltapsi = Math.log(Math.tan(phi2 / 2 + pi / 4) / Math.tan(phi1 / 2 + pi / 4));
        const q = Math.abs(deltapsi) > 10e-12 ? deltaphi / deltapsi : Math.cos(phi1);

        // distance is pythagoras on 'stretched' Mercator projection, √(deltaphi² + q²·deltalambda²)
        const delta = Math.sqrt(deltaphi*deltaphi + q*q * deltalambda*deltalambda); // angular distance in radians
        const d = delta * R;

        return d;
    }


    /**
     * Returns the bearing from ‘this’ point to destination point along a rhumb line.
     *
     * @param   {LatLon}    point - Latitude/longitude of destination point.
     * @returns {number}    Bearing in degrees from north.
     *
     * @example
     *   const p1 = new LatLon(51.127, 1.338);
     *   const p2 = new LatLon(50.964, 1.853);
     *   const d = p1.rhumbBearingTo(p2); // 116.7°
     */
    rhumbBearingTo(point) {
        if (!(point instanceof LatLonSpherical)) point = LatLonSpherical.parse(point); // allow literal forms
        if (this.equals(point)) return NaN; // coincident points

        const phi1 = this.lat.toRadians();
        const phi2 = point.lat.toRadians();
        let deltalambda = (point.lon - this.lon).toRadians();
        // if dLon over 180° take shorter rhumb line across the anti-meridian:
        if (Math.abs(deltalambda) > pi) deltalambda = deltalambda > 0 ? -(2 * pi - deltalambda) : (2 * pi + deltalambda);

        const deltapsi = Math.log(Math.tan(phi2 / 2 + pi / 4) / Math.tan(phi1 / 2 + pi / 4));

        const theta = Math.atan2(deltalambda, deltapsi);

        const bearing = theta.toDegrees();

        return Dms.wrap360(bearing);
    }


    /**
     * Returns the destination point having travelled along a rhumb line from ‘this’ point the given
     * distance on the given bearing.
     *
     * @param   {number} distance - Distance travelled, in same units as earth radius (default: metres).
     * @param   {number} bearing - Bearing in degrees from north.
     * @param   {number} [radius=6371e3] - (Mean) radius of earth (defaults to radius in metres).
     * @returns {LatLon} Destination point.
     *
     * @example
     *   const p1 = new LatLon(51.127, 1.338);
     *   const p2 = p1.rhumbDestinationPoint(40300, 116.7); // 50.9642°N, 001.8530°E
     */
    rhumbDestinationPoint(distance, bearing, radius=6371e3) {
        const phi1 = this.lat.toRadians(), lambda1 = this.lon.toRadians();
        const theta = Number(bearing).toRadians();

        const delta = distance / radius; // angular distance in radians

        const deltaphi = delta * Math.cos(theta);
        let phi2 = phi1 + deltaphi;

        // check for some daft bugger going past the pole, normalise latitude if so
        if (Math.abs(phi2) > pi / 2) phi2 = phi2 > 0 ? pi - phi2 : -pi - phi2;

        const deltapsi = Math.log(Math.tan(phi2 / 2 + pi / 4) / Math.tan(phi1 / 2 + pi / 4));
        const q = Math.abs(deltapsi) > 10e-12 ? deltaphi / deltapsi : Math.cos(phi1); // E-W course becomes ill-conditioned with 0/0

        const deltalambda = delta * Math.sin(theta) / q;
        const lambda2 = lambda1 + deltalambda;

        const lat = phi2.toDegrees();
        const lon = lambda2.toDegrees();

        return new LatLonSpherical(lat, lon);
    }


    /**
     * Returns the loxodromic midpoint (along a rhumb line) between ‘this’ point and second point.
     *
     * @param   {LatLon} point - Latitude/longitude of second point.
     * @returns {LatLon} Midpoint between this point and second point.
     *
     * @example
     *   const p1 = new LatLon(51.127, 1.338);
     *   const p2 = new LatLon(50.964, 1.853);
     *   const pMid = p1.rhumbMidpointTo(p2); // 51.0455°N, 001.5957°E
     */
    rhumbMidpointTo(point) {
        if (!(point instanceof LatLonSpherical)) point = LatLonSpherical.parse(point); // allow literal forms

        // see mathforum.org/kb/message.jspa?messageID=148837

        const phi1 = this.lat.toRadians(); let lambda1 = this.lon.toRadians();
        const phi2 = point.lat.toRadians(), lambda2 = point.lon.toRadians();

        if (Math.abs(lambda2 - lambda1) > pi) lambda1 += 2 * pi; // crossing anti-meridian

        const phi3 = (phi1 + phi2) / 2;
        const f1 = Math.tan(pi / 4 + phi1 / 2);
        const f2 = Math.tan(pi / 4 + phi2 / 2);
        const f3 = Math.tan(pi / 4 + phi3 / 2);
        let lambda3 = ((lambda2 - lambda1) * Math.log(f3) + lambda1 * Math.log(f2) - lambda2 * Math.log(f1)) / Math.log(f2 / f1);

        if (!isFinite(lambda3)) lambda3 = (lambda1 + lambda2) / 2; // parallel of latitude

        const lat = phi3.toDegrees();
        const lon = lambda3.toDegrees();

        return new LatLonSpherical(lat, lon);
    }


    /* Area - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */


    /**
     * Calculates the area of a spherical polygon where the sides of the polygon are great circle
     * arcs joining the vertices.
     *
     * @param   {LatLon[]} polygon - Array of points defining vertices of the polygon.
     * @param   {number}   [radius=6371e3] - (Mean) radius of earth (defaults to radius in metres).
     * @returns {number}   The area of the polygon in the same units as radius.
     *
     * @example
     *   const polygon = [new LatLon(0,0), new LatLon(1,0), new LatLon(0,1)];
     *   const area = LatLon.areaOf(polygon); // 6.18e9 m²
     */
    static areaOf(polygon, radius=6371e3) {
        // uses method due to Karney: osgeo-org.1560.x6.nabble.com/Area-of-a-spherical-polygon-td3841625.html;
        // for each edge of the polygon, tan(E/2) = tan(deltalambda/2)·(tan(phi₁/2)+tan(phi₂/2)) / (1+tan(phi₁/2)·tan(phi₂/2))
        // where E is the spherical excess of the trapezium obtained by extending the edge to the equator
        // (Karney's method is probably more efficient than the more widely known L’Huilier’s Theorem)

        const R = radius;

        // close polygon so that last point equals prime point
        const closed = polygon[0].equals(polygon[polygon.length-1]);
        if (!closed) polygon.push(polygon[0]);

        const nVertices = polygon.length - 1;

        let S = 0; // spherical excess in steradians
        for (let v=0; v<nVertices; v++) {
            const phi1 = polygon[v].lat.toRadians();
            const phi2 = polygon[v+1].lat.toRadians();
            const deltalambda = (polygon[v+1].lon - polygon[v].lon).toRadians();
            const E = 2 * Math.atan2(Math.tan(deltalambda/2) * (Math.tan(phi1/2)+Math.tan(phi2/2)), 1 + Math.tan(phi1/2)*Math.tan(phi2/2));
            S += E;
        }

        if (isPoleEnclosedBy(polygon)) S = Math.abs(S) - 2*pi;

        const A = Math.abs(S * R*R); // area in units of R

        if (!closed) polygon.pop(); // restore polygon to pristine condition

        return A;

        // returns whether polygon encloses pole: sum of course deltas around pole is 0° rather than
        // normal ±360°: blog.element84.com/determining-if-a-spherical-polygon-contains-a-pole.html
        function isPoleEnclosedBy(p) {
            // TODO: any better test than this?
            let sigmadelta = 0;
            let prevBrng = p[0].initialBearingTo(p[1]);
            for (let v=0; v<p.length-1; v++) {
                const initBrng = p[v].initialBearingTo(p[v+1]);
                const finalBrng = p[v].finalBearingTo(p[v+1]);
                sigmadelta += (initBrng - prevBrng + 540) % 360 - 180;
                sigmadelta += (finalBrng - initBrng + 540) % 360 - 180;
                prevBrng = finalBrng;
            }
            const initBrng = p[0].initialBearingTo(p[1]);
            sigmadelta += (initBrng - prevBrng + 540) % 360 - 180;
            // TODO: fix (intermittant) edge crossing pole - eg (85,90), (85,0), (85,-90)
            const enclosed = Math.abs(sigmadelta) < 90; // 0°-ish
            return enclosed;
        }
    }


    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */


    /**
     * Checks if another point is equal to ‘this’ point.
     *
     * @param   {LatLon} point - Point to be compared against this point.
     * @returns {bool}   True if points have identical latitude and longitude values.
     *
     * @example
     *   const p1 = new LatLon(52.205, 0.119);
     *   const p2 = new LatLon(52.205, 0.119);
     *   const equal = p1.equals(p2); // true
     */
    equals(point) {
        if (!(point instanceof LatLonSpherical)) point = LatLonSpherical.parse(point); // allow literal forms

        if (Math.abs(this.lat - point.lat) > Number.EPSILON) return false;
        if (Math.abs(this.lon - point.lon) > Number.EPSILON) return false;

        return true;
    }


    /**
     * Converts ‘this’ point to a GeoJSON object.
     *
     * @returns {Object} this point as a GeoJSON ‘Point’ object.
     */
    toGeoJSON() {
        return { type: 'Point', coordinates: [ this.lon, this.lat ] };
    }


    /**
     * Returns a string representation of ‘this’ point, formatted as degrees, degrees+minutes, or
     * degrees+minutes+seconds.
     *
     * @param   {string} [format=d] - Format point as 'd', 'dm', 'dms', or 'n' for signed numeric.
     * @param   {number} [dp=4|2|0] - Number of decimal places to use: default 4 for d, 2 for dm, 0 for dms.
     * @returns {string} Comma-separated formatted latitude/longitude.
     * @throws  {RangeError} Invalid format.
     *
     * @example
     *   const greenwich = new LatLon(51.47788, -0.00147);
     *   const d = greenwich.toString();                        // 51.4779°N, 000.0015°W
     *   const dms = greenwich.toString('dms', 2);              // 51°28′40.37″N, 000°00′05.29″W
     *   const [lat, lon] = greenwich.toString('n').split(','); // 51.4779, -0.0015
     */
    toString(format='d', dp=undefined) {
        // note: explicitly set dp to undefined for passing through to toLat/toLon
        if (![ 'd', 'dm', 'dms', 'n' ].includes(format)) throw new RangeError(`invalid format ‘${format}’`);

        if (format == 'n') { // signed numeric degrees
            if (dp == undefined) dp = 4;
            return `${this.lat.toFixed(dp)},${this.lon.toFixed(dp)}`;
        }
        const lat = Dms.toLat(this.lat, format, dp);
        const lon = Dms.toLon(this.lon, format, dp);
        return `${lat}, ${lon}`;
    }

}


/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */

export { LatLonSpherical as default, Dms };
