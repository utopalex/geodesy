/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */
/* Ordnance Survey Grid Reference functions                           (c) Chris Veness 2005-2019  */
/*                                                                                   MIT Licence  */
/* www.movable-type.co.uk/scripts/latlong-gridref.html                                            */
/* www.movable-type.co.uk/scripts/geodesy-library.html#osgridref                                  */
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */

import LatLonEllipsoidal, { Dms } from './latlon-ellipsoidal-datum.js';

/**
 * Ordnance Survey OSGB grid references provide geocoordinate references for UK mapping purposes.
 *
 * Formulation implemented here due to Thomas, Redfearn, etc is as published by OS, but is inferior
 * to Krüger as used by e.g. Karney 2011.
 *
 * www.ordnancesurvey.co.uk/docs/support/guide-coordinate-systems-great-britain.pdf.
 *
 * Note OSGB grid references cover Great Britain only; Ireland and the Channel Islands have their
 * own references.
 *
 * Note that these formulae are based on ellipsoidal calculations, and according to the OS are
 * accurate to about 4–5 metres – for greater accuracy, a geoid-based transformation (OSTN15) must
 * be used.
 */

/*
 * Converted 2015 to work with WGS84 by default, OSGB36 as option;
 * www.ordnancesurvey.co.uk/blog/2014/12/confirmation-on-changes-to-latitude-and-longitude
 */


/* OsGridRef  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */


/**
 * OS grid references with methods to parse and convert them to latitude/longitude points.
 *
 * @module osgridref
 */
class OsGridRef {

    /**
     * Creates an OsGridRef object.
     *
     * @param {number} easting - Easting in metres from OS false origin.
     * @param {number} northing - Northing in metres from OS false origin.
     *
     * @example
     *   import OsGridRef from '/js/geodesy/osgridref.js';
     *   const gridref = new OsGridRef(651409, 313177);
     */
    constructor(easting, northing) {
        this.easting = Number(easting);
        this.northing = Number(northing);

        if (isNaN(easting)  || this.easting<0  || this.easting>700e3) throw new RangeError(`invalid easting ‘${easting}’`);
        if (isNaN(northing) || this.northing<0 || this.northing>1300e3) throw new RangeError(`invalid northing ‘${northing}’`);
    }


    /**
     * Converts ‘this’ Ordnance Survey grid reference easting/northing coordinate to latitude/longitude
     * (SW corner of grid square).
     *
     * While OS grid references are based on OSGB-36, the Ordnance Survey have deprecated the use of
     * OSGB-36 for latitude/longitude coordinates (in favour of WGS-84), hence this function returns
     * WGS-84 by default, with OSGB-36 as an option. See www.ordnancesurvey.co.uk/blog/2014/12/2.
     *
     * Note formulation implemented here due to Thomas, Redfearn, etc is as published by OS, but is
     * inferior to Krüger as used by e.g. Karney 2011.
     *
     * @param   {LatLon.datum} [datum=WGS84] - Datum to convert grid reference into.
     * @returns {LatLon}       Latitude/longitude of supplied grid reference.
     *
     * @example
     *   const gridref = new OsGridRef(651409.903, 313177.270);
     *   const pWgs84 = gridref.toLatLon();                    // 52°39′28.723″N, 001°42′57.787″E
     *   // to obtain (historical) OSGB36 lat/lon point:
     *   const pOsgb = gridref.toLatLon(LatLon.datums.OSGB36); // 52°39′27.253″N, 001°43′04.518″E
     */
    toLatLon(datum=LatLonEllipsoidal.datums.WGS84) {
        const { easting: E, northing: N } = this;

        const a = 6377563.396, b = 6356256.909;             // Airy 1830 major & minor semi-axes
        const F0 = 0.9996012717;                            // NatGrid scale factor on central meridian
        const phi0 = (49).toRadians(), lambda0 = (-2).toRadians(); // NatGrid true origin is 49°N,2°W
        const N0 = -100e3, E0 = 400e3;                      // northing & easting of true origin, metres
        const e2 = 1 - (b*b)/(a*a);                         // eccentricity squared
        const n = (a-b)/(a+b), n2 = n*n, n3 = n*n*n;        // n, n², n³

        let phi=phi0, M=0;
        do {
            phi = (N-N0-M)/(a*F0) + phi;

            const Ma = (1 + n + (5/4)*n2 + (5/4)*n3) * (phi-phi0);
            const Mb = (3*n + 3*n*n + (21/8)*n3) * Math.sin(phi-phi0) * Math.cos(phi+phi0);
            const Mc = ((15/8)*n2 + (15/8)*n3) * Math.sin(2*(phi-phi0)) * Math.cos(2*(phi+phi0));
            const Md = (35/24)*n3 * Math.sin(3*(phi-phi0)) * Math.cos(3*(phi+phi0));
            M = b * F0 * (Ma - Mb + Mc - Md);               // meridional arc

        } while (Math.abs(N-N0-M) >= 0.00001);  // ie until < 0.01mm

        const cosphi = Math.cos(phi), sinphi = Math.sin(phi);
        const ny = a*F0/Math.sqrt(1-e2*sinphi*sinphi);             // nu = transverse radius of curvature
        const rho = a*F0*(1-e2)/Math.pow(1-e2*sinphi*sinphi, 1.5);     // rho = meridional radius of curvature
        const eta2 = ny/rho-1;                                   // eta = ?

        const tanphi = Math.tan(phi);
        const tan2phi = tanphi*tanphi, tan4phi = tan2phi*tan2phi, tan6phi = tan4phi*tan2phi;
        const secphi = 1/cosphi;
        const ny3 = ny*ny*ny, ny5 = ny3*ny*ny, ny7 = ny5*ny*ny;
        const VII = tanphi/(2*rho*ny);
        const VIII = tanphi/(24*rho*ny3)*(5+3*tan2phi+eta2-9*tan2phi*eta2);
        const IX = tanphi/(720*rho*ny5)*(61+90*tan2phi+45*tan4phi);
        const X = secphi/ny;
        const XI = secphi/(6*ny3)*(ny/rho+2*tan2phi);
        const XII = secphi/(120*ny5)*(5+28*tan2phi+24*tan4phi);
        const XIIA = secphi/(5040*ny7)*(61+662*tan2phi+1320*tan4phi+720*tan6phi);

        const dE = (E-E0), dE2 = dE*dE, dE3 = dE2*dE, dE4 = dE2*dE2, dE5 = dE3*dE2, dE6 = dE4*dE2, dE7 = dE5*dE2;
        phi = phi - VII*dE2 + VIII*dE4 - IX*dE6;
        const lambda = lambda0 + X*dE - XI*dE3 + XII*dE5 - XIIA*dE7;

        let point = new LatLon_OsGridRef(phi.toDegrees(), lambda.toDegrees(), 0, LatLonEllipsoidal.datums.OSGB36);

        if (datum != LatLonEllipsoidal.datums.OSGB36) {
            // if point is required in datum other than OSGB36, convert it
            point = point.convertDatum(datum);
            // convertDatum() gives us a LatLon: convert to LatLon_OsGridRef which includes toOsGrid()
            point = new LatLon_OsGridRef(point.lat, point.lon, point.height, point.datum);
        }

        return point;
    }


    /**
     * Parses grid reference to OsGridRef object.
     *
     * Accepts standard grid references (eg 'SU 387 148'), with or without whitespace separators, from
     * two-digit references up to 10-digit references (1m × 1m square), or fully numeric comma-separated
     * references in metres (eg '438700,114800').
     *
     * @param   {string}    gridref - Standard format OS grid reference.
     * @returns {OsGridRef} Numeric version of grid reference in metres from false origin (SW corner of
     *   supplied grid square).
     * @throws  {Error}     Invalid grid reference.
     *
     * @example
     *   const grid = OsGridRef.parse('TG 51409 13177'); // grid: { easting: 651409, northing: 313177 }
     */
    static parse(gridref) {
        gridref = String(gridref).trim();

        // check for fully numeric comma-separated gridref format
        let match = gridref.match(/^(\d+),\s*(\d+)$/);
        if (match) return new OsGridRef(match[1], match[2]);

        // validate format
        match = gridref.match(/^[A-Z]{2}\s*[0-9]+\s*[0-9]+$/i);
        if (!match) throw new Error(`invalid grid reference ‘${gridref}’`);

        // get numeric values of letter references, mapping A->0, B->1, C->2, etc:
        let l1 = gridref.toUpperCase().charCodeAt(0) - 'A'.charCodeAt(0);
        let l2 = gridref.toUpperCase().charCodeAt(1) - 'A'.charCodeAt(0);
        // shuffle down letters after 'I' since 'I' is not used in grid:
        if (l1 > 7) l1--;
        if (l2 > 7) l2--;

        // sanity check
        if (l1<8 || l1 > 18) throw new Error(`invalid grid reference ‘${gridref}’`);

        // convert grid letters into 100km-square indexes from false origin (grid square SV):
        const e100km = ((l1 - 2) % 5) * 5 + (l2 % 5);
        const n100km = (19 - Math.floor(l1 / 5) * 5) - Math.floor(l2 / 5);

        // skip grid letters to get numeric (easting/northing) part of ref
        let en = gridref.slice(2).trim().split(/\s+/);
        // if e/n not whitespace separated, split half way
        if (en.length == 1) en = [ en[0].slice(0, en[0].length / 2), en[0].slice(en[0].length / 2) ];

        // validation
        if (en[0].length != en[1].length) throw new Error(`invalid grid reference ‘${gridref}’`);

        // standardise to 10-digit refs (metres)
        en[0] = en[0].padEnd(5, '0');
        en[1] = en[1].padEnd(5, '0');

        const e = e100km + en[0];
        const n = n100km + en[1];

        return new OsGridRef(e, n);
    }


    /**
     * Converts ‘this’ numeric grid reference to standard OS grid reference.
     *
     * @param   {number} [digits=10] - Precision of returned grid reference (10 digits = metres);
     *   digits=0 will return grid reference in numeric format.
     * @returns {string} This grid reference in standard format.
     *
     * @example
     *   const gridref = new OsGridRef(651409, 313177).toString(8); // 'TG 5140 1317'
     *   const gridref = new OsGridRef(651409, 313177).toString(0); // '651409,313177'
     */
    toString(digits=10) {
        if (![ 0,2,4,6,8,10,12,14,16 ].includes(Number(digits))) throw new RangeError(`invalid precision ‘${digits}’`); // eslint-disable-line comma-spacing

        let { easting: e, northing: n } = this;

        // use digits = 0 to return numeric format (in metres) - note northing may be >= 1e7
        if (digits == 0) {
            const format = { useGrouping: false,  minimumIntegerDigits: 6, maximumFractionDigits: 3 };
            const ePad = e.toLocaleString('en', format);
            const nPad = n.toLocaleString('en', format);
            return `${ePad},${nPad}`;
        }

        // get the 100km-grid indices
        const e100km = Math.floor(e / 100000), n100km = Math.floor(n / 100000);

        // translate those into numeric equivalents of the grid letters
        let l1 = (19 - n100km) - (19 - n100km) % 5 + Math.floor((e100km + 10) / 5);
        let l2 = (19 - n100km) * 5 % 25 + e100km % 5;

        // compensate for skipped 'I' and build grid letter-pairs
        if (l1 > 7) l1++;
        if (l2 > 7) l2++;
        const letterPair = String.fromCharCode(l1 + 'A'.charCodeAt(0), l2 + 'A'.charCodeAt(0));

        // strip 100km-grid indices from easting & northing, and reduce precision
        e = Math.floor((e % 100000) / Math.pow(10, 5 - digits / 2));
        n = Math.floor((n % 100000) / Math.pow(10, 5 - digits / 2));

        // pad eastings & northings with leading zeros
        e = e.toString().padStart(digits/2, '0');
        n = n.toString().padStart(digits/2, '0');

        return `${letterPair} ${e} ${n}`;
    }

}


/* LatLon_OsGridRef - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */


/**
 * Extends LatLon class with method to convert LatLon point to OS grid reference.
 *
 * @extends LatLon
 */
class LatLon_OsGridRef extends LatLonEllipsoidal {

    /**
     * Converts latitude/longitude to Ordnance Survey grid reference easting/northing coordinate.
     *
     * @returns {OsGridRef} OS Grid Reference easting/northing.
     *
     * @example
     *   const grid = new LatLon(52.65798, 1.71605).toOsGrid(); // TG 51409 13177
     *   // for conversion of (historical) OSGB36 latitude/longitude point:
     *   const grid = new LatLon(52.65798, 1.71605).toOsGrid(LatLon.datums.OSGB36);
     */
    toOsGrid() {
        // if necessary convert to OSGB36 prime
        const point = this.datum == LatLonEllipsoidal.datums.OSGB36
            ? this
            : this.convertDatum(LatLonEllipsoidal.datums.OSGB36);

        const phi = point.lat.toRadians();
        const lambda = point.lon.toRadians();

        const a = 6377563.396, b = 6356256.909;              // Airy 1830 major & minor semi-axes
        const F0 = 0.9996012717;                             // NatGrid scale factor on central meridian
        const phi0 = (49).toRadians(), lambda0 = (-2).toRadians();  // NatGrid true origin is 49°N,2°W
        const N0 = -100000, E0 = 400000;                     // northing & easting of true origin, metres
        const e2 = 1 - (b*b)/(a*a);                          // eccentricity squared
        const n = (a-b)/(a+b), n2 = n*n, n3 = n*n*n;         // n, n², n³

        const cosphi = Math.cos(phi), sinphi = Math.sin(phi);
        const ny = a*F0/Math.sqrt(1-e2*sinphi*sinphi);            // nu = transverse radius of curvature
        const rho = a*F0*(1-e2)/Math.pow(1-e2*sinphi*sinphi, 1.5); // rho = meridional radius of curvature
        const eta2 = ny/rho-1;                                    // eta = ?

        const Ma = (1 + n + (5/4)*n2 + (5/4)*n3) * (phi-phi0);
        const Mb = (3*n + 3*n*n + (21/8)*n3) * Math.sin(phi-phi0) * Math.cos(phi+phi0);
        const Mc = ((15/8)*n2 + (15/8)*n3) * Math.sin(2*(phi-phi0)) * Math.cos(2*(phi+phi0));
        const Md = (35/24)*n3 * Math.sin(3*(phi-phi0)) * Math.cos(3*(phi+phi0));
        const M = b * F0 * (Ma - Mb + Mc - Md);              // meridional arc

        const cos3phi = cosphi*cosphi*cosphi;
        const cos5phi = cos3phi*cosphi*cosphi;
        const tan2phi = Math.tan(phi)*Math.tan(phi);
        const tan4phi = tan2phi*tan2phi;

        const I = M + N0;
        const II = (ny/2)*sinphi*cosphi;
        const III = (ny/24)*sinphi*cos3phi*(5-tan2phi+9*eta2);
        const IIIA = (ny/720)*sinphi*cos5phi*(61-58*tan2phi+tan4phi);
        const IV = ny*cosphi;
        const V = (ny/6)*cos3phi*(ny/rho-tan2phi);
        const VI = (ny/120) * cos5phi * (5 - 18*tan2phi + tan4phi + 14*eta2 - 58*tan2phi*eta2);

        const deltalambda = lambda-lambda0;
        const deltalambda2 = deltalambda*deltalambda, deltalambda3 = deltalambda2*deltalambda, deltalambda4 = deltalambda3*deltalambda, deltalambda5 = deltalambda4*deltalambda, deltalambda6 = deltalambda5*deltalambda;

        let N = I + II*deltalambda2 + III*deltalambda4 + IIIA*deltalambda6;
        let E = E0 + IV*deltalambda + V*deltalambda3 + VI*deltalambda5;

        N = Number(N.toFixed(3)); // round to mm precision
        E = Number(E.toFixed(3));

        try {
            return new OsGridRef(E, N); // note: gets truncated to SW corner of 1m grid square
        } catch (e) {
            throw new Error(`${e.message} from (${point.lat.toFixed(6)},${point.lon.toFixed(6)}).toOsGrid()`);
        }
    }

}


/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */

export { OsGridRef as default, LatLon_OsGridRef as LatLon, Dms };
