// Starmap.cpp
// rev 1 - May 2024 - shabaz - first version

#include "Starmap.h"
#include <stdio.h>
#include <string.h>


// defines




extern mapwindow skywin, telwin;		// Sky and Telescope map windows
extern mapwindow *mapwin[2];			// Window list
extern int numwins;						// Number of windows

//char *starCatResName[] = {		// Star catalogue resource names (don't translate these)
//	  "Yale"                              // 0: Yale bright star catalogue
//  };
//char *tzName[2] = { "", "" };




// *************************************************
// ********  Starmap class implementation  ********
// *************************************************

Starmap::Starmap()
{
  // Starmap object constructor
  mDummy = 0;
  plan_epli = 0.0;
  plan_Stheta = 0.0;
  plan_Sr = 0.0;
  skyLimitMag = 5.5;
  skyShowName = TRUE;
  skyNameMag = 2;
  skyShowBflam = FALSE;
  skyBflamMag = 3.5;
  skyShowDeep = TRUE;
  skyDeepMag = 4.5;
  skyShowConstellations = TRUE;
  skyShowConbounds = TRUE;
  skyShowConnames = TRUE;
  skyAlignConnames = FALSE;
  skyShowCoords = TRUE;
  skyShowPlanets = TRUE;
  skyShowTiming = TRUE;
  precessionCalculation = PrecAuto;
  firstTick = FALSE;
  nextTick = FALSE;
  bailedOut = FALSE;
  calculationNumber = 0;
  Flip = 1;
  pirefc = 0;
  starCatalogue = 0;
  starQuality = 2;
  showStarColours = FALSE;
  jdtime = 0.0;
  siteLat = 47;
  siteLon = 122;
  plutoPrecise = TRUE;
  draw_col_store = 0xffff;
  // default colors
  col_coord_grid = DEFAULT_COL_COORD_GRID;
  col_ecliptic = DEFAULT_COL_ECLIPTIC;
  col_constel = DEFAULT_COL_CONSTEL;
  col_stardim = DEFAULT_COL_STARDIM;
  col_starbright = DEFAULT_COL_STARBRIGHT;
  col_startext = DEFAULT_COL_STARTEXT;

  col_moon_bright = DEFAULT_COL_MOON_BRIGHT;
  col_moon_dim = DEFAULT_COL_MOON_DIM;
  col_moon_dark = DEFAULT_COL_MOON_DARK;
  col_moon_phtext = DEFAULT_COL_MOON_PHTEXT;

  do_constellation_text = 1;

  col_constel_text = DEFAULT_COL_TEXT_GENERIC;
  col_bright_st_text = DEFAULT_COL_TEXT_GENERIC;
  col_bayerf_text = DEFAULT_COL_TEXT_GENERIC;
  col_celest_eq_text = DEFAULT_COL_TEXT_GENERIC;
  col_ecliptic_text = DEFAULT_COL_TEXT_GENERIC;
  

}

double Starmap::jtime(tm_t *t)
{
  double j;
  //	if (t->tm_year==104)
  //	  sprintf(log2ram_buf, "222\n");
  //	while(1);
  j = ucttoj((t->tm_year) + 1900, t->tm_mon, t->tm_mday, t->tm_hour, t->tm_min, t->tm_sec);
  return j;
}

double Starmap::ucttoj(long year, int mon, int mday,
	      int hour, int min, int sec)
{

	// Algorithm as given in Meeus, Astronomical Algorithms, Chapter 7, page 61

	int a, b, m;
	long y;
	double j;

    // assert(mon  >= 0 && mon  < 12);
    // assert(mday >  0 && mday < 32);
    // assert(hour >= 0 && hour < 24);
    // assert(min  >= 0 && min  < 60);
    // assert(sec  >= 0 && sec  < 60);

    m = mon + 1;
    y = year;

	if (m <= 2) {
		y--;
		m += 12;
	}

	/* Determine whether date is in Julian or Gregorian calendar based on
	   canonical date of calendar reform. */

	if ((year < 1582) || ((year == 1582) && ((mon < 9) || (mon == 9 && mday < 5)))) {
		b = 0;
	} else {
		a = ((int) (y / 100));
		b = 2 - a + (a / 4);
	}

	j = (((long) (365.25 * (y + 4716))) + ((int) (30.6001 * (m + 1))) +
				mday + b - 1524.5) +
			((sec + 60L * (min + 60L * hour)) / 86400.0);
	return j;

}

void Starmap::jyear(double td, long *yy, int *mm, int *dd)
{
	double z, f, a, alpha, b, c, d, e;

	td += 0.5;
	z = floor(td);
	f = td - z;

	if (z < 2299161.0) {
		a = z;
	} else {
		alpha = floor((z - 1867216.25) / 36524.25);
		a = z + 1 + alpha - floor(alpha / 4);
	}

	b = a + 1524;
	c = floor((b - 122.1) / 365.25);
	d = floor(365.25 * c);
	e = floor((b - d) / 30.6001);

	*dd = (int) (b - d - floor(30.6001 * e) + f);
	*mm = (int) ((e < 14) ? (e - 1) : (e - 13));
	*yy = (long) ((*mm > 2) ? (c - 4716) : (c - 4715));
}

void Starmap::jhms(double j, int *h, int *m, int *s)
{
    long ij;

    j += 0.5;			      /* Astronomical to civil */
    ij = (long) ((j - floor(j)) * 86400.0);
    *h = (int) (ij / 3600L);
    *m = (int) ((ij / 60L) % 60L);
    *s = (int) (ij % 60L);
}

#ifdef REDUNDANT
double Starmap::kepler(double m, double ecc)
{
    double e, delta;
#define EPSILON 1E-6

    e = m = dtr(m);
    do {
	delta = e - ecc * sin(e) - m;
	e -= delta / (1 - ecc * cos(e));
    } while (abs(delta) > EPSILON);
    return e;
}
#endif

double Starmap::gmst(double jd)
{
    double t, theta0;

    /* Time, in Julian centuries of 36525 ephemeris days,
       measured from the epoch 1900 January 0.5 ET. */
mydeb=0xb0;
    t = floor(jd + 0.5) - 0.5;
    t=t - 2415020.0;
    
    t=t / JulianCentury;
mydeb=0xb5;
    theta0 = 6.6460656 + 2400.051262 * t;
    theta0 = theta0 + 0.00002581 * t * t;
mydeb=0xb9;
    t = jd + 0.5;
    t=t - floor(jd + 0.5);
mydeb=0xbc;
    theta0 += (t * 24.0) * 1.002737908;
mydeb=0xbe;
    theta0 = (theta0 - 24.0 * (floor(theta0 / 24.0)));
mydeb=0xbf;
    return theta0;
}

void Starmap::rgmst(double jd, double* ojd)
{
    double t, theta0;

    /* Time, in Julian centuries of 36525 ephemeris days,
       measured from the epoch 1900 January 0.5 ET. */
mydeb=0xb0;
    t = floor(jd + 0.5) - 0.5;
    t=t - 2415020.0;
    
    t=t / JulianCentury;
mydeb=0xb5;
    theta0 = 6.6460656 + 2400.051262 * t;
    theta0 = theta0 + 0.00002581 * t * t;
mydeb=0xb9;
    t = jd + 0.5;
    t=t - floor(jd + 0.5);
mydeb=0xbc;
    theta0 += (t * 24.0) * 1.002737908;
mydeb=0xbe;
    theta0 = (theta0 - 24.0 * (floor(theta0 / 24.0)));
mydeb=0xbf;
   // return theta0;
   *ojd=theta0;
}

double Starmap::phase(
  double  pdate,		      /* Date for which to calculate phase */
  double  *pphase,		  /* Illuminated fraction */
  double  *mage,		  /* Age of moon in days */
  double  *dist,		  /* Distance in kilometres */
  double  *angdia,		  /* Angular diameter in degrees */
  double  *sudist,		  /* Distance to Sun */
  double  *suangdia)      /* Sun's angular diameter */
{

    double Day, N, M, Ec, Lambdasun, ml, MM, MN, Ev, Ae, A3, MmP,
	   mEc, A4, lP, Varia, lPP, NP, y, x, Lambdamoon,
	   MoonAge, MoonPhase,
	   MoonDist, MoonDFrac, MoonAng,
	   F, SunDist, SunAng;

    /* Calculation of the Sun's position */

    Day = pdate - epoch;		    /* Date within epoch */
    N = fixangle((360 / 365.2422) * Day);   /* Mean anomaly of the Sun */
    M = fixangle(N + elonge - elongp);	    /* Convert from perigee
					       co-ordinates to epoch 1980.0 */
    Ec = kepler(M, eccent);		    /* Solve equation of Kepler */
    Ec = sqrt((1 + eccent) / (1 - eccent)) * tan(Ec / 2);
    Ec = 2 * rtd(atan(Ec));		  /* True anomaly */
    Lambdasun = fixangle(Ec + elongp);      /* Sun's geocentric ecliptic
					       longitude */
    /* Orbital distance factor */
    F = ((1 + eccent * cos(dtr(Ec))) / (1 - eccent * eccent));
    SunDist = SunSMAX / F;		    /* Distance to Sun in km */
    SunAng = F * sunangsiz;                 /* Sun's angular size in degrees */

    /* Calculation of the Moon's position */

    /* Moon's mean longitude */
    ml = fixangle(13.1763966 * Day + mmlong);

    /* Moon's mean anomaly */
    MM = fixangle(ml - 0.1114041 * Day - mmlongp);

    /* Moon's ascending node mean longitude */
    MN = fixangle(mlnode - 0.0529539 * Day);

    /* Evection */
    Ev = 1.2739 * sin(dtr(2 * (ml - Lambdasun) - MM));

    /* Annual equation */
    Ae = 0.1858 * sin(dtr(M));

    /* Correction term */
    A3 = 0.37 * sin(dtr(M));

    /* Corrected anomaly */
    MmP = MM + Ev - Ae - A3;

    /* Correction for the equation of the centre */
    mEc = 6.2886 * sin(dtr(MmP));

    /* Another correction term */
    A4 = 0.214 * sin(dtr(2 * MmP));

    /* Corrected longitude */
    lP = ml + Ev + mEc - Ae + A4;

    /* Variation */
    Varia = 0.6583 * sin(dtr(2 * (lP - Lambdasun)));

    /* True longitude */
    lPP = lP + Varia;

    /* Corrected longitude of the node */
    NP = MN - 0.16 * sin(dtr(M));

    /* Y inclination coordinate */
    y = sin(dtr(lPP - NP)) * cos(dtr(minc));

    /* X inclination coordinate */
    x = cos(dtr(lPP - NP));

    /* Ecliptic longitude */
    Lambdamoon = rtd(atan2(y, x));
    Lambdamoon += NP;

    /* Ecliptic latitude */
//    BetaM = rtd(asin(sin(dtr(lPP - NP)) * sin(dtr(minc))));

    /* Calculation of the phase of the Moon */

    /* Age of the Moon in degrees */
    MoonAge = lPP - Lambdasun;

    /* Phase of the Moon */
    MoonPhase = (1 - cos(dtr(MoonAge))) / 2;

    /* Calculate distance of moon from the centre of the Earth */

    MoonDist = (msmax * (1 - mecc * mecc)) /
	       (1 + mecc * cos(dtr(MmP + mEc)));

    /* Calculate Moon's angular diameter */

    MoonDFrac = MoonDist / msmax;
    MoonAng = mangsiz / MoonDFrac;

    /* Calculate Moon's parallax */

//    MoonPar = mparallax / MoonDFrac;

    *pphase = MoonPhase;
    *mage = SynMonth * (fixangle(MoonAge) / 360.0);
    *dist = MoonDist;
    *angdia = MoonAng;
    *sudist = SunDist;
    *suangdia = SunAng;
    return fixangle(MoonAge) / 360.0;
}

void Starmap::highmoon(double jd, double *Lambda, double *Beta)
{
	double t, t2, t3, lprime, m, mprime, d, f, omega, v, e, e2,
	       b, om1, om2;

	/* Time in Julian centuries from 1900 January 0.5 */
	t = (jd - 2415020.0) / JulianCentury;
	t2 = t * t;		   /* Square for frequent use */
	t3 = t2 * t;		   /* Cube for frequent use */

        /* Moon's mean longitude */
	lprime = 270.434164 + 481267.8831 * t - 0.001133 * t2 +
		 0.0000019 * t3;

        /* Sun's mean anomaly */
	m = 358.475833 + 35999.0498 * t - 0.000150 * t2 -
	    0.0000033 * t3;

        /* Moon's mean anomaly */
	mprime = 296.104608 + 477198.8491 * t + 0.009192 * t2 +
		 0.0000144 * t3;

        /* Moon's mean elongation */
	d = 350.737486 + 445267.1142 * t - 0.001436 * t2 +
	    0.0000019 * t3;

	/* Mean distance of the Moon from its ascending node */
	f = 11.250889 + 483202.0251 * t - 0.003211 * t2 -
	    0.0000003 * t3;

        /* Longitude of the Moon's ascending node */
	omega = 259.183275 - 1934.1420 * t + 0.002078 * t2 +
		0.0000022 * t3;

	/* Explicitly range-reduce base angles to reduce additive
	   round-off in subsequent calculations. */

	lprime = fixangle(lprime);
	m = fixangle(m);
	mprime = fixangle(mprime);
	d = fixangle(d);
	f = fixangle(f);
	omega = fixangle(omega);

	/* Additive terms:
		1782-year periodic variations  */
	v = sin(dtr(51.2 + 20.2 * t));
	lprime +=   0.000233 * v;
	m      +=  -0.001778 * v;
	mprime +=   0.000817 * v;
	d      +=   0.002011 * v;

	/*	 The Great Venus Term: 271 year period */
	v = 0.003964 *
	    sin(dtr(346.560 + 132.870 * t - 0.0091731 * t2));
	lprime += v;
	mprime += v;
	d += v;
	f += v;

	/*	 Corrections based on ascending node longitude */
	v = sin(dtr(omega));
	lprime +=  0.001964 * v;
	mprime +=  0.002541 * v;
	d      +=  0.001964 * v;
	f      += -0.024691 * v;
	f      += -0.004328 * sin(dtr(omega + 275.05 - 2.30 * t));

	/* Common multiplicative factors used in the expansions. */

	e = 1 - 0.002495 * t - 0.00000752 * t2;
	e2 = e * e;

	/* Geocentric longitude of the Moon */

	*Lambda = lprime
		    +	   6.288750 * dsin(mprime)
		    +	   1.274018 * dsin(2 * d - mprime)
		    +	   0.658309 * dsin(2 * d)
		    +	   0.213616 * dsin(2 * mprime)
		    - e *  0.185596 * dsin(m)
		    -	   0.114336 * dsin(2 * f)
		    +	   0.058793 * dsin(2 * d - 2 * mprime)
		    + e *  0.057212 * dsin(2 * d - m - mprime)
		    +	   0.053320 * dsin(2 * d + mprime)
		    + e *  0.045874 * dsin(2 * d - m)
		    + e *  0.041024 * dsin(mprime - m)
		    -	   0.034718 * dsin(d)
		    - e *  0.030465 * dsin(m + mprime)
		    +	   0.015326 * dsin(2 * d - 2 * f)
		    -	   0.012528 * dsin(2 * f + mprime)
		    -	   0.010980 * dsin(2 * f - mprime)
		    +	   0.010674 * dsin(4 * d - mprime)
		    +	   0.010034 * dsin(3 * mprime)
		    +	   0.008548 * dsin(4 * d - 2 * mprime)
		    - e *  0.007910 * dsin(m - mprime + 2 * d)
		    - e *  0.006783 * dsin(2 * d + m)
		    +	   0.005162 * dsin(mprime - d)
		    + e *  0.005000 * dsin(m + d)
		    + e *  0.004049 * dsin(mprime - m + 2 * d)
		    +	   0.003996 * dsin(2 * mprime + 2 * d)
		    +	   0.003862 * dsin(4 * d)
		    +	   0.003665 * dsin(2 * d - 3 * mprime)
		    + e *  0.002695 * dsin(2 * mprime - m)
		    +	   0.002602 * dsin(mprime - 2 * f - 2 * d)
		    + e *  0.002396 * dsin(2 * d - m - 2 * mprime)
		    -	   0.002349 * dsin(mprime + d)
		    + e2 * 0.002249 * dsin(2 * d - 2 * m)
		    - e *  0.002125 * dsin(2 * mprime + m)
		    - e2 * 0.002079 * dsin(2 * m)
		    + e2 * 0.002059 * dsin(2 * d - mprime - 2 * m)
		    -	   0.001773 * dsin(mprime + 2 * d - 2 * f)
		    -	   0.001595 * dsin(2 * f + 2 * d)
		    + e *  0.001220 * dsin(4 * d - m - mprime)
		    -	   0.001110 * dsin(2 * mprime + 2 * f)
		    +	   0.000892 * dsin(mprime - 3 * d)
		    - e *  0.000811 * dsin(m + mprime + 2 * d)
		    + e *  0.000761 * dsin(4 * d - m - 2 * mprime)
		    + e2 * 0.000717 * dsin(mprime - 2 * m)
		    + e2 * 0.000704 * dsin(mprime - 2 * m - 2 * d)
		    + e *  0.000693 * dsin(m - 2 * mprime + 2 * d)
		    + e *  0.000598 * dsin(2 * d - m - 2 * f)
		    +	   0.000550 * dsin(mprime + 4 * d)
		    +	   0.000538 * dsin(4 * mprime)
		    + e *  0.000521 * dsin(4 * d - m)
		    +	   0.000486 * dsin(2 * mprime - d);

	/* Geocentric latitude of the Moon */

	b =	     5.128189 * dsin(f)
	      +      0.280606 * dsin(mprime + f)
	      +      0.277693 * dsin(mprime - f)
	      +      0.173238 * dsin(2 * d - f)
	      +      0.055413 * dsin(2 * d + f - mprime)
	      +      0.046272 * dsin(2 * d - f - mprime)
	      +      0.032573 * dsin(2 * d + f)
	      +      0.017198 * dsin(2 * mprime + f)
	      +      0.009267 * dsin(2 * d + mprime - f)
	      +      0.008823 * dsin(2 * mprime - f)
	      + e *  0.008247 * dsin(2 * d - m - f)
	      +      0.004323 * dsin(2 * d - f - 2 * mprime)
	      +      0.004200 * dsin(2 * d + f + mprime)
	      + e *  0.003372 * dsin(f - m - 2 * d)
	      + e *  0.002472 * dsin(2 * d + f - m - mprime)
	      + e *  0.002222 * dsin(2 * d + f - m)
	      + e *  0.002072 * dsin(2 * d - f - m - mprime)
	      + e *  0.001877 * dsin(f - m + mprime)
	      +      0.001828 * dsin(4 * d - f - mprime)
	      - e *  0.001803 * dsin(f + m)
	      -      0.001750 * dsin(3 * f)
	      + e *  0.001570 * dsin(mprime - m - f)
	      -      0.001487 * dsin(f + d)
	      - e *  0.001481 * dsin(f + m + mprime)
	      + e *  0.001417 * dsin(f - m - mprime)
	      + e *  0.001350 * dsin(f - m)
	      +      0.001330 * dsin(f - d)
	      +      0.001106 * dsin(f + 3 * mprime)
	      +      0.001020 * dsin(4 * d - f)
	      +      0.000833 * dsin(f + 4 * d - mprime)
	      +      0.000781 * dsin(mprime - 3 * f)
	      +      0.000670 * dsin(f + 4 * d - 2 * mprime)
	      +      0.000606 * dsin(2 * d - 3 * f)
	      +      0.000597 * dsin(2 * d + 2 * mprime - f)
	      + e *  0.000492 * dsin(2 * d + mprime - m - f)
	      +      0.000450 * dsin(2 * mprime - f - 2 * d)
	      +      0.000439 * dsin(3 * mprime - f)
	      +      0.000423 * dsin(f + 2 * d + 2 * mprime)
	      +      0.000422 * dsin(2 * d - f - 3 * mprime)
	      - e *  0.000367 * dsin(m + f + 2 * d - mprime)
	      - e *  0.000353 * dsin(m + f + 2 * d)
	      +      0.000331 * dsin(f + 4 * d)
	      + e *  0.000317 * dsin(2 * d + f - m + mprime)
	      + e2 * 0.000306 * dsin(2 * d - 2 * m - f)
	      -      0.000283 * dsin(mprime + 3 * f);

	om1 = 0.0004664 * dcos(omega);
	om2 = 0.0000754 * dcos(omega + 275.05 - 2.30 * t);
	*Beta = b * (1 - om1 - om2);
}

double Starmap::obliqeq(double jd)
{
#define Asec(x)	((x) / 3600.0)

	static double oterms[10] = {
		Asec(-4680.93),
		Asec(   -1.55),
		Asec( 1999.25),
		Asec(  -51.38),
		Asec( -249.67),
		Asec(  -39.05),
		Asec(    7.12),
		Asec(   27.87),
		Asec(    5.79),
		Asec(    2.45)
	};

	double eps = 23 + (26 / 60.0) + (21.448 / 3600.0), u, v;
	int i;

	v = u = (jd - J2000) / (JulianCentury * 100);

    if (abs(u) < 1.0) {
		for (i = 0; i < 10; i++) {
			eps += oterms[i] * v;
			v *= u;
		}
	}
	return eps;
}

void Starmap::robliqeq(double jd, double *rdj)
{
#define Asec(x)	((x) / 3600.0)

	static double oterms[10] = {
		Asec(-4680.93),
		Asec(   -1.55),
		Asec( 1999.25),
		Asec(  -51.38),
		Asec( -249.67),
		Asec(  -39.05),
		Asec(    7.12),
		Asec(   27.87),
		Asec(    5.79),
		Asec(    2.45)
	};

	double eps = 23 + (26 / 60.0) + (21.448 / 3600.0), u, v;
	int i;

	v = u = (jd - J2000) / (JulianCentury * 100);

    if (abs(u) < 1.0) {
		for (i = 0; i < 10; i++) {
			eps += oterms[i] * v;
			v *= u;
		}
	}
	*rdj=eps;
	//return eps;
}

void Starmap::ecliptoeq(double jd, double Lambda, double Beta,
					   double *Ra, double *Dec)
{
	double eps;

    /* Obliquity of the ecliptic. */

    
    robliqeq(jd, &eps);
    eps=dtr(eps);
    //eps = dtr(obliqeq(jd));

    *Ra = fixangle(rtd(atan2((cos(eps) * sin(dtr(Lambda)) -
					     (tan(dtr(Beta)) * sin(eps))), cos(dtr(Lambda)))));
    *Dec = rtd(asin((sin(eps) * sin(dtr(Lambda)) * cos(dtr(Beta))) +
			     (sin(dtr(Beta)) * cos(eps))));

}

void Starmap::definePrecession(double targetEpoch)
{
	double t, t2, t3;

	t = (targetEpoch - 2000.0) / 100.0;
	t3 = (t2 = t * t) * t;
#define SecToR(x)	(dtr((x)) / 3600.0)	// Seconds of arc to radians
	preZeta = SecToR(2306.2181 * t + 0.30188 * t2 + 0.017998 * t3);
	preZ = SecToR(2306.2181 * t + 1.09468 * t2 + 0.018203 * t3);
	preTheta = SecToR(2004.3109 * t - 0.42665 * t2 - 0.041833 * t3);
}

void Starmap::precessObject(double ira, double idec, double *ora, double *odec)
{
	double rira = dtr(ira), ridec = dtr(idec), a, b, c;

	a = cos(ridec) * sin(rira + preZeta);
	b = cos(preTheta) * cos(ridec) * cos(rira + preZeta) - sin(preTheta) * sin(ridec);
	c = sin(preTheta) * cos(ridec) * cos(rira + preZeta) + cos(preTheta) * sin(ridec);

	*ora = rtd(atan2(a, b) + preZ);
	*odec = rtd((idec > 85.0) ? acos(sqrt(a * a + b * b)) : asin(c));
}

int Starmap::initxform(mapwindow *win)
{
	xfs_proj_mode = win->proj_mode;
		
	if (win->scale <= 0.0)
		return FALSE;
	if (win->height == 0)
		return FALSE;
		
	xfs_ra_cen = win->racen;
	xfs_dl_cen = win->dlcen;
	xf_xcen = win->x_offset + (win->width) / 2;
	xf_ycen = win->y_offset + (win->height) / 2;
	xf_ybot = win->y_offset;
		
	xf_north = (win->dlcen + win->scale / 2);
	xf_south = (win->dlcen - win->scale / 2);
		
	if (xf_north > 90.0)
		xf_north = 90.0;
	if (xf_south < -90.0)
		xf_south = -90.0;
		
	if (win->invert) {
		xfs_vinv = -1.0;
		xfs_hinv = 1;
		xf_bottom = xf_north;
	} else {
		xfs_vinv = 1.0;
		xfs_hinv = -1;
		xf_bottom = xf_south;
	}
	
  	xfs_wide_warn = FALSE;


	xf_w_left = win->x_offset;
	xf_w_right = win->x_offset + win->width;
	xf_w_bot = win->y_offset;
	xf_w_top = win->y_offset + win->height;

	switch (xfs_proj_mode) {
		case STEREOGR:
			sin_dlcen = DSIN(win->dlcen);
			cos_dlcen = DCOS(win->dlcen);
			chart_scale = win->scale * .0174532925199; /* Radians */
			break;
		default:
			break;
	}

	/* xf_c_scale is the size in degrees which one pixel occupies on the map */
	/* xfs_scale is the conversion factor for size of the picture
	     (= R in some formulas for stereographic,
	     gnomonic and orthographic projections) */
	     
	xfs_scale = MIN(win->height, win->width) / (4.0 * DTAN(win->scale / 2.0));
	xf_c_scale = win->c_scale = 1.0 / (2.0 * DTAN(0.5) * xfs_scale);
	
	/* initialize gnomonic transform function */
	init_gt(win);
	
	return TRUE;
}

void Starmap::xform(double lat, double lon, int *xloc, int *yloc, int *inregion)
{
  double theta, actheta, rac_l;
  double denom;
  double Dcoslat, Dsinlat, Dcosrac_l, Dsinrac_l;
  /* Dcoslat, Dsinlat: of object latitude in degrees = phi
     Dcosrac_l, Dsinrac_l: of object ra - longditude of center = d(lambda) */

  switch (xfs_proj_mode) {

  case STEREOGR:
	/* Stereographic projection */
    rac_l = lon - xfs_ra_cen;
    Dsinlat = DSIN(lat);
    Dcoslat = DCOS(lat);
    Dcosrac_l = DCOS(rac_l);
    Dsinrac_l = DSIN(rac_l);
    actheta = sin_dlcen * Dsinlat + cos_dlcen * Dcoslat * Dcosrac_l;
    if (actheta > 1.0)
    	theta = 0.0;
    else if (actheta < -1.0)
    	theta = 3.14159265358979323846;
    else theta = acos(actheta);

    *inregion = (theta <= chart_scale);
    if (*inregion) {
      denom = (1 + sin_dlcen * Dsinlat + cos_dlcen * Dcoslat * Dcosrac_l) / xfs_scale;
      *xloc = (int) (xf_xcen - 2 * xfs_hinv * Dcoslat * Dsinrac_l / denom + 0.5);
      *yloc = (int) (xf_ycen + 2 * xfs_vinv * (cos_dlcen * Dsinlat
                                 - sin_dlcen * Dcoslat * Dcosrac_l) / denom + 0.5);
    }
    break;

  default:
    break;
  }
}

void Starmap::init_gt(mapwindow *win)
{
	double adj;

	gt_use_boundaries = TRUE;

	gt_sin_dlcen = DSIN(win->dlcen);
	gt_cos_dlcen = DCOS(win->dlcen);
	gt_chart_scale = win->scale * .0174532925199; /* Radians */

	/* gt_scale is the conversion factor for size of the picture ( = R) */
	gt_scale = MIN(win->height, win->width) / (2.0 * DTAN(win->scale));

	adj = xf_c_scale * 0.9;			/* use boundaries slightly
	                                   more restricted than full plot */

	/* calculate boundaries of region */
	switch (xfs_proj_mode) {
		case STEREOGR:
			gt_r = MIN(win->height, win->width) / 2.0 - 1;
			break;

		default:                      /* error */
			break;
	}
}

/* Note, returns xloc and yloc as doubles */
void Starmap::do_gt(double lat, double lon, double *xloc, double *yloc, double *r_theta)
{
	double theta, rac_l;
	double denom;
	double Dcoslat, Dsinlat, Dcosrac_l, Dsinrac_l;
	/* Dcoslat, Dsinlat: of object latitude in degrees = phi
	 Dcosrac_l, Dsinrac_l: of object ra - longditude of center = d(lambda) */

	rac_l = lon - xfs_ra_cen;
	Dsinlat = DSIN(lat);
	Dcoslat = DCOS(lat);
	Dcosrac_l = DCOS(rac_l);
	Dsinrac_l = DSIN(rac_l);

	*r_theta =
		theta = acos(gt_sin_dlcen*Dsinlat + gt_cos_dlcen * Dcoslat * Dcosrac_l);

	if (theta <= 1.57) { /* avoid wrapping */
		denom = (gt_sin_dlcen * Dsinlat + gt_cos_dlcen * Dcoslat * Dcosrac_l) / gt_scale;
		*yloc = xfs_vinv *
		    (gt_cos_dlcen * Dsinlat - gt_sin_dlcen * Dcoslat * Dcosrac_l) / denom;
		*xloc = xfs_hinv * (- Dcoslat * Dsinrac_l / denom);
	};
}

void Starmap::inv_gt(double x, double y, double *latp, double *lonp)
{
	double rho;
	double R, theta;
	double l, m, n;
	double l_, m_, n_;
    
    x *= xfs_hinv;
	y *= xfs_vinv;

	*latp = 0.0;
	*lonp = 0.0;

	/* Calculate lat and lon */
	R = sqrt((double) ((((long) x) * x) + (((long) y) * y)));
	theta = atan2((double) y, (double) x);

	/* rho is the angle from the center of the display to the object on the
	 unit sphere. */
	rho = atan(R / gt_scale);

	/* transform from (rho, theta) to l m n direction cosines */
	l = sin(rho) * cos(theta);      /* rho and theta are in radians */
	m = sin(rho) * sin(theta);
	n = cos(rho);

	/* transform to new declination at center
	 new axes rotated about x axis (l) */
	l_ = l;
	m_ = m * gt_sin_dlcen - n * gt_cos_dlcen;
	n_ = m * gt_cos_dlcen + n * gt_sin_dlcen;

	/* calculate lon and lat */
	*lonp = atan2(l_, m_) / 0.0174532925199 + xfs_ra_cen - 180.0;
	if (n_ > 1) n_ = 1;
	if (n_ < -1) n_ = -1;
	*latp = 90 - acos(n_) / 0.0174532925199;

	if (*lonp >= 360.0)
		*lonp -= 360.0;
	if (*lonp < 0.0)
		*lonp += 360.0;
}

void Starmap::quadrat(double a, double b, double c, double *x_1, double *x_2, int *n)
{
	double t;

	if (a == 0) {
		*n = 0;
	} else {
		t = b * b - 4 * a * c;
		if (t < 0) {
		  *n = 0;
		} else if (t == 0) {
		  *x_1 = -b/(2*a);
		  *n = 1;
		} else {
		  *x_1 = (-b + sqrt(t)) / (2 * a);
		  *x_2 = (-b - sqrt(t)) / (2 * a);
		  *n = 2;
		};
	};
}

void Starmap::gcmidpoint(double lat1, double lon1, double lat2, double lon2,
					   double *pmlat, double *pmlon)
{
	double l1, m1, n1;
	double l2, m2, n2;
	double l3, m3, n3;

	/* transform from (ra, dec) to l m n direction cosines */
	l1 = DCOS(lat1) * DCOS(lon1);
	m1 = DCOS(lat1) * DSIN(lon1);
	n1 = DSIN(lat1);

	l2 = DCOS(lat2) * DCOS(lon2);
	m2 = DCOS(lat2) * DSIN(lon2);
	n2 = DSIN(lat2);

	l3 = l1 + l2;
	m3 = m1 + m2;
	n3 = n1 + n2;
	n3 /= sqrt(l3 * l3 + m3 * m3 + n3 * n3);

	*pmlon = atan2(m3, l3) / 0.0174532925199;
	if ((*pmlon < 0) && (lon1 > 0) && (lon2 > 0))
		*pmlon += 360.0;
	*pmlat = asin(n3) / 0.0174532925199;
}

void Starmap::circ_intersect(double x_1, double y_1, double x_2, double y_2,
						   double r, double *x1, double *y1, int *int_1,
						   			 double *x2, double *y2, int *int_2)
{
	double c, d;
	double xroot1, xroot2;
	double yr1, yr2, r1, r2;
	int n;
	int xt1, yt1;
	double lat_1, lon_1;
	int in;

	if (fabs(x_2 - x_1) < 1e-5) {         /* Line has infinite slope */
		xroot1 = xroot2 = *x1 = *x2 = x_1;
		if (fabs(r) > fabs(x_1)) {
		  yr1 = sqrt(r * r - x_1 * x_1);
		  yr2 = -yr1;
		  n = 2;
		  *int_1 = *int_2 = TRUE;
		} else if (fabs(r) == fabs(x_1)) {
		  yr1 = 0;
		  n = 1;
		  *int_1 = *int_2 = TRUE;
		} else {
		  n = 0;
		  *int_1 = *int_2 = FALSE;
		};
		} else {                              /* Line slope may be calculated */
		c = (y_2 - y_1)/(x_2 - x_1);
		d = y_1 - c * x_1;
		quadrat(1 + c * c, 2 * c * d, d * d - r * r,
		        &xroot1, &xroot2, &n);
		if (n == 0) {               /* No intersection */
		  *int_1 = *int_2 = FALSE;
		} else if (n == 1) {        /* One intersection */
		  *x1 = xroot1;
		  *y1 = c * xroot1 + d;
		  *int_1 = TRUE;
		  *int_2 = FALSE;
		} else {                    /* Two intersections */
		  yr1 = c*xroot1 + d;
		  yr2 = c*xroot2 + d;
		  *int_1 = *int_2 = TRUE;
		};
	};

	if (n == 2) {
		r1 = (xroot1 - x_1) * (xroot1 - x_1) + (yr1 - y_1) * (yr1 - y_1)
		      + (xroot1 - x_2) * (xroot1 - x_2) + (yr1 - y_2) * (yr1 - y_2);
		r2 = (xroot2 - x_1) * (xroot2 - x_1) + (yr2 - y_1) * (yr2 - y_1)
		      + (xroot2 - x_2) * (xroot2 - x_2) + (yr2 - y_2) * (yr2 - y_2);
		if (r1 > r2) {
		  *x1 = xroot2;
		  *y1 = yr2;
		  *x2 = xroot1;
		  *y2 = yr1;
		  *int_1 = *int_2 = TRUE;
		} else {
		  *x1 = xroot1;
		  *y1 = yr1;
		  *x2 = xroot2;
		  *y2 = yr2;
		  *int_1 = *int_2 = TRUE;
		}
	}

	if (*int_1)
		if ((((y_1 <= *y1) && (*y1 <= y_2)) || ((y_2 <= *y1) && (*y1 <= y_1)))
		    && (((x_1 <= *x1) && (*x1 <= x_2)) || ((x_2 <= *x1) && (*x1 <= x_1))))
		  *int_1 = TRUE;
		else
		  *int_1 = FALSE;

		if (*int_1) {
		inv_gt(*x1, *y1, &lat_1, &lon_1);
		xform(lat_1, lon_1, &xt1, &yt1, &in);
		if (!in)
			*int_1 = FALSE;
	}

	if (*int_2)
		if ((((y_1 <= *y2) && (*y2 <= y_2)) || ((y_2 <= *y2) && (*y2 <= y_1)))
		    && (((x_1 <= *x2) && (*x2 <= x_2)) || ((x_2 <= *x2) && (*x2 <= x_1))))
		  *int_2 = TRUE;
		else
		  *int_2 = FALSE;

		if (*int_2) {
			inv_gt(*x2, *y2, &lat_1, &lon_1);
			xform(lat_1, lon_1, &xt1, &yt1, &in);
			if (!in)
				*int_2 = FALSE;
		}
}

int Starmap::clipr_xform(double lat1, double lon1, double lat2, double lon2,
				int *xloc1, int *yloc1, int *xloc2, int *yloc2,
				int great_circle,
				double *plat1, double *plon1, double *plat2, double *plon2)
{
	int Lisin, Risin;
	double x_1, y_1, x_2, y_2;
	double theta_1, theta_2;
	int int_w, int_e, int_n, int_s, int_r1, int_r2;
	double xr1, xr2, yr1, yr2;

	double Llat, Llon, Rlat, Rlon, Mlat, Mlon;
	int Lx, Ly, Rx, Ry, Mx, My;
	int inL, inR, inM;


	*plon1 = lon1;
	*plon2 = lon2;
	*plat1 = lat1;
	*plat2 = lat2;
	clip_at1 = clip_at2 = NO_CLIP;
	xform(lat1, lon1, xloc1, yloc1, &Lisin);
	xform(lat2, lon2, xloc2, yloc2, &Risin);
	if (Lisin && Risin)           /* is already ok: case 1 */
		return TRUE;

	if (great_circle && gt_use_boundaries) {
	/* Transform to gnomonic */
	do_gt(lat1, lon1, &x_1, &y_1, &theta_1);
	do_gt(lat2, lon2, &x_2, &y_2, &theta_2);

	if ((theta_1 > 1.57) || (theta_2 > 1.57)) /* out of field, skip */
	  return FALSE;

	/* Find intersections with boundaries */
	switch (xfs_proj_mode) {
		case STEREOGR:
		  circ_intersect(x_1, y_1, x_2, y_2, gt_r, &xr1, &yr1, &int_r1,
		                  &xr2, &yr2, &int_r2);
		  int_w = int_n = int_s = int_e = FALSE;

		default:                          /* error */
		  break;
	};


	if (!(!Lisin && !Risin)) {         /* case 2 */
if (int_r1) {
	    x_1 = xr1; y_1 = yr1;
	    if (Risin)
	      clip_at1 = RADIUS_CLIP;
	    else
	      clip_at2 = RADIUS_CLIP;
	  } else {
	/*      fprintf(stderr, "Error drawing vector\n");
	    fprintf(stderr, "from (%.3f %.3f) to (%.3f %.3f)\n",
	            lat1, lon1, lat2, lon2);*/
	    return FALSE;
	  };
	  if (Lisin) {                    /* Need to find new point 2 */
	    inv_gt(x_1, y_1, plat2, plon2);
	    xform(*plat2, *plon2, xloc2, yloc2, &inM);
	  } else {                        /* Need to find new point 1 */
	    inv_gt(x_1, y_1, plat1, plon1);
	    xform(*plat1, *plon1, xloc1, yloc1, &inM);
	  };
	} else {                          /* case 3 */
	  return FALSE;

//	  inv_gt(x_1, y_1, plat1, plon1);
//	  inv_gt(x_2, y_2, plat2, plon2);
//	  xform(*plat1, *plon1, xloc1, yloc1, &inM);
//	  xform(*plat2, *plon2, xloc2, yloc2, &inM);
	}
	return TRUE;
	} else {                            /* find boundaries by bisection */

	if (!Lisin && !Risin)        /* is hopeless */
	  return FALSE;

	/* Now, one side is in, and the other out.  Make sure we won't have
	   problems with crossing 0h */
	/* If the difference between lon1 and lon2 is greater than
	   the difference if you subtract 360 from the larger,
	   then shift the larger by 360 degrees */

	if (fabs(MAX(lon1,lon2) - MIN(lon1,lon2))
	    > fabs(MAX(lon1,lon2) - 360.0 - MIN(lon1,lon2)))
	  if (lon2 > 180.0) lon2 -= 360.0;
	  else lon1 -= 360.0;

	Llat = lat1;
	Llon = lon1;
	Rlat = lat2;
	Rlon = lon2;
	xform(Llat, Llon, &Lx, &Ly, &inL);
	xform(Rlat, Rlon, &Rx, &Ry, &inR);


	/* One endpoint is in.
	   Now use bisection to find point at edge */
	do {
	  if (great_circle) {
	    gcmidpoint(Llat, Llon, Rlat, Rlon, &Mlat, &Mlon);
	  } else {
	    Mlat = (Llat + Rlat) / 2.0;
	    Mlon = (Llon + Rlon) / 2.0;
	  };
	  xform(Mlat, Mlon, &Mx, &My, &inM);

	  if (inL)                   /* L in R out */
	    if (inM) {               /* between M and R */
	      Llat = Mlat;
	      Llon = Mlon;
	      inL = inM;
	      Lx = Mx;
	      Ly = My;
	    } else {                 /* between M and L */
	      Rlat = Mlat;
	      Rlon = Mlon;
	      inR = inM;
	      Rx = Mx;
	      Ry = My;
	    }
	  else                       /* L out R in */
	    if (inM) {               /* between M and L */
	      Rlat = Mlat;
	      Rlon = Mlon;
	      inR = inM;
	      Rx = Mx;
	      Ry = My;
	    } else {                 /* between M and R */
	      Llat = Mlat;
	      Llon = Mlon;
	      inL = inM;
	      Lx = Mx;
	      Ly = My;
	    };
	} while ((fabs((Llat - Rlat)) > xf_c_scale) ||
	         (fabs((Llon - Rlon)) > xf_c_scale));

	if (Lisin) {                 /* Left point is in,
	                                bisection found right point */
	  *xloc2 = Lx;               /* Use Lx, Ly, since they're inside bounds */
	  *yloc2 = Ly;
	  *plon2 = Llon;
	  *plat2 = Llat;
	} else {                     /* Bisection found left point */
	  *xloc1 = Rx;               /* Use Rx, Ry, since they're inside bounds */
	  *yloc1 = Ry;
	  *plon1 = Rlon;
	  *plat1 = Rlat;
	}

	return TRUE;
	}
}

void Starmap::drawcurveline(double  lat1, double lon1, double lat2, double lon2,
				   int xloc1, int yloc1, int xloc2, int yloc2,
				   int line_style, int great_circle, int clevel)
{
  double mlat, mlon;    /* midpoint lat and long */
  int mxloc, myloc;     /* transformed */
  int mpx, mpy;         /* from given x,y */
  int inregion;

/* ra difference should be less than 180 degrees: take shortest path */
  if ((xfs_proj_mode == STEREOGR) || (xfs_proj_mode == GNOMONIC)
      || (xfs_proj_mode == ORTHOGR))
    if ((lon1 - lon2) > 180.0) lon1 -= 360.0;
    else if ((lon2 - lon1) > 180.0) lon2 -= 360.0;

  if (great_circle) {
    gcmidpoint(lat1, lon1, lat2, lon2, &mlat, &mlon);
  } else {
    mlat = (lat1 + lat2)/2;
    mlon = (lon1 + lon2)/2;
  };



  xform(mlat, mlon, &mxloc, &myloc, &inregion);

  mpx = (xloc1 + xloc2) / 2;
  mpy = (yloc1 + yloc2) / 2;

  if (((abs(mpx - mxloc) + abs(mpy - myloc)) > 0) && (clevel < 100)) {
    /* center is not where it should be */
    drawcurveline(lat1, lon1, mlat, mlon, xloc1, yloc1,
                  mxloc, myloc, line_style, great_circle, ++clevel);
    drawcurveline(mlat, mlon, lat2, lon2, mxloc, myloc,
                  xloc2, yloc2, line_style, great_circle, ++clevel);
  } else {
	MoveTo(xloc1, yloc1);
	LineTo(xloc2, yloc2);
  }
}

double Starmap::evalPoly(double a0, double a1, double a2, double a3, double t)
{
    return(a0 + a1 * t + a2 * t * t + a3 * t * t * t);
}

double Starmap::aint(double z)
{
	int trunk;

	trunk = (int) z;
	z = (double) trunk;
	return z;
}

double Starmap::range(double val)
{
	while (val < 0.0) {
		val = val + 360.0;
	}
	if (val > 360.0) {
		val = val - (aint(val / 360.0) * 360.0);
	}
	return val;
}

double Starmap::kepler(double e, double M)
{
	double corr, e0, E0, E1;
	e0 = e / radn;
	corr = 1;
	M = range(M);
	E0 = M;
	while (corr > 0.000001) {
		corr = (M + e0 * sin(E0 * radn) - E0)/(1 - e * cos(E0 * radn));
		E1 = E0 + corr;
		if (corr < 0) {
			corr = -1.0 * corr ;
		}
		E0 = E1;
	}

	return E1;
}

double Starmap::truean(double e, double E)
{
	double nu;

	nu = 2.0 * atan(sqrt((1 + e)/(1 - e)) * tan(E * radn / 2));
	nu = nu / radn;
	if (nu < 0.0) {
		nu = nu + 360.0;
	}
	if (fabs(nu - E) > 90.0) {
		if (nu > 180.0) {
			nu = nu - 180.0;
		} else {
			nu = nu + 180.0;
		}
	}
	return nu;
}

double Starmap::longi(double w2, double i, double u)
{
	double x, y, l;

	y = cos(i * radn) * sin(u * radn);
	x = cos(u * radn);
	l = atan2(y, x);
	l = l / radn;
	if (l < 0.0) {
		l = l + 360.0;
	}
	return l+ w2;
}

double Starmap::lati(double u, double i)
{
	double b;

	b = asin(sin(u * radn) * sin(i * radn)) / radn;
	return b;
}

void Starmap::speak(int which, double ra, double dec, double dis, int mag,
	struct planet *planet_info)
{
	if ( ra < 0) {
		ra = ra + 360.0;
	}
	planet_info[which].ra = ra;
	planet_info[which].dec = dec;
	planet_info[which].dist = dis;
	planet_info[which].mag = mag / 100.0;
}

void Starmap::trans(int which, double r, double b, double ll, double Stheta, double Sr,
		  double epli, int mag, struct planet *planet_info)
{
	double N, D, lambda, delta, beta, RA, DECLI;

	N = r * cos(b * radn) * sin((ll - Stheta) * radn);
	D = r * cos(b * radn) * cos((ll - Stheta) * radn) + Sr;
	lambda = atan2(N, D) / radn;
	if (lambda < 0.0) {
		lambda = lambda + 360.0;
	}
	lambda = range(lambda + Stheta);
	delta = sqrt(N * N + D * D + (r * sin(b * radn)) * (r * sin(b * radn)));
	beta = asin(r * sin(b * radn) / delta) / radn;
	N = sin(lambda * radn) * cos(epli * radn)
	    - tan(beta * radn) * sin(epli * radn);
	D = cos(lambda * radn);
	RA = atan2(N, D) / radn;
	DECLI = asin(sin(beta * radn) * cos(epli * radn)
	     + cos(beta * radn) * sin(epli * radn) * sin(lambda * radn)) / radn;

	planet_info[which].hlong = range(ll);
	planet_info[which].hlat = b;
	planet_info[which].hrv = r;
	speak(which, RA, DECLI, delta, mag, planet_info);
}

void Starmap::do_mercury(double T0, struct planet *planet_info)
{

	double a0, a1, a2, a3;

	pie = PI;
	radn = pie / 180.0;
	
	/* Calculate orbital elements for Mercury */
	a0=178.179078;
	a1=149474.07078;
	a2=0.0003011;
	a3=0.0;
	plan_l = evalPoly(a0,a1,a2,a3,T0);
	plan_l = range(plan_l);
	plan_a = 0.3870986;
	a0 = 0.20561421;
	a1 = 0.00002046;
	a2 = -0.000000030;
	plan_e = evalPoly(a0,a1,a2,a3,T0);
	a0 = 7.002881;
	a1 = 0.0018608;
	a2 = -0.0000183;
	plan_i = evalPoly(a0,a1,a2,a3,T0);
	a0 = 28.753753;
	a1 = 0.3702806;
	a2 = 0.0001208;
	plan_w1 = evalPoly(a0,a1,a2,a3,T0);
	plan_w1 = range(plan_w1);
	a0 = 47.145944;
	a1 = 1.1852083;
	a2= 0.0001739;
	plan_w2 = evalPoly(a0,a1,a2,a3,T0);
	plan_w2 = range(plan_w2);
	plan_M1 = 102.27938 + 149472.51529*T0 + 0.000007*T0*T0;
	plan_M1 = range(plan_M1);

	/* solving Kepler find the eccentric anomaly ECC    */

	plan_ECC =kepler(plan_e,plan_M1);
	plan_nu = truean(plan_e,plan_ECC);
	plan_r = plan_a*(1.0 - plan_e * cos(plan_ECC*pie/180.0));
	plan_u = plan_l + plan_nu - plan_M1 - plan_w2;
	plan_ll = longi(plan_w2,plan_i,plan_u);
	plan_b = lati(plan_u,plan_i);

	/* Now make corrections due to perturbations */
	plan_M2 = 212.60322 + 58517.80387*T0 +0.001286*T0*T0;
	plan_M2 = range(plan_M2);
	plan_M4 = 319.51913 + 19139.85475*T0 + 0.000181*T0*T0;
	plan_M4 = range(plan_M4);
	plan_M5 = 225.32833 + 3034.69202*T0 - 0.000722*T0*T0;
	plan_M5 = range(plan_M5);
	plan_M6 = 175.46622 +1221.55147*T0 - 0.000502*T0*T0;
	plan_M6 = range(plan_M6);
	plan_lonpert = 0.00204*cos((5*plan_M2-2*plan_M1+12.220)*radn)
		 +0.00103*cos((2*plan_M2-plan_M1-160.692)*radn)
		 +0.00091*cos((2*plan_M5-plan_M1-37.003)*radn)
		 +0.00078*cos((5*plan_M2-3*plan_M1+10.137)*radn);

	plan_radpert = 0.000007525*cos((2*plan_M5-plan_M1+53.013)*radn)
		 +0.000006802*cos((5*plan_M2-3*plan_M1-259.918)*radn)
		 +0.000005457*cos((2*plan_M2-2*plan_M1-71.188)*radn)
		 +0.000003569*cos((5*plan_M2-plan_M1-77.75)*radn);

	plan_r = plan_r + plan_radpert;
	plan_ll = plan_ll + plan_lonpert;
	/* Transformation of coordinates on Mercury and output */
	trans(1, plan_r,plan_b,plan_ll,plan_Stheta,plan_Sr,plan_epli, MAGMER, planet_info);	
}

void Starmap::do_venus(double T0, struct planet *planet_info)
{

	double a0, a1, a2, a3;

	pie = PI;
	radn = pie / 180.0;	
	
	/* Now start on Venus */
	a0 = 342.767053;
	a1 = 58519.21191;
	a2 = 0.0003097;
	a3 = 0.0;
	plan_l = evalPoly(a0,a1,a2,a3,T0);
	plan_l = range(plan_l);
	plan_a = 0.7233316;
	a0 = 0.00682069;
	a1 = -0.00004774;
	a2 = 0.000000091;
	plan_e = evalPoly(a0,a1,a2,a3,T0);
	a0 = 3.393631;
	a1 = 0.0010058;
	a2 = -0.0000010;
	plan_i = evalPoly(a0,a1,a2,a3,T0);
	a0 = 54.384186;
	a1 = 0.5081861;
	a2 = -0.0013864;
	plan_w1 = evalPoly(a0,a1,a2,a3,T0);
	plan_w1 = range(plan_w1);
	a0 = 75.779647;
	a1 = 0.8998500;
	a2 = 0.0004100;
	plan_w2 = evalPoly(a0,a1,a2,a3,T0);
	plan_w2 = range(plan_w2);
	/* Venus has a long period pert that needs to be added before Kelper */
	plan_lonpert = 0.00077 *sin((237.24 + 150.27*T0)*radn);
	plan_l = plan_l + plan_lonpert;
	plan_M0 = plan_M2 + plan_lonpert;
	plan_ECC = kepler(plan_e,plan_M0);
	plan_nu = truean(plan_e,plan_ECC);
	plan_r = plan_a*(1.0 - plan_e * cos(plan_ECC*radn));
	plan_u = plan_l + plan_nu - plan_M0- plan_w2;
	plan_ll = longi(plan_w2,plan_i,plan_u);
	plan_b = lati(plan_u,plan_i);

	/* now Venus perturbations */

	plan_lonpert = 0.00313*cos((2*plan_M-2*plan_M2 -148.225)*radn)
		 +0.00198*cos((3*plan_M-3*plan_M2 +2.565)*radn)
		 +0.00136*cos((plan_M-plan_M2-119.107)*radn)
		 +0.00096*cos((3*plan_M-2*plan_M2-135.912)*radn)
		 +0.00082*cos((plan_M5-plan_M2-208.087)*radn);
	plan_ll = plan_ll + plan_lonpert;

	plan_radpert = 0.000022501 * cos((2*plan_M-2*plan_M2-58.208)*radn)
		 +0.000019045 * cos((3*plan_M-3*plan_M2+92.577)*radn)
		 +0.000006887 * cos((plan_M5-plan_M2-118.090)*radn)
		 +0.000005172 * cos((plan_M-plan_M2-29.110)*radn)
		 +0.000003620 * cos((5*plan_M-4*plan_M2-104.208)*radn)
		 +0.000003283 * cos((4*plan_M-4*plan_M2+63.513)*radn)
		 +0.000003074 * cos((2*plan_M5-2*plan_M2-55.167)*radn);
	plan_r = plan_r + plan_radpert;
	trans(2, plan_r,plan_b,plan_ll,plan_Stheta,plan_Sr,plan_epli, MAGVEN, planet_info);	
}

void Starmap::do_mars(double T0, struct planet *planet_info)
{
	double a0, a1, a2, a3;

	pie = PI;
	radn = pie / 180.0;
		
	/* Now	start the planet Mars */
	a0 = 293.737334;
	a1 = 19141.69551;
	a2 = 0.0003107;
	a3 = 0.0;
	plan_l = evalPoly(a0,a1,a2,a3,T0);
	plan_l = range(plan_l);
	plan_a = 1.5236883;
	a0 = 0.09331290;
	a1 = 0.000092064;
	a2 = -0.000000077;
	plan_e = evalPoly(a0,a1,a2,a3,T0);
	a0 = 1.850333;
	a1 = -0.0006750;
	a2 = 0.0000126;
	plan_i = evalPoly(a0,a1,a2,a3,T0);
	a0 = 285.431761;
	a1 = 1.0697667;
	a2 = 0.0001313;
	a3 = 0.00000414;
	plan_w1 = evalPoly(a0,a1,a2,a3,T0);
	plan_w1 = range(plan_w1);
	a0 = 48.786442;
	a1 = 0.7709917;
	a2 = -0.0000014;
	a3 = -0.00000533;
	plan_w2 = evalPoly(a0,a1,a2,a3,T0);
	plan_w2 = range(plan_w2);

	/* Mars has a long period perturbation */
	plan_lonpert = -0.01133*sin((3*plan_M5-8*plan_M4 +4*plan_M)*radn)
		  -0.00933*cos((3*plan_M5-8*plan_M4 +4*plan_M)*radn);
	plan_l = plan_l + plan_lonpert;
	plan_M0 = plan_M4 + plan_lonpert;
	plan_ECC = kepler(plan_e,plan_M0);
	plan_nu = truean(plan_e,plan_ECC);
	plan_r = plan_a*(1.0 - plan_e * cos(plan_ECC*radn));
	plan_u = plan_l + plan_nu - plan_M0- plan_w2;
	plan_ll = longi(plan_w2,plan_i,plan_u);
	plan_b = lati(plan_u,plan_i);

	/* Now Mars Perturbations */
	plan_lonpert = 0.00705*cos((plan_M5-plan_M4-48.958)*radn)
		 +0.00607*cos((2*plan_M5-plan_M4-188.350)*radn)
		 +0.00445*cos((2*plan_M5-2*plan_M4-191.897)*radn)
		 +0.00388*cos((plan_M-2*plan_M4+20.495)*radn)
		 +0.00238*cos((plan_M-plan_M4+35.097)*radn)
		 +0.00204*cos((2*plan_M-3*plan_M4+158.638)*radn)
		 +0.00177*cos((3*plan_M4-plan_M2-57.602)*radn)
		 +0.00136*cos((2*plan_M-4*plan_M4+154.093)*radn)
		 +0.00104*cos((plan_M5+17.618)*radn);
	plan_ll = plan_ll + plan_lonpert;

	plan_radpert = 0.000053227*cos((plan_M5-plan_M4+41.1306)*radn)
		 +0.000050989*cos((2*plan_M5-2*plan_M4-101.9847)*radn)
		 +0.000038278*cos((2*plan_M5-plan_M4-98.3292)*radn)
		 +0.000015996*cos((plan_M-plan_M4-55.555)*radn)
		 +0.000014764*cos((2*plan_M-3*plan_M4+68.622)*radn)
		 +0.000008966*cos((plan_M5-2*plan_M4+43.615)*radn)
		 +0.000007914*cos((3*plan_M5-2*plan_M4-139.737)*radn)
		 +0.000007004*cos((2*plan_M5-3*plan_M4-102.888)*radn)
		 +0.000006620*cos((plan_M-2*plan_M4+113.202)*radn)
		 +0.000004930*cos((3*plan_M5-3*plan_M4-76.243)*radn)
		 +0.000004693*cos((3*plan_M-5*plan_M4+190.603)*radn)
		 +0.000004571*cos((2*plan_M-4*plan_M4+244.702)*radn)
		 +0.000004409*cos((3*plan_M5-plan_M4-115.828)*radn);
	plan_r = plan_r + plan_radpert;
	trans(4, plan_r,plan_b,plan_ll,plan_Stheta,plan_Sr,plan_epli, MAGMAR, planet_info);	
}

void Starmap::do_jupiter(double T0, struct planet *planet_info)
{
	double a0, a1, a2, a3;

	pie = PI;
	radn = pie / 180.0;

	/* Now start Jupiter */
	a0 = 238.049257;
	a1 = 3036.301986;
	a2 = 0.0003347;
	a3 = -0.00000165;
	plan_l = evalPoly(a0,a1,a2,a3,T0);
	plan_l = range(plan_l);
	plan_a = 5.202561;
	a0 = 0.04833475;
	a1 = 0.000164180;
	a2 = -0.0000004676;
	a3 = -0.0000000017;
	plan_e = evalPoly(a0,a1,a2,a3,T0);
	a0 = 1.308736;
	a1 = -0.0056961;
	a2 = 0.0000039;
	a3 = 0.0;
	plan_i = evalPoly(a0,a1,a2,a3,T0);
	a0 = 273.277558;
	a1 = 0.5994317;
	a2 = 0.00070405;
	a3 = 0.00000508;
	plan_w1 = evalPoly(a0,a1,a2,a3,T0);
	plan_w1 = range(plan_w1);
	a0 = 99.443414;
	a1 = 1.0105300;
	a2 = 0.00035222;
	a3 = -0.00000851;
	plan_w2 = evalPoly(a0,a1,a2,a3,T0);
	plan_w2 = range(plan_w2);

	/* Now start perturbation calculations */

	plan_nu2 = T0/5.0 +0.1;
	plan_P = 237.47555 +3034.9061*T0;
	plan_Q = 265.91650 + 1222.1139*T0;
	plan_S = 243.51721 + 428.4677*T0;
	plan_V = 5.0*plan_Q -2.0*plan_P;
	plan_W = 2.0*plan_P - 6.0*plan_Q +3.0*plan_S;
	plan_ze = plan_Q - plan_P;
	plan_psi = plan_S - plan_Q;

	plan_P = range(plan_P)*radn;
	plan_Q = range(plan_Q)*radn;
	plan_S = range(plan_S)*radn;
	plan_V = range(plan_V)*radn;
	plan_W = range(plan_W)*radn;
	plan_ze = range(plan_ze)*radn;
	plan_psi = range(plan_psi)*radn;

	plan_l1pert = (0.331364 - 0.010281*plan_nu2 - 0.004692*plan_nu2*plan_nu2)*sin(plan_V)
		+(0.003228 - 0.064436*plan_nu2 + 0.002075*plan_nu2*plan_nu2)*cos(plan_V)
		-(0.003083 + 0.000275*plan_nu2 - 0.000489*plan_nu2*plan_nu2)*sin(2*plan_V)
		+0.002472*sin(plan_W)
		+0.013619*sin(plan_ze)
		+0.018472*sin(2*plan_ze)
		+0.006717*sin(3*plan_ze)
		+0.002775*sin(4*plan_ze)
		+(0.007275 - 0.001253*plan_nu2)*sin(plan_ze)*sin(plan_Q)
		+0.006417*sin(2*plan_ze)*sin(plan_Q)
		+0.002439*sin(3*plan_ze)*sin(plan_Q);

	plan_l1pert = plan_l1pert -(0.033839 + 0.001125*plan_nu2)*cos(plan_ze)*sin(plan_Q)
		-0.003767*cos(2*plan_ze)*sin(plan_Q)
		-(0.035681 + 0.001208*plan_nu2)*sin(plan_ze)*cos(plan_Q)
		-0.004261*sin(2*plan_ze)*cos(plan_Q)
		+0.002178*cos(plan_Q)
		+(-0.006333 + 0.001161*plan_nu2)*cos(plan_ze)*cos(plan_Q)
		-0.006675*cos(2*plan_ze)*cos(plan_Q)
		-0.002664*cos(3*plan_ze)*cos(plan_Q)
		-0.002572*sin(plan_ze)*sin(2*plan_Q)
		-0.003567*sin(2*plan_ze)*sin(2*plan_Q)
		+0.002094*cos(plan_ze)*cos(2*plan_Q)
		+0.003342*cos(2*plan_ze)*cos(2*plan_Q);

	plan_epert =  (.0003606 + .0000130*plan_nu2 - .0000043*plan_nu2*plan_nu2)*sin(plan_V)
		+(.0001289 - .0000580*plan_nu2)*cos(plan_V)
		-.0006764*sin(plan_ze)*sin(plan_Q)
		-.0001110*sin(2*plan_ze)*sin(plan_Q)
		-.0000224*sin(3*plan_ze)*sin(plan_Q)
		-.0000204*sin(plan_Q)
		+(.0001284 + .0000116*plan_nu2)*cos(plan_ze)*sin(plan_Q)
		+.0000188*cos(2*plan_ze)*sin(plan_Q)
		+(.0001460 + .0000130*plan_nu2)*sin(plan_ze)*cos(plan_Q)
		+.0000224*sin(2*plan_ze)*cos(plan_Q)
		-.0000817*cos(plan_Q);

	plan_epert = plan_epert	+.0006074*cos(plan_ze)*cos(plan_Q)
		+.0000992*cos(2*plan_ze)*cos(plan_Q)
		+.0000508*cos(3*plan_ze)*cos(plan_Q)
		+.0000230*cos(4*plan_ze)*cos(plan_Q)
		+.0000108*cos(5*plan_ze)*cos(plan_Q)
		-(.0000956 + .0000073*plan_nu2)*sin(plan_ze)*sin(2*plan_Q)
		+.0000448*sin(2*plan_ze)*sin(2*plan_Q)
		+.0000137*sin(3*plan_ze)*sin(2*plan_Q)
		+(-.0000997 + .0000108*plan_nu2)*cos(plan_ze)*sin(2*plan_Q)
		+.0000480*cos(2*plan_ze)*sin(2*plan_Q);

	plan_epert = plan_epert	+.0000148*cos(3*plan_ze)*sin(2*plan_Q)
		+(-.0000956 +.0000099*plan_nu2)*sin(plan_ze)*cos(2*plan_Q)
		+.0000490*sin(2*plan_ze)*cos(2*plan_Q)
		+.0000158*sin(3*plan_ze)*cos(2*plan_Q)
		+.0000179*cos(2*plan_Q)
		+(.0001024 + .0000075*plan_nu2)*cos(plan_ze)*cos(2*plan_Q)
		-.0000437*cos(2*plan_ze)*cos(2*plan_Q)
		-.0000132*cos(3*plan_ze)*cos(2*plan_Q);

	plan_w1pert = (0.007192 - 0.003147*plan_nu2)*sin(plan_V)
		+(-0.020428 - 0.000675*plan_nu2 + 0.000197*plan_nu2*plan_nu2)*cos(plan_V)
		+(0.007269 + 0.000672*plan_nu2)*sin(plan_ze)*sin(plan_Q)
		-0.004344*sin(plan_Q)
		+0.034036*cos(plan_ze)*sin(plan_Q)
		+0.005614*cos(2*plan_ze)*sin(plan_Q)
		+0.002964*cos(3*plan_ze)*sin(plan_Q)
		+0.037761*sin(plan_ze)*cos(plan_Q);

	plan_w1pert = plan_w1pert
		+0.006158*sin(2*plan_ze)*cos(plan_Q)
		-0.006603*cos(plan_ze)*cos(plan_Q)
		-0.005356*sin(plan_ze)*sin(2*plan_Q)
		+0.002722*sin(2*plan_ze)*sin(2*plan_Q)
		+0.004483*cos(plan_ze)*sin(2*plan_Q)
		-0.002642*cos(2*plan_ze)*sin(2*plan_Q)
		+0.004403*sin(plan_ze)*cos(2*plan_Q)
		-0.002536*sin(2*plan_ze)*cos(2*plan_Q)
		+0.005547*cos(plan_ze)*cos(2*plan_Q)
		-0.002689*cos(2*plan_ze)*cos(2*plan_Q);

	plan_lonpert = plan_l1pert -(plan_w1pert/plan_e);
	plan_l = plan_l + plan_l1pert;
	plan_M0 = plan_M5 + plan_lonpert;
	plan_e = plan_e + plan_epert;
	plan_w1 = plan_w1 + plan_w1pert;

	plan_apert =  -.000263*cos(plan_V)
		+.000205*cos(plan_ze)
		+.000693*cos(2*plan_ze)
		+.000312*cos(3*plan_ze)
		+.000147*cos(4*plan_ze)
		+.000299*sin(plan_ze)*sin(plan_Q)
		+.000181*cos(2*plan_ze)*sin(plan_Q)
		+.000204*sin(2*plan_ze)*cos(plan_Q)
		+.000111*sin(3*plan_ze)*cos(plan_Q)
		-.000337*cos(plan_ze)*cos(plan_Q)
		-.000111*cos(2*plan_ze)*cos(plan_Q);

	plan_a = plan_a + plan_apert;
	plan_ECC = kepler(plan_e,plan_M0);
	plan_nu = truean(plan_e,plan_ECC);
	plan_r = plan_a*(1.0 - plan_e * cos(plan_ECC*radn));
	plan_u = plan_l + plan_nu - plan_M0 - plan_w2;
	plan_ll = longi(plan_w2,plan_i,plan_u);
	plan_b = lati(plan_u,plan_i);

	trans(5, plan_r,plan_b,plan_ll,plan_Stheta,plan_Sr,plan_epli, MAGJUP, planet_info);	
}

void Starmap::do_saturn(double T0, struct planet *planet_info)
{
	double a0, a1, a2, a3;

	pie = PI;
	radn = pie / 180.0;	
	
	/* Now start Saturn */

	a0 = 266.564377;
	a1 = 1223.509884;
	a2 = 0.0003245;
	a3 = -0.0000058;
	plan_l = evalPoly(a0,a1,a2,a3,T0);
	plan_l = range(plan_l);
	plan_a = 9.554747;
	a0 = 0.05589232;
	a1 = -0.00034550;
	a2 = -0.000000728;
	a3 = 0.00000000074;
	plan_e = evalPoly(a0,a1,a2,a3,T0);
	a0 = 2.492519;
	a1 = -0.0039189;
	a2 = -0.00001549;
	a3 = 0.00000004;
	plan_i = evalPoly(a0,a1,a2,a3,T0);
	a0 = 338.307800;
	a1 = 1.0852207;
	a2 = 0.00097854;
	a3 = 0.00000992;
	plan_w1 = evalPoly(a0,a1,a2,a3,T0);
	plan_w1 = range(plan_w1);
	a0 = 112.790414;
	a1 = 0.8731951;
	a2 = -0.00015218;
	a3 = -0.00000531;
	plan_w2 = evalPoly(a0,a1,a2,a3,T0);
	plan_w2 = range(plan_w2);

	/* Now Saturn's perturbations */

	plan_l1pert = (-0.814181 + 0.018150*plan_nu2 + 0.016714*plan_nu2*plan_nu2)*sin(plan_V)
		 +(-0.010497 + 0.160906*plan_nu2 - 0.004100*plan_nu2*plan_nu2)*cos(plan_V)
		 +0.007581*sin(2*plan_V)
		 -0.007986*sin(plan_W)
		 -0.148811*sin(plan_ze)
		 -0.040786*sin(2*plan_ze)
		 -0.015208*sin(3*plan_ze)
		 -0.006339*sin(4*plan_ze)
		 -0.006244*sin(plan_Q);
	plan_l1pert = plan_l1pert
		+(0.008931 + 0.002728*plan_nu2)*sin(plan_ze)*sin(plan_Q)
		-0.016500*sin(2*plan_ze)*sin(plan_Q)
		-0.005775*sin(3*plan_ze)*sin(plan_Q)
		+(0.081344 + 0.003206*plan_nu2)*cos(plan_ze)*sin(plan_Q)
		+0.015019*cos(2*plan_ze)*sin(plan_Q)
		+(0.085581 + 0.002494*plan_nu2)*sin(plan_ze)*cos(plan_Q)
		+(0.025328 - 0.003117*plan_nu2)*cos(plan_ze)*cos(plan_Q);
	plan_l1pert = plan_l1pert
		+0.014394*cos(2*plan_ze)*cos(plan_Q)
		+0.006319*cos(3*plan_ze)*cos(plan_Q)
		+0.006369*sin(plan_ze)*sin(2*plan_Q)
		+0.009156*sin(2*plan_ze)*sin(2*plan_Q)
		+0.007525*sin(3*plan_psi)*sin(2*plan_Q)
		-0.005236*cos(plan_ze)*cos(2*plan_Q)
		-0.007736*cos(2*plan_ze)*cos(2*plan_Q)
		-0.007528*cos(3*plan_psi)*cos(2*plan_Q);

	plan_epert = (-.0007927 + .0002548*plan_nu2 +.0000091*plan_nu2*plan_nu2)*sin(plan_V)
		+(.0013381 + .0001226*plan_nu2 -.0000253*plan_nu2*plan_nu2)*cos(plan_V)
		+(.0000248 - .0000121*plan_nu2)*sin(2*plan_V)
		-(.0000305 + .0000091*plan_nu2)*cos(2*plan_V)
		+.0000412*sin(2*plan_ze)
		+.0012415*sin(plan_Q)
		+(.0000390 -.0000617*plan_nu2)*sin(plan_ze)*sin(plan_Q)
		+(.0000165 - .0000204*plan_nu2)*sin(2*plan_ze)*sin(plan_Q)
		+.0026599*cos(plan_ze)*sin(plan_Q)
		-.0004687*cos(2*plan_ze)*sin(plan_Q);
	plan_epert = plan_epert
		-.0001870*cos(3*plan_ze)*sin(plan_Q)
		-.0000821*cos(4*plan_ze)*sin(plan_Q)
		-.0000377*cos(5*plan_ze)*sin(plan_Q)
		+.0000497*cos(2*plan_psi)*sin(plan_Q)
		+(.0000163 - .0000611*plan_nu2)*cos(plan_Q)
		-.0012696*sin(plan_ze)*cos(plan_Q)
		-.0004200*sin(2*plan_ze)*cos(plan_Q)
		-.0001503*sin(3*plan_ze)*cos(plan_Q)
		-.0000619*sin(4*plan_ze)*cos(plan_Q)
		-.0000268*sin(5*plan_ze)*cos(plan_Q);
	plan_epert = plan_epert
		-(.0000282 + .0001306*plan_nu2)*cos(plan_ze)*cos(plan_Q)
		+(-.0000086 + .0000230*plan_nu2)*cos(2*plan_ze)*cos(plan_Q)
		+.0000461*sin(2*plan_psi)*cos(plan_Q)
		-.0000350*sin(2*plan_Q)
		+(.0002211 - .0000286*plan_nu2)*sin(plan_ze)*sin(2*plan_Q)
		-.0002208*sin(2*plan_ze)*sin(2*plan_Q)
		-.0000568*sin(3*plan_ze)*sin(2*plan_Q)
		-.0000346*sin(4*plan_ze)*sin(2*plan_Q)
		-(.0002780 + .0000222*plan_nu2)*cos(plan_ze)*sin(2*plan_Q)
		+(.0002022 + .0000263*plan_nu2)*cos(2*plan_ze)*sin(2*plan_Q);
	plan_epert = plan_epert
		+.0000248*cos(3*plan_ze)*sin(2*plan_Q)
		+.0000242*sin(3*plan_psi)*sin(2*plan_Q)
		+.0000467*cos(3*plan_psi)*sin(2*plan_Q)
		-.0000490*cos(2*plan_Q)
		-(.0002842 + .0000279*plan_nu2)*sin(plan_ze)*cos(2*plan_Q)
		+(.0000128 + .0000226*plan_nu2)*sin(2*plan_ze)*cos(2*plan_Q)
		+.0000224*sin(3*plan_ze)*cos(2*plan_Q)
		+(-.0001594 + .0000282*plan_nu2)*cos(plan_ze)*cos(2*plan_Q)
		+(.0002162 - .0000207*plan_nu2)*cos(2*plan_ze)*cos(2*plan_Q)
		+.0000561*cos(3*plan_ze)*cos(2*plan_Q);
	plan_epert = plan_epert
		+.0000343*cos(4*plan_ze)*cos(2*plan_Q)
		+.0000469*sin(3*plan_psi)*cos(2*plan_Q)
		-.0000242*cos(3*plan_psi)*cos(2*plan_Q)
		-.0000205*sin(plan_ze)*sin(3*plan_Q)
		+.0000262*sin(3*plan_ze)*sin(3*plan_Q)
		+.0000208*cos(plan_ze)*cos(3*plan_Q)
		-.0000271*cos(3*plan_ze)*cos(3*plan_Q)
		-.0000382*cos(3*plan_ze)*sin(4*plan_Q)
		-.0000376*sin(3*plan_ze)*cos(4*plan_Q);

	plan_w1pert = (0.077108 + 0.007186*plan_nu2 - 0.001533*plan_nu2*plan_nu2)*sin(plan_V)
		+(0.045803 - 0.014766*plan_nu2 - 0.000536*plan_nu2*plan_nu2)*cos(plan_V)
		-0.007075*sin(plan_ze)
		-0.075825*sin(plan_ze)*sin(plan_Q)
		-0.024839*sin(2*plan_ze)*sin(plan_Q)
		-0.008631*sin(3*plan_ze)*sin(plan_Q)
		-0.072586*cos(plan_Q)
		-0.150383*cos(plan_ze)*cos(plan_Q)
		+0.026897*cos(2*plan_ze)*cos(plan_Q)
		+0.010053*cos(3*plan_ze)*cos(plan_Q);
	plan_w1pert = plan_w1pert
		-(0.013597 +0.001719*plan_nu2)*sin(plan_ze)*sin(2*plan_Q)
		+(-0.007742 + 0.001517*plan_nu2)*cos(plan_ze)*sin(2*plan_Q)
		+(0.013586 - 0.001375*plan_nu2)*cos(2*plan_ze)*sin(2*plan_Q)
		+(-0.013667 + 0.001239*plan_nu2)*sin(plan_ze)*cos(2*plan_Q)
		+0.011981*sin(2*plan_ze)*cos(2*plan_Q)
		+(0.014861 + 0.001136*plan_nu2)*cos(plan_ze)*cos(2*plan_Q)
		-(0.013064 + 0.001628*plan_nu2)*cos(2*plan_ze)*cos(2*plan_Q);

	plan_lonpert = plan_l1pert -(plan_w1pert/plan_e);
	plan_l = plan_l + plan_l1pert;
	plan_M0 = plan_M6 + plan_lonpert;
	plan_e = plan_e + plan_epert;
	plan_w1 = plan_w1 + plan_w1pert;

	plan_apert = .000572*sin(plan_V) -.001590*sin(2*plan_ze)*cos(plan_Q)
		+.002933*cos(plan_V) -.000647*sin(3*plan_ze)*cos(plan_Q)
		+.033629*cos(plan_ze) -.000344*sin(4*plan_ze)*cos(plan_Q)
		-.003081*cos(2*plan_ze) +.002885*cos(plan_ze)*cos(plan_Q)
		-.001423*cos(3*plan_ze) +(.002172 + .000102*plan_nu2)*cos(2*plan_ze)*cos(plan_Q)
		-.000671*cos(4*plan_ze) +.000296*cos(3*plan_ze)*cos(plan_Q)
		-.000320*cos(5*plan_ze) -.000267*sin(2*plan_ze)*sin(2*plan_Q);
	plan_apert = plan_apert
		+.001098*sin(plan_Q) -.000778*cos(plan_ze)*sin(2*plan_Q)
		-.002812*sin(plan_ze)*sin(plan_Q) +.000495*cos(2*plan_ze)*sin(2*plan_Q)
		+.000688*sin(2*plan_ze)*sin(plan_Q) +.000250*cos(3*plan_ze)*sin(2*plan_Q);
	plan_apert = plan_apert
		-.000393*sin(3*plan_ze)*sin(plan_Q)
		-.000228*sin(4*plan_ze)*sin(plan_Q)
		+.002138*cos(plan_ze)*sin(plan_Q)
		-.000999*cos(2*plan_ze)*sin(plan_Q)
		-.000642*cos(3*plan_ze)*sin(plan_Q)
		-.000325*cos(4*plan_ze)*sin(plan_Q)
		-.000890*cos(plan_Q)
		+.002206*sin(plan_ze)*cos(plan_Q);
	plan_apert = plan_apert
		-.000856*sin(plan_ze)*cos(2*plan_Q)
		+.000441*sin(2*plan_ze)*cos(2*plan_Q)
		+.000296*cos(2*plan_ze)*cos(2*plan_Q)
		+.000211*cos(3*plan_ze)*cos(2*plan_Q)
		-.000427*sin(plan_ze)*sin(3*plan_Q)
		+.000398*sin(3*plan_ze)*sin(3*plan_Q)
		+.000344*cos(plan_ze)*cos(3*plan_Q)
		-.000427*cos(3*plan_ze)*cos(3*plan_Q);

	plan_a = plan_a + plan_apert;
	plan_ECC = kepler(plan_e,plan_M0);
	plan_nu = truean(plan_e,plan_ECC);
	plan_r = plan_a*(1.0 - plan_e * cos(plan_ECC*radn));
	plan_u = plan_l + plan_nu - plan_M0 - plan_w2;
	plan_ll = longi(plan_w2,plan_i,plan_u);
	plan_b = lati(plan_u,plan_i);

	plan_b = plan_b
		+0.000747*cos(plan_ze)*sin(plan_Q)
		+0.001069*cos(plan_ze)*cos(plan_Q)
		+0.002108*sin(2*plan_ze)*sin(2*plan_Q)
		+0.001261*cos(2*plan_ze)*sin(2*plan_Q)
		+0.001236*sin(2*plan_ze)*cos(2*plan_Q)
		-0.002075*cos(2*plan_ze)*cos(2*plan_Q);

	trans(6, plan_r,plan_b,plan_ll,plan_Stheta,plan_Sr,plan_epli, MAGSAT, planet_info);	
}

void Starmap::do_uranus(double T0, struct planet *planet_info)
{
	double a0, a1, a2, a3;

	pie = PI;
	radn = pie / 180.0;
		
	/* Now Start on Uranus */
	a0 = 244.197470;
	a1 = 429.863546;
	a2 = 0.0003160;
	a3 = -0.00000060;
	plan_l = evalPoly(a0,a1,a2,a3,T0);
	plan_l = range(plan_l);
	plan_a = 19.21814;
	a0 = 0.0463444;
	a1 = -0.00002658;
	a2 = 0.000000077;
	a3 = 0.0;
	plan_e = evalPoly(a0,a1,a2,a3,T0);
	a0 = 0.772464;
	a1 = 0.0006253;
	a2 = 0.0000395;
	plan_i = evalPoly(a0,a1,a2,a3,T0);
	a0 = 98.071581;
	a1 = 0.9857650;
	a2 = -0.0010745;
	a3 = -0.00000061;
	plan_w1 = evalPoly(a0,a1,a2,a3,T0);
	plan_w1 = range(plan_w1);
	a0 = 73.477111;
	a1 = 0.4986678;
	a2 = 0.0013117;
	a3 = 0.0;
	plan_w2 = evalPoly(a0,a1,a2,a3,T0);
	plan_w2 = range(plan_w2);
	plan_M7 = 72.64878 + 428.37911*T0 + 0.000079*T0*T0;
	plan_M7 = range(plan_M7);
	/* now perturbation corrections for Uranus */
	plan_G = 83.76922 + 218.4901*T0;
	plan_S = plan_S/radn;
	plan_P = plan_P/radn;
	plan_Q = plan_Q/radn;
	plan_H = 2.0*plan_G - plan_S;
	plan_ze = plan_S - plan_P;
	plan_eta = plan_S - plan_Q;
	plan_th = plan_G - plan_S;
	plan_S = plan_S*radn;
	plan_G = range(plan_G)*radn;
	plan_P = plan_P*radn;
	plan_Q = plan_Q*radn;
	plan_H = range(plan_H)*radn;
	plan_ze = range(plan_ze)*radn;
	plan_eta = range(plan_eta)*radn;
	plan_th = range(plan_th)*radn;

	plan_l1pert = (0.864319 - 0.001583*plan_nu2)*sin(plan_H)
		+(0.082222 - 0.006833*plan_nu2)*cos(plan_H)
		+0.036017*sin(2*plan_H)
		-0.003019*cos(2*plan_H)
		+0.008122*sin(plan_W);

	plan_epert = (-.0003349 + .0000163*plan_nu2)*sin(plan_H)
		+.0020981*cos(plan_H)
		+.0001311*cos(plan_H);

	plan_w1pert = 0.120303*sin(plan_H)
		+(0.019472 - 0.000947*plan_nu2)*cos(plan_H)
		+0.006197*sin(2*plan_H);

	plan_lonpert = plan_l1pert -(plan_w1pert/plan_e);
	plan_l = plan_l + plan_l1pert;
	plan_M0 = plan_M7 + plan_lonpert;
	plan_e = plan_e + plan_epert;
	plan_w1 = plan_w1 + plan_w1pert;

	plan_a = plan_a - 0.003825*cos(plan_H);
	plan_ECC = kepler(plan_e,plan_M0);
	plan_nu = truean(plan_e,plan_ECC);
	plan_r = plan_a*(1.0 - plan_e * cos(plan_ECC*radn));
	plan_u = plan_l + plan_nu - plan_M0 - plan_w2;
	plan_ll = longi(plan_w2,plan_i,plan_u);
	plan_b = lati(plan_u,plan_i);

	plan_ll = plan_ll
		+(0.010122 - 0.000988*plan_nu2)*sin(plan_S+plan_eta)
		+(-0.038581 + 0.002031*plan_nu2 - 0.001910*plan_nu2*plan_nu2)*cos(plan_S+plan_eta)
		+(0.034964 - 0.001038*plan_nu2 + 0.000868*plan_nu2*plan_nu2)*cos(2*plan_S+plan_eta)
		+0.005594*sin(plan_S +3*plan_th);
	plan_ll = plan_ll
		-0.014808*sin(plan_ze)
		-0.005794*sin(plan_eta)
		+0.002347*cos(plan_eta)
		+0.009872*sin(plan_th)
		+0.008803*sin(2*plan_th)
		-0.004308*sin(3*plan_th);

	plan_b = plan_b
		+(0.000458*sin(plan_eta) - 0.000642*cos(plan_eta) - 0.000517*cos(4*plan_th))
		*sin(plan_S)
		-(0.000347*sin(plan_eta) + 0.000853*cos(plan_eta) + 0.000517*sin(4*plan_eta))
		*cos(plan_S)
		+0.000403*(cos(2*plan_th)*sin(2*plan_S) + sin(2*plan_th)*cos(2*plan_S));

	plan_r = plan_r
		-.025948
		+.004985*cos(plan_ze)
		-.001230*cos(plan_S)
		+.003354*cos(plan_eta)
		+(.005795*cos(plan_S) - .001165*sin(plan_S) + .001388*cos(2*plan_S))*sin(plan_eta)
		+(.001351*cos(plan_S) + .005702*sin(plan_S) + .001388*sin(2*plan_S))*cos(plan_eta)
		+.000904*cos(2*plan_th)
		+.000894*(cos(plan_th) - cos(3*plan_th));
	trans(7, plan_r,plan_b,plan_ll,plan_Stheta,plan_Sr,plan_epli, MAGURA, planet_info);	
}

void Starmap::do_neptune(double T0, struct planet *planet_info)
{
	double a0, a1, a2, a3;

	pie = PI;
	radn = pie / 180.0;	
	
	/* Now start Neptune */
	a0 = 84.457994;
	a1 = 219.885914;
	a2 = 0.0003205;
	a3 = -0.00000060;
	plan_l = evalPoly(a0,a1,a2,a3,T0);
	plan_l = range(plan_l);
	plan_a = 30.10957;
	a0 = 0.00899704;
	a1 = 0.000006330;
	a2 = -0.000000002;
	a3 = 0.0;
	plan_e = evalPoly(a0,a1,a2,a3,T0);
	a0 = 1.779242;
	a1 = -0.0095436;
	a2 = -0.0000091;
	plan_i = evalPoly(a0,a1,a2,a3,T0);
	a0 = 276.045975;
	a1 = 0.3256394;
	a2 = 0.00014095;
	a3 = 0.000004113;
	plan_w1 = evalPoly(a0,a1,a2,a3,T0);
	plan_w1 = range(plan_w1);
	a0 = 130.681389;
	a1 = 1.0989350;
	a2 = 0.00024987;
	a3 = -0.000004718;
	plan_w2 = evalPoly(a0,a1,a2,a3,T0);
	plan_w2 = range(plan_w2);
	plan_M8 = 37.73063 + 218.46134*T0 -0.000070*T0*T0;

	/* now perturbation corrections for neptune */
	plan_G = plan_G/radn;
	plan_P = plan_P/radn;
	plan_Q = plan_Q/radn;
	plan_S = plan_S/radn;
	plan_ze = plan_G - plan_P;
	plan_eta = plan_G - plan_Q;
	plan_th = plan_G - plan_S;
	plan_G = plan_G*radn;
	plan_P = plan_P*radn;
	plan_Q = plan_Q*radn;
	plan_S = plan_S*radn;
	plan_ze = range(plan_ze)*radn;
	plan_eta = range(plan_eta)*radn;
	plan_th = range(plan_th)*radn;

	plan_l1pert = (-0.589833 + 0.001089*plan_nu2)*sin(plan_H)
		+(-0.056094 + 0.004658*plan_nu2)*cos(plan_H)
		-0.024286*sin(2*plan_H);

	plan_epert = .0004389*sin(plan_H)
		+.0004262*cos(plan_H)
		+.0001129*sin(2*plan_H)
		+.0001089*cos(2*plan_H);

	plan_w1pert = 0.024039*sin(plan_H)
		-0.025303*cos(plan_H)
		+0.006206*sin(2*plan_H)
		-0.005992*cos(2*plan_H);


	plan_lonpert = plan_l1pert -(plan_w1pert/plan_e);
	plan_l = plan_l + plan_l1pert;
	plan_M0 = plan_M8 + plan_lonpert;
	plan_e = plan_e + plan_epert;
	plan_w1 = plan_w1 + plan_w1pert;

	plan_a = plan_a - 0.000817*sin(plan_H)
		+0.008189*cos(plan_H)
		+0.000781*cos(2*plan_H);

	plan_ECC = kepler(plan_e,plan_M0);
	plan_nu = truean(plan_e,plan_ECC);
	plan_r = plan_a*(1.0 - plan_e * cos(plan_ECC*radn));
	plan_u = plan_l + plan_nu - plan_M0 - plan_w2;
	plan_ll = longi(plan_w2,plan_i,plan_u);
	plan_b = lati(plan_u,plan_i);

	plan_ll = plan_ll
		-0.009556*sin(plan_ze)
		-0.005178*sin(plan_eta)
		+0.002572*sin(2*plan_th)
		-0.002972*cos(2*plan_th)*sin(plan_G)
		-0.002833*sin(2*plan_th)*cos(plan_G);

	plan_b = plan_b
		+0.000336*cos(2*plan_th)*sin(plan_G)
		+0.000364*sin(2*plan_th)*cos(plan_G);

	plan_r = plan_r
		-.040596
		+.004992*cos(plan_ze)
		+.002744*cos(plan_eta)
		+.002044*cos(plan_th)
		+.001051*cos(2*plan_th);
	trans(8, plan_r,plan_b,plan_ll,plan_Stheta,plan_Sr,plan_epli, MAGNEP, planet_info);	
}

int Starmap::pluto(double jd, double *l, double *b, double *r)
{
	double t, j, s, p, cl, cb, cr, alpha, salpha, calpha;
	int i;

	/* This expansion for Pluto's position is valid only for years
	   between 1885 and 2099.  If the date given is outside that
	   range, zero the result cells and return FALSE. */

	if (plutoPrecise && (jd < 2409543.0 || jd >= 2488070.0)) {
		*l = *b = *r = 0.0;
		return FALSE;
	}

	t = (jd - J2000) / JulianCentury;	// Time in Julian centuries since the Epoch J2000.0

	j = dtr( 34.35 + 3034.9057 * t);
	s = dtr( 50.08 + 1222.1138 * t);
	p = dtr(238.96 +  144.9600 * t);

	cl = 238.956785 + 144.96 * t;
	cb = -3.908202;
	cr = 40.7247248;

	for (i = 0; i < ELEMENTS(pt); i++) {
		alpha = pt[i].j * j + pt[i].s * s + pt[i].p * p;
		salpha = sin(alpha);
		calpha = cos(alpha);
		cl += pt[i].longA * salpha * 0.000001 + pt[i].longB * calpha * 0.000001;
		cb += pt[i].latA * salpha * 0.000001  +  pt[i].latB * calpha * 0.000001;
		cr += pt[i].radA * salpha * 0.0000001  +  pt[i].radB * calpha * 0.0000001;
	}

	*l = cl;
	*b = cb;
	*r = cr;
	return TRUE;
}


void Starmap::planets(double jd, int which, struct planet *planet_info)
{
	double T0;
	

	pie = PI;
	radn = pie / 180.0;

/*  calculate time T0 from 1900 January 0.5 ET */
	T0 = (jd - 2415020.0) / JulianCentury;

	if (Planet(0))
	{
	/* Find the Longitude and radius vector of the sun */

	plan_M= 358.47583 + 35999.04975*T0 - 0.000150*T0*T0
	   -0.0000033*T0*T0*T0;
	plan_M = range(plan_M);
	plan_esun = 0.01675104 - 0.0000418*T0 -0.000000126*T0*T0;
	plan_Lsun = 279.69668 + 36000.76892*T0 + 0.0003025*T0*T0;
	plan_Lsun = range(plan_Lsun);
	plan_Cen= (1.919460-0.004789*T0-0.000014*T0*T0)*sin(plan_M*radn)
	     +(0.020094 - 0.000100*T0)*sin(2*plan_M*radn)
	     +0.000293*sin(3*plan_M*radn);

	plan_Stheta= plan_Lsun + plan_Cen;
	plan_Stheta = range(plan_Stheta);
	plan_Snu = plan_M + plan_Cen;
	plan_Sr = 1.0000002*(1.0 - plan_esun*plan_esun)/(1.0 + plan_esun*cos(plan_Snu*radn));

	plan_omeg= 259.18 - 1934.142*T0;
	plan_thapp = plan_Stheta - 0.00569 - 0.00479* sin(plan_omeg*radn);
	plan_epli = 23.452294 - 0.0130125*T0
	       -0.00000164*T0*T0 + 0.000000503*T0*T0*T0
		+0.00256*cos(plan_omeg*radn);
	plan_N = cos(plan_epli*radn)*sin(plan_thapp*radn);
	plan_D = cos(plan_thapp*radn);
	plan_RA = atan2(plan_N,plan_D)/radn;
	plan_DEC = asin(sin(plan_epli*radn)*sin(plan_thapp*radn))/radn;
	speak(0, plan_RA,plan_DEC,plan_Sr, MAGSOL, planet_info);
	planet_info[0].hlong = plan_Stheta;
	planet_info[0].hlat = 0;
	planet_info[0].hrv = 0;
	}
	/*assert(plan_epli > 0.0);	*/			// Make sure Sun calculated before any other planets, because
	                                    // some of the values calculated for the sun are used by the
	                                    // other planets. These values are stored as statics, so that
	                                    // if the function re-enters, they are already calculated.
	if (Planet(1))
	{
		do_mercury(T0, planet_info);
	}
	if (Planet(2))
	{
		do_venus(T0, planet_info);
	}
	// Planet 3 is actually the moon - do it later.
	if (Planet(4))
	{
		do_mars(T0, planet_info);
	}
	if (Planet(5))
	{
		do_jupiter(T0, planet_info);
	}
	if (Planet(6))
	{
		do_saturn(T0, planet_info);
	}
	if (Planet(7))
	{
		do_uranus(T0, planet_info);
	}
	if (Planet(8))
	{
		do_neptune(T0, planet_info);
	}
	/* Invoke the special periodic term approximation for Pluto, segregated
	   in another module to preserve the purity of this file which is based
	   on genuine theories of planetary motion valid for millenia. */
    if (Planet(9)) {
		pluto(jd, &plan_ll, &plan_b, &plan_r);
		trans(9, plan_r, plan_b, plan_ll, plan_Stheta, plan_Sr, plan_epli, MAGPLU, planet_info);
	}
	/* Fill in the information for the Moon in slot number 3. */
    if (Planet(3)) {
		highmoon(jd, &plan_ll, &plan_b);
		ecliptoeq(jd, plan_ll, plan_b, &planet_info[3].ra, &planet_info[3].dec);
	}
}

void Starmap::calcPlanet(double jd)
{
	int i;
	double igmst, latsin, latcos;

	planets(jd, 0xFFFF, planet_info);
	//igmst = gmst(jd);
	rgmst(jd, &igmst);
	latsin = sin(dtr(siteLat));
	latcos = cos(dtr(siteLat));
	for (i = 0; i <= 9; i++) {
		planet_info[i].lha = dtr(fixangle((igmst * 15) - siteLon - planet_info[i].ra));
		planet_info[i].az = rtd(atan2(sin(planet_info[i].lha), cos(planet_info[i].lha) * latsin -
								tan(dtr(planet_info[i].dec)) * latcos));
		planet_info[i].alt = rtd(asin(latsin * sin(dtr(planet_info[i].dec)) +
								latcos * cos(dtr(planet_info[i].dec)) * cos(planet_info[i].lha)));
	}
}

void Starmap::plotLine(double fdec, double fra, double tdec, double tra)
{
	int vx, vy, vx2, vy2;
	double cdec1, cra1, cdec2, cra2;

	if (clipr_xform(fdec, fra, tdec, tra, &vx, &vy, &vx2, &vy2, FALSE,
		&cdec1, &cra1, &cdec2, &cra2)) {
		sprintf(log2ram_buf, "drawing line\n");
		drawcurveline(cdec1, cra1, cdec2, cra2, vx, vy, vx2, vy2, 0, FALSE, 0);
	}
}

void Starmap::paintSky(double limag, rect_s* br)
{
	int i, vx, vy, vis, vx2, vy2, smhand = -1, pflip, precess;
	void *smap;
	double igmst, lham, ra, dec, mag;
	unsigned char smex[4];
	int ii;
	
	int skyTel=TRUE; // always true
	
//	jdtime=0;

mydeb=0xa1;    
	if ((precessionCalculation == PrecAlways) ||
		((precessionCalculation == PrecAuto) && (abs(PrecEpoch - jdtime) >
			(PrecYears * 365.25)))) {
		definePrecession(2000.0 + ((jdtime - PrecEpoch) / 365.25));
		precess = TRUE;
	} else {
		precess = FALSE;
	}
mydeb=0xa2;
#define Prd(x, y) if (precess) { precessObject(x, y, &x, &y); }

	//igmst = gmst(jdtime);
	rgmst(jdtime, &igmst);
mydeb=0xa3;
    lham = skyLham = dtr((igmst * 15) - siteLon);
mydeb=0xa4;
	skywin.width = br->right - br->left;
	skywin.height = br->bottom - br->top;
	skywin.x_offset = br->left;
	skywin.y_offset = br->top;
	skywin.proj_mode = STEREOGR;
	pflip = skywin.invert = Flip > 0;
	//igmst = gmst(jdtime);
	rgmst(jdtime, &igmst);
mydeb=0xa5;
	skywin.racen = fixangle((igmst * 15) - siteLon);
mydeb=0xa6;
	skywin.dlcen = siteLat;
	skywin.scale = 90.0;
	initxform(&skywin);
mydeb=0xa7;
	// Draw coordinate grid markers, if requested

	if (skyShowCoords) {
		int i;
		int cpen, open;
		double epsilon, esin, ecos, eqra, eqdec, eqlat, pera, pedec;
		int ofont, cfont;
#define tickWid	0.75
mydeb=0xa8;
        cfont = 0; //CreateFont;
        ofont = cfont; //SelectObject(hDC, cfont);
    	//SetBkMode(hDC, TRANSPARENT);
    	//SetTextAlign(hDC, TA_CENTER | TA_BASELINE | TA_NOUPDATECP);
        //SetTextColor(hDC, RGB(0, 128, 128));

		open = 0;
		cpen=0;
		 //SelectObject(hDC, (cpen = CreatePen(PS_SOLID, 0, RGB(0, 128, 128))));
		set_col(col_coord_grid);
		for (i = -1; i <= 1; i += 2) {
			// Equinoctual colures
			plotLine(88.0 * i, 0.0, 90.0 * i, 0.0);
			plotLine(88.0 * i, 180.0, 90.0 * i, 180.0);

			// Solstitial colures
			plotLine(88.0 * i, 90.0, 90.0 * i, 90.0);
			plotLine(88.0 * i, 270.0, 90.0 * i, 270.0);
		}
mydeb=0x11;		
		// Celestial equator
		for (i = 0; i < 360; i += 15) {		// Draw the tick marks at hour angle intervals
			char hlab[4];
			int vx, vy, vis;

            plotLine(0.0, (double) i, 0.0, i + 15.0);
			plotLine(-tickWid, (double) i, tickWid, (double) i);
			sprintf(hlab, "%d'", i / 15);  // Format 42
			xform(tickWid * (pflip ? 1 : -1), (double) i, &vx, &vy, &vis);
			if (vis) {
				TextOut(vx, vy, hlab, strlen(hlab), TEXT_TYPE_CELEST_EQ);
			}
mydeb=mydeb+1;

		}
mydeb=0x70;
		//SelectObject(hDC, open);
		//DeleteObject(cpen);
		//ticker();

		// Ecliptic

		open = cpen; //SelectObject(hDC, (cpen = CreatePen(PS_SOLID, 0, RGB(128, 0, 0))));
		set_col(col_ecliptic);
mydeb=0x71;
        robliqeq(jdtime, &epsilon);
        epsilon=dtr(epsilon);
//		epsilon = dtr(obliqeq(jdtime));	// Get current obliquity of ecliptic
mydeb=0x72;
		esin = sin(epsilon);
mydeb=0x73;
		ecos = cos(epsilon);
mydeb=0x74;
        //SetTextColor(hDC, RGB(128, 0, 0));
		pera = pedec = 0.0;		// Dirty trick: ecliptic intersects equator at 0 longitude !
mydeb=0x75;
		for (i = 1; i <= 32; i++) {			// Draw the ecliptic itself
			eqlat = ((PI * 2) / 32.0) * i;
			eqra = fixangle(rtd(atan2(ecos * sin(eqlat), cos(eqlat))));
			eqdec = rtd(asin(esin * sin(eqlat)));
			plotLine(pedec, pera, eqdec, eqra);
			pera = eqra;
			pedec = eqdec;
mydeb++;
		}
mydeb=0x20;
		for (i = 0; i < 360; i += 15) {		// Draw the tick marks at 15 degree intervals
			char hlab[6];
			int vx, vy, vis;

			eqlat = ((PI * 2) / 360.0) * i;
			pera = fixangle(rtd(atan2((ecos * sin(eqlat) -
					     	(tan(dtr(-tickWid)) * esin)), cos(eqlat))));
			pedec = rtd(asin((esin * sin(eqlat) * cos(dtr(-tickWid))) +
			     		(sin(dtr(-tickWid)) * ecos)));
			eqra = fixangle(rtd(atan2((ecos * sin(eqlat) -
					     	(tan(dtr(tickWid)) * esin)), cos(eqlat))));
			eqdec = rtd(asin((esin * sin(eqlat) * cos(dtr(tickWid))) +
			     		(sin(dtr(tickWid)) * ecos)));
			plotLine(pedec, pera, eqdec, eqra);
			sprintf(hlab, "%d''", i); // Format(43)
			xform(!pflip ? pedec : eqdec, !pflip ? pera : eqra, &vx, &vy, &vis);
			if (vis) {
				TextOut(vx, vy, hlab, strlen(hlab), TEXT_TYPE_ECLIPTIC);
			}

		}
mydeb=0x21;		
#undef tickWid

		//SelectObject(hDC, open);
		//DeleteObject(cpen);
		//SelectObject(hDC, ofont);
		//DeleteObject(cfont);
	}
	//ticker();
	//if (bailout()) {  // this is the check for keyboard events, to bail out of calculations early
	//	return;
	//}

	// Draw constellation names, if requested

	if (do_constellation_text) 
	{
		int cntick = 0;		
		char *vfmt;// = Format(40);

		if (1) { 
			char cl[80], cn[80];
			unsigned int tra;
			int tdec;
			int ofont, cfont;
			unsigned char *const_ptr;
			char constel[20];
			char constel_eng[20];
			
			const_ptr=(unsigned char *)constellation_array;

			if (!skyAlignConnames) {
	            cfont = 0; // normal font?
	            ofont = cfont; //SelectObject(hDC, cfont);
            }
        	//SetBkMode(hDC, TRANSPARENT);
           	//SetTextColor(hDC, RGB(255, 255, 0));
        	//SetTextAlign(hDC, TA_CENTER | TA_BASELINE | TA_NOUPDATECP);
			for (i=0; i<89; i++) // 89 constellation names
			{			
			    tra = (*const_ptr++)<<24;
			    tra |= ((*const_ptr++)<<16);
			    tra |= ((*const_ptr++)<<8);
			    tra |= *const_ptr++;
			    
			    tdec = (*const_ptr++)<<24;
			    tdec |= ((*const_ptr++)<<16);
			    tdec |= ((*const_ptr++)<<8);
			    tdec |= *const_ptr++;
			    
			    
			    strncpy(constel, (const char*)const_ptr, 20);
			    const_ptr=const_ptr+20;
			    strncpy(constel_eng, (const char*)const_ptr, 20);
			    const_ptr=const_ptr+20;
				
				strncpy(cn, constel, 20); // use normal constellation name (not english name) for now
				
				ra = tra / (24000.0 / 360.0);
	    		dec = tdec / (9000.0 / 90.0);
	    		Prd(ra, dec);
				xform(dec, ra, &vx, &vy, &vis);
				if (vis)
				{
					if (1) //skyAlignConnames
					{
						int iangle;
						double eangle;
	
						vx2 = vx - (br->left + (br->right - br->left) / 2);
						vy2 = vy - (br->top + (br->bottom - br->top) / 2);
						if (vx2 == vy2)
						{
							eangle = 0;
						}
						else
						{
							eangle = rtd(atan2((double) vy2, (double) vx2));
							if (eangle >= 0)
							{
								eangle -= 90;
							}
							else
							{
								eangle += 90;
							}
						}
						eangle = - eangle;
						iangle = (int) (eangle * 10); 
			            cfont = 0; //CreateFont
			            ofont = 0; //SelectObject(hDC, cfont);
		            }
					TextOut(vx, vy, cn, strlen(cn), TEXT_TYPE_CONSTEL);
					if (1) //skyAlignConnames
					{
						//SelectObject(hDC, ofont);
						//DeleteObject(cfont);
					}
					if ((++cntick % 2) == 0)
					{
						//ticker();
						//if (bailout()) { this is for keyboard events
							//break;
						//}
					}
				}
			}
			//closeRSC();
			//if (!skyAlignConnames) {			
			//	SelectObject(hDC, ofont);
			//	DeleteObject(cfont);
			//}
		}
	}
	//ticker();
	//if (bailout()) { keyboard events
	//	return;
	//}

	// Draw constellation outlines, if requested

	if (skyShowConstellations) {
		int conoTick = 0;		
		char *vfmt; // = Format(36);
		short* cl_ptr;
		unsigned char* txt_ptr;
		unsigned short *ra_ptr;
		char cl[80], cn[10];
	    unsigned int fra, tra;
	    int fdec, tdec;
		double ra2, dec2, cra1, cdec1, cra2, cdec2;
		int cpen, open;
		cl_ptr = (short*)constellation_lines_array;

		if (1)
		{


            cl[3]=0; // end of string
			open = 0; //SelectObject(hDC, (cpen = CreatePen(PS_SOLID, 0, RGB(128, 128, 128))));
			set_col(col_constel);
			for (i=0; i<646; i++) // 646 lines
			{
				//cl[strlen(cl) - 1] = EOS;
				txt_ptr=(unsigned char *)cl_ptr;
				cl[0] = *txt_ptr++;
				cl[1] = *txt_ptr++;
				cl[2] = *txt_ptr;
				cl_ptr += 2; // move pointer along by 4 bytes (2 shorts)
				
				
				ra_ptr=(unsigned short *)cl_ptr++;
				fra=(unsigned int)*ra_ptr;
				fdec=(int)*cl_ptr++;
				ra_ptr=(unsigned short *)cl_ptr++;
				tra=(unsigned int)*ra_ptr;
				tdec=(int)*cl_ptr++;
				

	    		ra = fra / (24000.0 / 360.0);
	    		dec = fdec / (9000.0 / 90.0);
	    		Prd(ra, dec);
	    		ra2 = tra / (24000.0 / 360.0);
	    		dec2 = tdec / (9000.0 / 90.0);
	    		Prd(ra2, dec2);
				if (clipr_xform(dec, ra, dec2, ra2, &vx, &vy, &vx2, &vy2, FALSE,
					&cdec1, &cra1, &cdec2, &cra2)) {
					MoveTo(vx, vy);
					LineTo(vx2, vy2);
				}
			}
			//SelectObject(hDC, open);
			//DeleteObject(cpen);
			//closeRSC();
		}
	}

	// *** plot stars ****
	// smap = LockResource(starMapRes);
	uint32_t addr = ADDR_YALE;
	if (1) { //(smap != NULL) {
		struct starMapData sd[7]; // we misuse this to store data read from storage
		struct starMapData *smp; //  = (smap == NULL) ? &sd : smap;
		int imag;
		short ipmra, ipmdec;
		unsigned int ilimag;
		int starTick = 0;
		char* test_end;
		int read_result;
		
		// struct is 7 bytes, but there are optional fields for smex, and impra/impdec and also the star name.
		// we should read at least 7+3+4 bytes (14) plus the length of star name, lets assume 32 bytes for that,
		// so total is 46 bytes, which should fit in sd[7] since the struct is of size 7 bytes.
		read_result = storage_read(addr, (char*)sd, 46);
		smp = &sd[0];

		// if (fChildPreview) {
		//	limag =  2.0;
		//}
        ilimag = (unsigned int) ((limag + 1.5) * 100); // For quick tests against limiting magnitude
		while (read_result>=0) {
			char *sname; // used as a pointer to after the struct

			// check if end of star database
			test_end = (char*)&sd;
			if ((*test_end=='Y') && (*(test_end+1)=='Y') && (*(test_end+2)=='Y')) {
				break;
			}
			
			// smex is 3 bytes after the struct, so read that if bitmap indicates it is present
          	memset(smex, 0, sizeof smex);
          	sname = ((char *) smp) + 7; //sizeof(struct starMapData);
          	if (smp->mag & 0x2000) {
          		memcpy(smex, sname, 3);
          		sname += 3;
          	}
			// read impra and impdec (total 4 bytes) if bitmap indicates they are present
          	if (smp->mag & 0x1000) {
          		memcpy(&ipmra, sname, sizeof ipmra);
          		memcpy(&ipmdec, sname + sizeof ipmra, sizeof ipmdec);
          		sname += sizeof ipmra + sizeof ipmdec;
          	} else {
          		ipmra = ipmdec = 0;
          	}
    		ra = smp->lon / (65536.0 / 360.0) +
    			(((jdtime - PrecEpoch) / (JulianCentury / 100)) * (ipmra / (1000.0 * 60 * 60)));
			dec = smp->lat / (32767.0 / 90.0) +
				(((jdtime - PrecEpoch) / (JulianCentury / 100)) * (ipmdec / (1000.0 * 60 * 60)));
    		Prd(ra, dec);

			xform(dec, ra, &vx, &vy, &vis);
			if (vis) {
				double mbase = 3.5;
                
                if ((++starTick % 40) == 0) {
                	starTick = 0;
                	//ticker();
                	//if (bailout()) {
                	//	break;
                	//}
                } 
	    		mag = ((smp->mag & 0xFFF) / 100.0) - 1.5;
				if (1) { //!fChildPreview) {
					int named = FALSE;

					/* If star is bright enough to meet the "show name" criterion
					   and it has a name, display it. */

					if ((smp->mag & 0x8000) && skyShowName && (mag < skyNameMag)) {
						//SelectObject(hDC, cfont);
	       		     	TextOut(vx + 8, vy + 2, sname, strlen(sname), TEXT_TYPE_BRIGHT_ST);
	       		     	named = TRUE;
					}

					// See if the Bayer/Flamsteed code should be drawn

					if ((smex[0] != 0) && skyShowBflam && (mag < skyBflamMag)) {
						//SelectObject(hDC, smex[0] < 32 ? gfont : cfont);
						if (smex[0] < 32) {
							smex[0] += ('a' - 1);
						}
						TextOut(vx + (named ? -16 : 8), vy + 2,
							(char*)smex, strlen((char*)smex), TEXT_TYPE_BAYERF_CODE);

					}
				}

				// Finally, paint the star according to its magnitude and the quality setting
#ifdef WITH_ICONS
				imag = ((int) ELEMENTS(starIcons)) - ((int) (mag - BrightestMag));
				if (limag < (ELEMENTS(starIcons) - 2)) {
					imag -= (ELEMENTS(starIcons) - ((int) limag)) - 1;
				}
				imag = min(max(imag, 0), ELEMENTS(starIcons) - 1);
				DrawIcon(hDC, vx - 16, vy - 16, starIcons[imag]);
#else
				// handle 9 different star magnitudes
				imag = 9 - ((int) (mag - BrightestMag));
				if (limag < 7) {
					imag -= (9-((int)limag)) -1;
				}
				imag = imag -1;
				if (imag<0) imag = 0;
				if (imag>8) imag = 8;
				DrawStar(vx, vy, imag);
				// draw a cross for now
				//draw_line(vx-imag, vy-imag, vx+imag, vy+imag, DEFAULT_COL_STAR);
				//draw_line(vx-imag, vy+imag, vx+imag, vy-imag, DEFAULT_COL_STAR);
#endif 
			}

			/* Quit if we hit the end of file or exceed the limiting magnitude.
			   Note the assumption that the catalogue is sorted by magnitude! */

    		if ((smp->mag & 0x4000) || ((smp->mag & 0xFFF) > ((int) ilimag))) {
    			break;
    		}

			// increment addr to the next star entry
			addr += (sname - ((char *) smp)) + ((smp->mag & 0x8000) ? strlen(sname) + 1 : 0);
			// read from storage
			read_result = storage_read(addr, (char*)sd, 46);
			smp = &sd[0];

    		//if (smap != NULL) {
    		//	smp = (struct starMapData *) (sname + ((smp->mag & 0x8000) ? (strlen(sname) + 1) : 0));
    		//}
		}


	}

	// *** end of plot stars ****

	// calc planet positions
	calcPlanet(jdtime);
	// get moon phase
	double p, aom, cphase, cdist, cangdia, csund, csuang;
	p = phase(jdtime, &cphase, &aom, &cdist, &cangdia, &csund, &csuang);
    int mage = (int) aom;
	int ph = (int)(cphase*100);
	//printf("phase cphase: %d\n", ph);
	

    if (skyShowPlanets) {
		for (i = 0; i <= 9; i++) {
			if (((planet_info[i].alt > 0.0) || !skyTel) && ((i < 9) || (planet_info[i].hrv > 0))) {
	    		ra = dtr(planet_info[i].ra);
	    		dec = dtr(planet_info[i].dec);

	    		/* Note that since planetary positions are computed for the equinox
	    		   of the date, they are *not* precessed to the current equinox. */

				xform(rtd(dec), rtd(ra), &vx, &vy, &vis);
				if (!vis) {
					continue;
				}

				if (i == 3) {
					// Use the moon icon to show current phase
					//DrawIcon(hDC, vx - 16, vy - 13, micon);
					DrawMoon(vx, vy, mage);
				} else {
					// DrawIcon(hDC, vx - 16, vy - 16, PlanetIcon(i));
				}
			}
		}
	}




	//ticker();
	//if (bailout()) {
	//	return;
	//}
  //********************************
  return;
}

void Starmap::draw_line(int x0, int y0, int x1, int y1, uint16_t color)
{
    uint16_t pix = color; // could put a colour mapping function here
    int dx, dy;
    int stepx, stepy;
    int length, extras, incr2;
    int c, incr1, d, i;
    
    dy = y1 - y0;
    dx = x1 - x0;
    
    if (dy < 0) { dy = -dy;  stepy = -1; } else { stepy = 1; }
    if (dx < 0) { dx = -dx;  stepx = -1; } else { stepx = 1; }

    plot_pixel(pix, x0, y0);
    plot_pixel(pix, x1, y1);
    if (dx > dy) {
        length = (dx - 1) >> 2;
        extras = (dx - 1) & 3;
        incr2 = (dy << 2) - (dx << 1);
        if (incr2 < 0) {
            c = dy << 1;
            incr1 = c << 1;
            d =  incr1 - dx;
            for (i = 0; i < length; i++) {
                x0 += stepx;
                x1 -= stepx;
                if (d < 0) {						// Pattern:
                    plot_pixel(pix, x0, y0);			//
                    plot_pixel(pix, x0 += stepx, y0);	//  x o o
                    plot_pixel(pix, x1, y1);			//
                    plot_pixel(pix, x1 -= stepx, y1);
                    d += incr1;
                } else {
                    if (d < c) {							// Pattern:
                        plot_pixel(pix, x0, y0);				//      o
                        plot_pixel(pix, x0 += stepx, y0 += stepy);		//  x o
                        plot_pixel(pix, x1, y1);				//
                        plot_pixel(pix, x1 -= stepx, y1 -= stepy);
                    } else {
                        plot_pixel(pix, x0, y0 += stepy);			// Pattern:
                        plot_pixel(pix, x0 += stepx, y0);			//    o o 
                        plot_pixel(pix, x1, y1 -= stepy);			//  x
                        plot_pixel(pix, x1 -= stepx, y1);			//
                    }
                    d += incr2;
                }
            }
            if (extras > 0) {
                if (d < 0) {
                    plot_pixel(pix, x0 += stepx, y0);
                    if (extras > 1) plot_pixel(pix, x0 += stepx, y0);
                    if (extras > 2) plot_pixel(pix, x1 -= stepx, y1);
                } else
                if (d < c) {
                    plot_pixel(pix, x0 += stepx, y0);
                    if (extras > 1) plot_pixel(pix, x0 += stepx, y0 += stepy);
                    if (extras > 2) plot_pixel(pix, x1 -= stepx, y1);
                } else {
                    plot_pixel(pix, x0 += stepx, y0 += stepy);
                    if (extras > 1) plot_pixel(pix, x0 += stepx, y0);
                    if (extras > 2) plot_pixel(pix, x1 -= stepx, y1 -= stepy);
                }
            }
        } else {
            c = (dy - dx) << 1;
            incr1 = c << 1;
            d =  incr1 + dx;
            for (i = 0; i < length; i++) {
                x0 += stepx;
                x1 -= stepx;
                if (d > 0) {
                    plot_pixel(pix, x0, y0 += stepy);			// Pattern:
                    plot_pixel(pix, x0 += stepx, y0 += stepy);		//      o
                    plot_pixel(pix, x1, y1 -= stepy);			//    o
                    plot_pixel(pix, x1 -= stepx, y1 -= stepy);		//  x
                    d += incr1;
                } else {
                    if (d < c) {
                        plot_pixel(pix, x0, y0);				// Pattern:
                        plot_pixel(pix, x0 += stepx, y0 += stepy);       //      o
                        plot_pixel(pix, x1, y1);                         //  x o
                        plot_pixel(pix, x1 -= stepx, y1 -= stepy);       //
                    } else {
                        plot_pixel(pix, x0, y0 += stepy);			// Pattern:
                        plot_pixel(pix, x0 += stepx, y0);			//    o o
                        plot_pixel(pix, x1, y1 -= stepy);			//  x
                        plot_pixel(pix, x1 -= stepx, y1);			//
                    }
                    d += incr2;
                }
            }
            if (extras > 0) {
                if (d > 0) {
                    plot_pixel(pix, x0 += stepx, y0 += stepy);
                    if (extras > 1) plot_pixel(pix, x0 += stepx, y0 += stepy);
                    if (extras > 2) plot_pixel(pix, x1 -= stepx, y1 -= stepy);
                } else
                if (d < c) {
                    plot_pixel(pix, x0 += stepx, y0);
                    if (extras > 1) plot_pixel(pix, x0 += stepx, y0 += stepy);
                    if (extras > 2) plot_pixel(pix, x1 -= stepx, y1);
                } else {
                    plot_pixel(pix, x0 += stepx, y0 += stepy);
                    if (extras > 1) plot_pixel(pix, x0 += stepx, y0);
                    if (extras > 2) {
                        if (d > c)
                            plot_pixel(pix, x1 -= stepx, y1 -= stepy);
                        else
                            plot_pixel(pix, x1 -= stepx, y1);
                    }
                }
            }
        }
    } else {
        length = (dy - 1) >> 2;
        extras = (dy - 1) & 3;
        incr2 = (dx << 2) - (dy << 1);
        if (incr2 < 0) {
            c = dx << 1;
            incr1 = c << 1;
            d =  incr1 - dy;
            for (i = 0; i < length; i++) {
                y0 += stepy;
                y1 -= stepy;
                if (d < 0) {
                    plot_pixel(pix, x0, y0);
                    plot_pixel(pix, x0, y0 += stepy);
                    plot_pixel(pix, x1, y1);
                    plot_pixel(pix, x1, y1 -= stepy);
                    d += incr1;
                } else {
                    if (d < c) {
                        plot_pixel(pix, x0, y0);
                        plot_pixel(pix, x0 += stepx, y0 += stepy);
                        plot_pixel(pix, x1, y1);
                        plot_pixel(pix, x1 -= stepx, y1 -= stepy);
                    } else {
                        plot_pixel(pix, x0 += stepx, y0);
                        plot_pixel(pix, x0, y0 += stepy);
                        plot_pixel(pix, x1 -= stepx, y1);
                        plot_pixel(pix, x1, y1 -= stepy);
                    }
                    d += incr2;
                }
            }
            if (extras > 0) {
                if (d < 0) {
                    plot_pixel(pix, x0, y0 += stepy);
                    if (extras > 1) plot_pixel(pix, x0, y0 += stepy);
                    if (extras > 2) plot_pixel(pix, x1, y1 -= stepy);
                } else
                if (d < c) {
                    plot_pixel(pix, stepx, y0 += stepy);
                    if (extras > 1) plot_pixel(pix, x0 += stepx, y0 += stepy);
                    if (extras > 2) plot_pixel(pix, x1, y1 -= stepy);
                } else {
                    plot_pixel(pix, x0 += stepx, y0 += stepy);
                    if (extras > 1) plot_pixel(pix, x0, y0 += stepy);
                    if (extras > 2) plot_pixel(pix, x1 -= stepx, y1 -= stepy);
                }
            }
        } else {
            c = (dx - dy) << 1;
            incr1 = c << 1;
            d =  incr1 + dy;
            for (i = 0; i < length; i++) {
                y0 += stepy;
                y1 -= stepy;
                if (d > 0) {
                    plot_pixel(pix, x0 += stepx, y0);
                    plot_pixel(pix, x0 += stepx, y0 += stepy);
                    plot_pixel(pix, x1 -= stepy, y1);
                    plot_pixel(pix, x1 -= stepx, y1 -= stepy);
                    d += incr1;
                } else {
                    if (d < c) {
                        plot_pixel(pix, x0, y0);
                        plot_pixel(pix, x0 += stepx, y0 += stepy);
                        plot_pixel(pix, x1, y1);
                        plot_pixel(pix, x1 -= stepx, y1 -= stepy);
                    } else {
                        plot_pixel(pix, x0 += stepx, y0);
                        plot_pixel(pix, x0, y0 += stepy);
                        plot_pixel(pix, x1 -= stepx, y1);
                        plot_pixel(pix, x1, y1 -= stepy);
                    }
                    d += incr2;
                }
            }
            if (extras > 0) {
                if (d > 0) {
                    plot_pixel(pix, x0 += stepx, y0 += stepy);
                    if (extras > 1) plot_pixel(pix, x0 += stepx, y0 += stepy);
                    if (extras > 2) plot_pixel(pix, x1 -= stepx, y1 -= stepy);
                } else
                if (d < c) {
                    plot_pixel(pix, x0, y0 += stepy);
                    if (extras > 1) plot_pixel(pix, x0 += stepx, y0 += stepy);
                    if (extras > 2) plot_pixel(pix, x1, y1 -= stepy);
                } else {
                    plot_pixel(pix, x0 += stepx, y0 += stepy);
                    if (extras > 1) plot_pixel(pix, x0, y0 += stepy);
                    if (extras > 2) {
                        if (d > c)
                            plot_pixel(pix, x1 -= stepx, y1 -= stepy);
                        else
                            plot_pixel(pix, x1, y1 -= stepy);
                    }
                }
            }
        }
    }
}

void Starmap::text_out(int x, int y, char* lab, unsigned char len, char type)
{
    int i;
	int col;

	switch (type) {
		case TEXT_TYPE_CONSTEL:
			col = col_constel_text;
			break;
		case TEXT_TYPE_BRIGHT_ST:
			col = col_bright_st_text;
			break;
		case TEXT_TYPE_BAYERF_CODE:
			col = col_bayerf_text;
			break;
		case TEXT_TYPE_CELEST_EQ:
			col = col_celest_eq_text;
			break;
		case TEXT_TYPE_ECLIPTIC:
			col = col_ecliptic_text;
			break;
		default:
			col = DEFAULT_COL_TEXT_GENERIC;
			break;
	}

	for(i=0; i<len; i++) {
		plot_char(lab[i], x, y, col);
		x = x + 6;
	}
}

void Starmap::MoveTo(int x, int y)
{
    draw_x_store=x;
    draw_y_store=y;
}

void Starmap::LineTo(int x, int y)
{
    draw_line(draw_x_store, draw_y_store, x, y, draw_col_store);
}

void Starmap::set_col(uint16_t color)
{
    draw_col_store=color;
}

void Starmap::TextOut(int x, int y, char* lab, unsigned char len, char type)
{
    text_out(x, y, lab, len, type);
}

void Starmap::plot_pixel(uint16_t color, int x, int y)
{
    // to do
    if (x>SCREEN_W) return;
    if (x<0) return;
    if (y>SCREEN_H) return;
    if (y<0) return;
    
    if (color)
        display_ram[y][x]=0xff;
    else
        display_ram[y][x]=0x00;
    
    
}

void Starmap::plot_char(char c, int x, int y, int color)
{
    int i, j;
    // check that the character is in the font
    if (c < 0x20 || c > 0x7f)
    {
        return;
    }

    for (i = 0; i < 5; i++)
    {
        for (j = 0; j < 7; j++)
        {
            if (font5_7[c - 0x20][i] & (1 << j))
            {
                plot_pixel(color, x + i, y - j);
            }
        }
    }
}

int Starmap::OnScreen(int x, int y)
{
	if ((x<0) || (y<0)) return(-1);
	if ((x>=skywin.width) || (y>=skywin.height)) return(-1);
	return(1);
}

void Starmap::DrawMoon(int x, int y, int age)
{
	int phase;
	int i, j;
	
	phase = (age + 2) % 30;
    phase = (phase * 8) / 30;

    for (i=0; i<7; i++) {
        for (j=0; j<7; j++) {
            if (moonphtext[phase][i] & (1 << j)) {
                plot_pixel(col_moon_phtext, x+i-3, y-j+3);
            } else {
                if (moonbright[i] & (1 << j)) {
                    plot_pixel(col_moon_bright, x+i-3, y-j+3);
                }
                if (moondim[i] & (1 << j)) {
                    plot_pixel(col_moon_dim, x+i-3, y-j+3);
                }
                if (moondark[i] & (1 << j)) {
                    plot_pixel(col_moon_dark, x+i-3, y-j+3);
                }
            }
        }
    }
}

void Starmap::DrawStar(int x, int y, int level)
{
	int col;
	// plot the center pixel, for level 0 and 1
	if (OnScreen(x, y)) {
		if (level==0) {
			plot_pixel(col_stardim, x, y);
		} else {
			plot_pixel(col_starbright, x, y);
		}
	}
	if (level <= 1) return;
	// plot the pixels for levels 2 and 3
	if (level==2) {
		col = col_stardim;
	} else {
		col = col_starbright;
	}
	if (OnScreen(x, y+1)) plot_pixel(col, x, y+1);
	if (OnScreen(x, y-1)) plot_pixel(col, x, y-1);
	if (OnScreen(x+1, y)) plot_pixel(col, x+1, y);
	if (OnScreen(x-1, y)) plot_pixel(col, x-1, y);
	if (level <= 3) return;
	// plot the pixels for levels 4 and 5
	if (level==4) {
		col = col_stardim;
	} else {
		col = col_starbright;
	}
	if (OnScreen(x-1, y-1)) plot_pixel(col, x-1, y-1);
	if (OnScreen(x+1, y-1)) plot_pixel(col, x+1, y-1);
	if (OnScreen(x-1, y+1)) plot_pixel(col, x-1, y+1);
	if (OnScreen(x+1, y+1)) plot_pixel(col, x+1, y+1);
	if (level <= 5) return;
	// plot the pixels for levels 6 and 7
	if (level==6) {
		col = col_stardim;
	} else {
		col = col_starbright;
	}
	if (OnScreen(x, y-2)) plot_pixel(col, x, y-2);
	if (OnScreen(x, y+2)) plot_pixel(col, x, y+2);
	if (OnScreen(x+2, y)) plot_pixel(col, x+2, y);
	if (OnScreen(x-2, y)) plot_pixel(col, x-2, y);
	if (level <= 7) return;
	// plot the pixels for level 9
	col = col_stardim;
	if (OnScreen(x-2, y-1)) plot_pixel(col, x-2, y-1);
	if (OnScreen(x-2, y+1)) plot_pixel(col, x-2, y+1);
	if (OnScreen(x+2, y-1)) plot_pixel(col, x+2, y-1);
	if (OnScreen(x+2, y+1)) plot_pixel(col, x+2, y+1);

	if (OnScreen(x+1, y-2)) plot_pixel(col, x+1, y-2);
	if (OnScreen(x+1, y+2)) plot_pixel(col, x+1, y+2);
	if (OnScreen(x-1, y-2)) plot_pixel(col, x-1, y-2);
	if (OnScreen(x-1, y+2)) plot_pixel(col, x-1, y+2);
}

int Starmap::storage_read(uint32_t addr, char* data, uint16_t len)
{
	return(-1);
}


// *************************************************
// ****   End of Starmap class implementation  ****
// *************************************************

