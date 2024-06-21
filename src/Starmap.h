// Starmap.h
// rev 1 - May 2024 - shabaz - first revision

#ifndef STARMAP_H
#define STARMAP_H

#include "Arduino.h"
#include <math.h>
#include "arrays.h"


// defines
/* storage addresses */
#define ADDR_ID 0x00000
#define ADDR_CONSTEL_ARRAY 0x02000
#define ADDR_CONSTEL_LINES_ARRAY 0x03400
#define ADDR_CONSTEL_BOUND_ARRAY 0x05400
#define ADDR_YALE 0x08c00
#define ADDR_YALE_END 0x24259
/* colors */
#define DEFAULT_COL_COORD_GRID 0xff00
#define DEFAULT_COL_ECLIPTIC   0x00ff
#define DEFAULT_COL_CONSTEL    0xffff
#define DEFAULT_COL_STARDIM    0x03ff
#define DEFAULT_COL_STARBRIGHT 0x0fff
#define DEFAULT_COL_STARTEXT 0xffff
#define DEFAULT_COL_MOON_BRIGHT 0xff11
#define DEFAULT_COL_MOON_DIM 0xff22
#define DEFAULT_COL_MOON_DARK 0x01ff
#define DEFAULT_COL_MOON_PHTEXT 0x00ff
#define DEFAULT_COL_TEXT_GENERIC 0x00ff
/*  Astronomical constants  */
#define epoch	    2444238.5	   /* 1980 January 0.0 */
/*  Constants defining the Sun's apparent orbit  */
#define elonge	    278.833540	   /* Ecliptic longitude of the Sun at epoch 1980.0 */
#define elongp	    282.596403	   /* Ecliptic longitude of the Sun at perigee */
#define eccent      0.016718       /* Eccentricity of Earth's orbit */
#define sunangsiz   0.533128       /* Sun's angular size, degrees, at semi-major axis distance */
/*  Elements of the Moon's orbit, epoch 1980.0  */
#define mmlong      64.975464      /* Moon's mean longitude at the epoch */
#define mmlongp     349.383063	   /* Mean longitude of the perigee at the epoch */
#define mlnode	    151.950429	   /* Mean longitude of the node at the epoch */
#define minc        5.145396       /* Inclination of the Moon's orbit */
#define mecc        0.054900       /* Eccentricity of the Moon's orbit */
#define mangsiz     0.5181         /* Moon's angular size at distance a from Earth */
#define msmax       384401.0       /* Semi-major axis of Moon's orbit in km */
#define mparallax   0.9507	   /* Parallax at distance a from Earth */
#define dsin(x) (sin(dtr((x)))) 			/* Sin from deg */
#define dcos(x) (cos(dtr((x)))) 			/* Cos from deg */
#define DCOS(x) (cos((x) * .0174532925199))
#define DSIN(x) (sin((x) * .0174532925199))
#define DTAN(x) (tan((x) * .0174532925199))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define Fuz 0.1
#define Planet(x)	(which && (1 << (x)))
#define BrightestMag -1

#define TRUE  1
#define FALSE 0
// screen dimensions, was 160x160, set to 4x4 (unusable),
// because we expect the coder to have their own plot_pixel function,
// so display_ram[] is not used.
#define SCREEN_W 4
#define SCREEN_H 4
//	Frequently used astronomical constants
#define J2000				2451545.0		// Julian day of J2000 epoch
#define JulianCentury		36525.0			// Days in Julian century
#define AstronomicalUnit	149597870.0		// Astronomical unit in kilometres
#define SunSMAX	 (AstronomicalUnit * 1.000001018) // Semi-major axis of Earth's orbit
#define EarthRad			6378.14			// Earth's equatorial radius, km (IAU 1976)
#define LunatBase			2423436.0		/* Base date for E. W. Brown's numbered series of lunations (1923 January 16) */
#define SynMonth			29.53058868		// Synodic month (mean time from new Moon to new Moon)
//	Precession calculation modes
#define PrecAuto		0					// Precess if more then PrecYears from PrecEpoch
#define PrecAlways		1					// Always precess
#define PrecNever		2					// Never correct for precession
#define PrecEpoch		J2000				// All databases are epoch J2000.0
#define PrecYears		25					// Consider databases valid for this time around epoch
//
#ifndef PI
#define PI 3.14159265358979323846
#endif
#define sgn(x) (((x) < 0) ? -1 : ((x) > 0 ? 1 : 0))       /* Extract sign */
#ifndef abs
#define abs(x) ((x) < 0 ? (-(x)) : (x))                   /* Absolute val */
#endif
#define fixangle(a) ((a) - 360.0 * (floor((a) / 360.0)))  /* Fix angle    */
#define fixangr(a)  ((a) - (PI*2) * (floor((a) / (PI*2))))/* Fix angle in radians*/
#define dtr(x) ((x) * (PI / 180.0))                       /* Degree->Radian */
#define rtd(x) ((x) / (PI / 180.0))                       /* Radian->Degree */
#define NeedToCalculate(qty, time)  ((abs((qty##Last) - (time)) > qty##Interval) ? \
                                     (qty##Last = (time), TRUE) : FALSE)
#define InvalidateCalculation(qty)      qty##Last = -1e10
#define ValidateCalculation(qty, time) qty##Last = (time)
#define CalculationInterval(x)  ((x) / (24.0 * 60 * 60))
#define RepositionTime	10		// Reposition image every RepositionTime minutes
#define EOS     '\0'
//  Determine number of elements in an array
#define ELEMENTS(array) (sizeof(array)/sizeof((array)[0]))
//
/* Projection Modes */
#define SANSONS 1
#define STEREOGR 2
#define GNOMONIC 3
#define ORTHOGR 4
#define RECTANGULAR 5
#define OTHERPROJ 6
//
#define NO_CLIP 0
#define WEST_CLIP 1
#define EAST_CLIP 2
#define NORTH_CLIP 4
#define SOUTH_CLIP 8
#define RADIUS_CLIP 16
#define NE_CORNER 6
#define SE_CORNER 10
#define NW_CORNER 5
#define SW_CORNER 9
/* crude magnitudes of planets (x100) for output filtering by brightness */
#define MAGSOL 0
#define MAGMER 150
#define MAGVEN 0
#define MAGMAR 100
#define MAGJUP 0
#define MAGSAT 100
#define MAGURA 590
#define MAGNEP 800
#define MAGPLU 1400
//
#define TEXT_TYPE_CONSTEL 1
#define TEXT_TYPE_BRIGHT_ST 2
#define TEXT_TYPE_BAYERF_CODE 3
#define TEXT_TYPE_CELEST_EQ 4
#define TEXT_TYPE_ECLIPTIC 5

// Yale star data (based on Yale.bin)

// format:  Note: BYTE-PACKED!!! So don't use this struct!!!
struct starMapData
{
  unsigned short lon;
  short lat;
  unsigned short mag;
  unsigned char spectral;
};
//
// mag format: 0x f          fff
//                bitmap     mag
// bitmap:
// mag & 0x2000? --> smex (bayer/flams) is next 3 bytes, so read it now
// mag & 0x1000? --> impra, impdec are the next 4 bytes - read now
// mag & 0x8000? --> star name exists, read it.
//
// smex[0]<32? --> add 96 to smex[0], and use the greek alphabet
//                 otherwise use conventional alphabet



struct planet {                 // Planet information entry
	double hlong;				// Heliocentric longitude
	double hlat;				// Heliocentric latitude
	double hrv;					// Heliocentric radius vector
    double ra;                  // Current right ascension
    double dec;                 // Current declination
    double dist;				// Distance from the Earth
    double mag;					// Approximate magnitude
    double lha;					// Local hour angle
    double alt;					// Altitude above (-below) horizon
    double az;					// Azimuth from South: West positive, East negative
};
typedef struct tm_s
{
    int    tm_sec;   //seconds [0,61]
    int    tm_min;   //minutes [0,59]
    int    tm_hour;  //hour [0,23]
    int    tm_mday;  //day of month [1,31]
    int    tm_mon;   //month of year [0,11]
    int    tm_year;  //years since 1900
    int    tm_wday;  //day of week [0,6] (Sunday = 0)
    int    tm_yday;  //day of year [0,365]
    int    tm_isdst; //daylight savings flag
} tm_t;

typedef struct rect_t
{
    int left;
    int right;
    int top;
    int bottom;
} rect_s;

typedef struct {
  int width, height, x_offset, y_offset; /* Size and position,
                                            in integer device coords */
  /* The next several variables may be set by the driver, but the main routines
    may reset them (and the driver routines may then override that) */
  int proj_mode;                   
  /* Projection mode for this map */
  int invert;                       /* Invert (flip north south) */
  /* The following are set by the main routines */
  double racen, dlcen, scale;       /* R.A. and decl. of center, scale in degrees */
  double c_scale;                   /* One second of arc in display units */
} mapwindow;

// class definition
class Starmap
{
  public:
    Starmap();
    double jtime(tm_t *t);
    void set_col(uint16_t color);
    void paintSky(double limag, rect_s *br);

    //
    void highmoon(double jd, double *Lambda, double *Beta);
    void ecliptoeq(double jd, double Lambda, double Beta,
					   double *Ra, double *Dec);
    double gmst(double jd);
    void rgmst(double jd, double* ojd);
    void robliqeq(double jd, double *rdj);
    int initxform(mapwindow *win);
    void xform(double lat, double lon, int *xloc, int *yloc, int *inregion);
    int clipr_xform(double lat1, double lon1, double lat2, double lon2,
			int *xloc1, int *yloc1, int *xloc2, int *yloc2, int great_circle,
			double *plat1, double *plon1, double *plat2, double *plon2);
    void drawcurveline(double  lat1, double lon1, double lat2, double lon2,
			int xloc1, int yloc1, int xloc2, int yloc2,
			int line_style, int great_circle, int clevel);
    void planets(double jd, int which, struct planet *planet_info);
    double phase(double  pdate, double  *pphase, double  *mage, double  *dist,	double  *angdia, double  *sudist,	double  *suangdia);
    void jhms(double j, int *h, int *m, int *s);
    void jyear(double td, long *yy, int *mm, int *dd);
    void plot_char(char c, int x, int y, int color);
    void DrawMoon(int x, int y, int age);

    virtual void draw_line(int x0, int y0, int x1, int y1, uint16_t color);
    virtual void plot_pixel(uint16_t color, int x, int y);
    virtual void text_out(int x, int y, char* lab, unsigned char len, char type);
    virtual int storage_read(uint32_t addr, char* data, uint16_t len);


    char log2ram_buf[255];
    char mydeb;
    double jdtime;
    double siteLat;
    double siteLon;
    char dateSep[2], timeSep[2], amPM[2][5];
    int dateFormat, timeFormat;
    char display_ram[SCREEN_H][SCREEN_W];
    int col_coord_grid;
    int col_ecliptic;
    int col_constel;
    int col_stardim;
    int col_starbright;
    int col_startext;
    int col_moon_bright;
    int col_moon_dim;
    int col_moon_dark;
    int col_moon_phtext;

    int col_constel_text;
    int col_bright_st_text;
    int col_bayerf_text;
    int col_celest_eq_text;
    int col_ecliptic_text;

    // flags
    char do_constellation_text;
    int skyShowBflam;


  private:
  double ucttoj(long year, int mon, int mday, int hour, int min, int sec);
  void init_gt(mapwindow *win);
  void do_gt(double lat, double lon, double *xloc, double *yloc, double *r_theta);
  void inv_gt(double x, double y, double *latp, double *lonp);
  void plotLine(double fdec, double fra, double tdec, double tra);
  void MoveTo(int x, int y);
  void LineTo(int x, int y);
  void TextOut(int x, int y, char* lab, unsigned char len, char type);
  int OnScreen(int x, int y);
  void DrawStar(int x, int y, int level);
  double obliqeq(double jd);
  void definePrecession(double targetEpoch);
  void precessObject(double ira, double idec, double *ora, double *odec);
  double evalPoly(double a0, double a1, double a2, double a3, double t);
  double aint(double z);
  double range(double val);
  double kepler(double e, double M);
  double truean(double e, double E);
  double longi(double w2, double i, double u);
  double lati(double u, double i);
  void quadrat(double a, double b, double c, double *x_1, double *x_2, int *n);
  void gcmidpoint(double lat1, double lon1, double lat2, double lon2, double *pmlat, double *pmlon);
  void circ_intersect(double x_1, double y_1, double x_2, double y_2, double r, double *x1, double *y1, int *int_1, double *x2, double *y2, int *int_2);
  void speak(int which, double ra, double dec, double dis, int mag, struct planet *planet_info);
  void trans(int which, double r, double b, double ll, double Stheta, double Sr, double epli, int mag, struct planet *planet_info);
  void do_mercury(double T0, struct planet *planet_info);
  void do_venus(double T0, struct planet *planet_info);
  void do_mars(double T0, struct planet *planet_info);
  void do_jupiter(double T0, struct planet *planet_info);
  void do_saturn(double T0, struct planet *planet_info);
  void do_uranus(double T0, struct planet *planet_info);
  void do_neptune(double T0, struct planet *planet_info);
  int pluto(double jd, double *l, double *b, double *r);
  void calcPlanet(double jd);

  
  char mDummy;
  double preZeta;
  double preZ;
  double preTheta;
  double xf_north, xf_south, xf_bottom;
  int xf_xcen, xf_ycen, xf_ybot;
  int xf_w_left, xf_w_right, xf_w_top, xf_w_bot;
  double xf_c_scale;
  int xfs_proj_mode;
  double xfs_ra_cen, xfs_dl_cen, sin_dlcen, cos_dlcen, chart_scale;
  double xfs_scale;                /* Formerly yscale */
  double xfs_vinv, xfs_hinv;
  int xfs_wide_warn;
  double gt_sin_dlcen, gt_cos_dlcen, gt_chart_scale;
  double gt_scale;
  double gt_wx1, gt_wy1, gt_wx2, gt_wy2;
  double gt_ex1, gt_ey1, gt_ex2, gt_ey2;
  double gt_ny0, gt_na, gt_nb;
  double gt_sy0, gt_sa, gt_sb;
  double gt_r;
  int gt_use_boundaries;
  int clip_at1, clip_at2;
  double pie, radn;
  double plan_l, plan_a, plan_e, plan_i, plan_w1, plan_w2;
  double plan_M, plan_M0, plan_M1, plan_M2, plan_M4, plan_M5, plan_M6, plan_M7, plan_M8;
  double plan_RA, plan_DEC;
  double plan_ECC, plan_nu, plan_r, plan_u, plan_ll, plan_b, plan_lonpert, plan_radpert;
  double plan_esun, plan_Lsun, plan_Cen, plan_Snu;
  double plan_N, plan_D, plan_thapp, plan_omeg;
  double plan_nu2, plan_P, plan_Q, plan_S, plan_V, plan_W;
  double plan_ze, plan_l1pert, plan_epert, plan_w1pert, plan_apert;
  double plan_psi, plan_H, plan_G, plan_eta, plan_th;
  double plan_epli, plan_Stheta, plan_Sr;
  double skyLimitMag;
  int skyShowName;
  double skyNameMag;
  double skyBflamMag;
  int skyShowDeep;
  double skyDeepMag;
  int skyShowConstellations;
  int skyShowConbounds;
  int skyShowConnames;
  int skyAlignConnames;
  int skyShowCoords;
  int skyShowPlanets;
  int skyShowTiming;
  int precessionCalculation;
  struct planet planet_info[11];
  mapwindow skywin;
  char firstTick;
  char nextTick;
  char bailedOut;
  long calculationNumber;
  double skyLham;
  int Flip;
  long tickerLast;
  int pirefc;
  int starCatalogue;
  int starQuality;
  int showStarColours;
  int draw_x_store, draw_y_store;
  uint16_t draw_col_store;
  int plutoPrecise;




};

#endif // STARMAP_H
