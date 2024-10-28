// StarmapProjectESP32.ino
// rev 1 - June 2024 - shabaz
// This code implements a star chart in an ESP32

// ******** includes ********
#include <Starmap.h>
#include <Flash_SST25VF.h>
//#include <Teseo.h>
#include <Arduino_GFX_Library.h>

// ******** defines ********
#define NO_GPS
#define FOREVER 1
#define DELAY_MS delay
// time offset, example: 1 hour ahead of UTC (e.g. British Summer Time) is 1
#define DISPLAYED_TIME_OFFSET 1
// GPS (uses Serial1, requires an Arduino which supports that)
// Note: set GNSS_BAUDRATE to 9600 for a new GPS module.
// I'm using 115200 because I've configured my GPS module to that baudrate
#define GNSS_BAUDRATE 115200
// screen dimensions
#define TFT_W 480
#define TFT_H 480
// default co-ordinates, lat: deg N, lon: deg W
#define DEFAULT_LAT 47.0
#define DEFAULT_LON 122.0
// hardware SPI CS must be set to the same as TFT_CS
#define FLASH_CS A0
// Flash memory command to read data
#define FLASH_READ 0x03
#define DATA_BUF_SIZE 256
// ***** flags *********
#define FLAG_USA 0
#define FLAG_EU 1
#define FLAG_CHINA 2
#define FLAG_RUSSIA 3
// flag colors in RGB565 format
#define FLAG_BORDER 0x0000
// USA
#define US_RED 0xf800
#define US_WHITE 0xffff
#define US_BLUE 0x001f
// Russia
#define RU_RED 0xf800
#define RU_WHITE 0xffff
#define RU_BLUE 0x001f
// China
#define CN_RED 0xf800
#define CN_YELLOW 0xffe0
// EU
#define EU_BLUE 0x001f
#define EU_YELLOW 0xffe0
// flag icons
const uint16_t flag_usa[7][7] = {
    {US_RED, US_WHITE, US_RED, US_BLUE, US_BLUE, US_BLUE, US_BLUE},
    {US_RED, US_WHITE, US_RED, US_BLUE, US_WHITE, US_BLUE, US_WHITE},
    {US_RED, US_WHITE, US_RED, US_BLUE, US_BLUE, US_BLUE, US_BLUE},
    {US_RED, US_WHITE, US_RED, US_BLUE, US_WHITE, US_BLUE, US_WHITE},
    {US_RED, US_WHITE, US_RED, US_WHITE, US_RED, US_WHITE, US_RED},
    {US_RED, US_WHITE, US_RED, US_WHITE, US_RED, US_WHITE, US_RED},
    {US_RED, US_WHITE, US_RED, US_WHITE, US_RED, US_WHITE, US_RED}
};
const uint16_t flag_russia[7][6] = {
    {RU_RED, RU_RED, RU_BLUE, RU_BLUE, RU_WHITE, RU_WHITE},
    {RU_RED, RU_RED, RU_BLUE, RU_BLUE, RU_WHITE, RU_WHITE},
    {RU_RED, RU_RED, RU_BLUE, RU_BLUE, RU_WHITE, RU_WHITE},
    {RU_RED, RU_RED, RU_BLUE, RU_BLUE, RU_WHITE, RU_WHITE},
    {RU_RED, RU_RED, RU_BLUE, RU_BLUE, RU_WHITE, RU_WHITE},
    {RU_RED, RU_RED, RU_BLUE, RU_BLUE, RU_WHITE, RU_WHITE},
    {RU_RED, RU_RED, RU_BLUE, RU_BLUE, RU_WHITE, RU_WHITE}
};
const uint16_t flag_china[7][7] = {
    {CN_RED, CN_RED, CN_RED, CN_YELLOW, CN_RED, CN_YELLOW, CN_YELLOW},
    {CN_RED, CN_RED, CN_RED, CN_RED, CN_RED, CN_RED, CN_YELLOW},
    {CN_RED, CN_RED, CN_RED, CN_RED, CN_YELLOW, CN_RED, CN_RED},
    {CN_RED, CN_RED, CN_RED, CN_RED, CN_RED, CN_RED, CN_YELLOW},
    {CN_RED, CN_RED, CN_RED, CN_RED, CN_RED, CN_RED, CN_RED},
    {CN_RED, CN_RED, CN_RED, CN_RED, CN_RED, CN_RED, CN_RED},
    {CN_RED, CN_RED, CN_RED, CN_RED, CN_RED, CN_RED, CN_RED}
};
const uint16_t flag_eu[7][7] = {
    {EU_BLUE, EU_BLUE, EU_BLUE, EU_YELLOW, EU_BLUE, EU_BLUE, EU_BLUE},
    {EU_BLUE, EU_YELLOW, EU_BLUE, EU_BLUE, EU_BLUE, EU_YELLOW, EU_BLUE},
    {EU_BLUE, EU_BLUE, EU_BLUE, EU_BLUE, EU_BLUE, EU_BLUE, EU_BLUE},
    {EU_YELLOW, EU_BLUE, EU_BLUE, EU_BLUE, EU_BLUE, EU_BLUE, EU_YELLOW},
    {EU_BLUE, EU_BLUE, EU_BLUE, EU_BLUE, EU_BLUE, EU_BLUE, EU_BLUE},
    {EU_BLUE, EU_YELLOW, EU_BLUE, EU_BLUE, EU_BLUE, EU_YELLOW, EU_BLUE},
    {EU_BLUE, EU_BLUE, EU_BLUE, EU_YELLOW, EU_BLUE, EU_BLUE, EU_BLUE}
};
//10x12 font for N,E,S,W characters only
const uint16_t font10_12[4][12] = {
    {0x0fff, 0x0fff, 0x0700, 0x03c0, 0x01e0, 0x0078, 0x003c, 0x000e, 0x0fff, 0x0fff},
    {0x0fff, 0x0fff, 0x0c63, 0x0c63, 0x0c63, 0x0c63, 0x0c63, 0x0c63, 0x0c03, 0x0c03},
    {0x038c, 0x07ce, 0x0ee7, 0x0c63, 0x0c63, 0x0c63, 0x0c63, 0x0e77, 0x073e, 0x031c},
    {0x0f80, 0x0ff8, 0x00ff, 0x0007, 0x007e, 0x007e, 0x0007, 0x00ff, 0x0ff8, 0x0f80}
};
// TFT configuration commands array
static const uint8_t st7701_2_8_aliexpress_operations[] = {
    BEGIN_WRITE, WRITE_COMMAND_8, 0xFF, WRITE_BYTES, 5, 0x77, 0x01, 0x00, 0x00, 0x13, WRITE_C8_D8, 0xEF, 0x08, WRITE_COMMAND_8, 0xFF,
    WRITE_BYTES, 5, 0x77, 0x01, 0x00, 0x00, 0x10, WRITE_C8_D16, 0xC0, 0x3B, 0x00, WRITE_C8_D16, 0xC1, 0x10, 0x0C, WRITE_C8_D16, 0xC2, 0x07, 0x0A, 
    WRITE_C8_D8, 0xC7, 0x04, WRITE_C8_D8, 0xCC, 0x10, WRITE_C8_D8, 0xCD, 0x08, WRITE_COMMAND_8, 0xB0, WRITE_BYTES, 16, 0x05, 0x12, 0x98, 0x0E,
    0x0F, 0x07, 0x07, 0x09, 0x09, 0x23, 0x05, 0x52, 0x0F, 0x67, 0x2C, 0x11, WRITE_COMMAND_8, 0xB1, WRITE_BYTES, 16, 0x0B, 0x11, 0x97, 0x0C,
    0x12, 0x06, 0x06, 0x08, 0x08, 0x22, 0x03, 0x51, 0x11, 0x66, 0x2B, 0x0F, WRITE_COMMAND_8, 0xFF, WRITE_BYTES, 5, 0x77, 0x01, 0x00, 0x00, 0x11,
    WRITE_C8_D8, 0xB0, 0x5D, WRITE_C8_D8, 0xB1, 0x3e, WRITE_C8_D8, 0xB2, 0x81, WRITE_C8_D8, 0xB3, 0x80, WRITE_C8_D8, 0xB5, 0x4E,
    WRITE_C8_D8, 0xB7, 0x85, WRITE_C8_D8, 0xB8, 0x20, WRITE_C8_D8, 0xC1, 0x78, WRITE_C8_D8, 0xC2, 0x78, WRITE_C8_D8, 0xD0, 0x88,
    WRITE_COMMAND_8, 0xE0, WRITE_BYTES, 3, 0x00, 0x00, 0x02, WRITE_COMMAND_8, 0xE1, WRITE_BYTES, 11, 0x06, 0x30, 0x08, 0x30, 0x05, 0x30, 0x07, 0x30,
    0x00, 0x33, 0x33, WRITE_COMMAND_8, 0xE2, WRITE_BYTES, 12, 0x11, 0x11, 0x33, 0x33, 0xF4, 0x00, 0x00, 0x00, 0xF4, 0x00, 0x00, 0x00,
    WRITE_COMMAND_8, 0xE3, WRITE_BYTES, 4, 0x00, 0x00, 0x11, 0x11, WRITE_C8_D16, 0xE4, 0x44, 0x44, WRITE_COMMAND_8, 0xE5, WRITE_BYTES, 16,
    0x0D, 0xF5, 0x30, 0xF0, 0x0F, 0xF7, 0x30, 0xF0, 0x09, 0xF1, 0x30, 0xF0, 0x0B, 0xF3, 0x30, 0xF0, WRITE_COMMAND_8, 0xE6, WRITE_BYTES, 4, 
    0x00, 0x00, 0x11, 0x11, WRITE_C8_D16, 0xE7, 0x44, 0x44, WRITE_COMMAND_8, 0xE8, WRITE_BYTES, 16, 0x0C, 0xF4, 0x30, 0xF0, 0x0E, 0xF6, 0x30, 0xF0,
    0x08, 0xF0, 0x30, 0xF0, 0x0A, 0xF2, 0x30, 0xF0, WRITE_C8_D16, 0xE9, 0x36, 0x01, WRITE_COMMAND_8, 0xEB, WRITE_BYTES, 7, 0x00, 0x01, 0xE4, 0xE4,
    0x44, 0x88, 0x40, WRITE_COMMAND_8, 0xED, WRITE_BYTES, 16, 0xFF, 0x10, 0xAF, 0x76, 0x54, 0x2B, 0xCF, 0xFF, 0xFF, 0xFC, 0xB2, 0x45,
    0x67, 0xFA, 0x01, 0xFF, WRITE_COMMAND_8, 0xEF, WRITE_BYTES, 6, 0x08, 0x08, 0x08, 0x45, 0x3F, 0x54, WRITE_COMMAND_8, 0xFF, WRITE_BYTES, 5, 
    0x77, 0x01, 0x00, 0x00, 0x00, WRITE_COMMAND_8, 0x11, END_WRITE, DELAY, 120, BEGIN_WRITE, WRITE_C8_D8, 0x3A, 0x66, WRITE_C8_D8, 0x36, 0x00,
    WRITE_C8_D8, 0x35, 0x00, WRITE_COMMAND_8, 0x29, END_WRITE};

#define NORTH_SYMBOL 0
#define EAST_SYMBOL 1
#define SOUTH_SYMBOL 2
#define WEST_SYMBOL 3
#define NESW_COLOR 0xf800
// misc
#ifndef PI
#define PI 3.14159265358979323846
#endif

// ***** class based on Starmap ******
class SM : public Starmap {
  // you need to implement plot_pixel and/or draw_line
  // if you want text output, then implement either plot_pixel or text_out
  void plot_pixel(uint16_t color, int x, int y);
  void draw_line(int x0, int y0, int x1, int y1, uint16_t color);
  // optionally you can also implement text_out
  // void text_out(int x, int y, char* lab, unsigned char len, char type);
  // if star plotting is required, then storage_read is implemented
  int storage_read(uint32_t addr, char* data, uint16_t len);
};

// ***** class based on Flash_SST25VF ********
class Flash : public Flash_SST25VF
{
  public:
    Flash(uint8_t pin_CS) : Flash_SST25VF(pin_CS) {};
    // override/release functions to play nice with
    // other devices on the SPI bus that may be using
    // the default CS pin
    void override_default_cs();  // override the default CS pin if required
    void release_default_cs();  // release the default CS pin if required
};

// ******** global variables ********
SM starmap;
Flash flash(FLASH_CS); // create an instance of the Flash class
#ifndef NO_GPS
Teseo gnss;
#endif

Arduino_XCA9554SWSPI *expander = new Arduino_XCA9554SWSPI(PCA_TFT_RESET, PCA_TFT_CS, PCA_TFT_SCK, PCA_TFT_MOSI, &Wire, 0x3F);
Arduino_ESP32RGBPanel *rgbpanel = new Arduino_ESP32RGBPanel(TFT_DE, TFT_VSYNC, TFT_HSYNC, TFT_PCLK, TFT_R1, TFT_R2, TFT_R3, TFT_R4, TFT_R5,
    TFT_G0, TFT_G1, TFT_G2, TFT_G3, TFT_G4, TFT_G5, TFT_B1, TFT_B2, TFT_B3, TFT_B4, TFT_B5, 1 , 50, 2, 44, 1, 16, 2, 18);
Arduino_RGB_Display *gfx = new Arduino_RGB_Display(480, 480, rgbpanel, 0 /* rotation */, true /* auto_flush */,
    expander, GFX_NOT_DEFINED /* RST */, st7701_2_8_aliexpress_operations, sizeof(st7701_2_8_aliexpress_operations));


double mag;
rect_s screen_rect;
tm_t tm;
char first_fix_done = 0;
char flash_present = 0; // set to 1 if the flash is present and valid
char flash_data[DATA_BUF_SIZE]; // buffer used for storing data read from Flash
uint32_t flash_data_addr; // starting address of the contents of flash_data
uint32_t starmap_update_period = 5 * 60000; // multiply by 60000 to convert minutes to milliseconds
uint8_t sat_id[24]; // buffer to store the satellite ID


// ******** function prototypes ************
void draw_flag(int x, int y, int flagtype);
void ccw_azimuth_elevation_to_xy(double az, double el, int *x, int *y);
void get_first_fix(void);
#ifndef NO_GPS
int copy_gnss_time_to_starmap(void);
int copy_gnss_loc_to_starmap(void);
#endif
void invalidate_displayed_sat_list(void);
int check_and_add_sat(uint8_t id);
void plot_char_10_12(char c, int x, int y, int color); // only supports N, E, S, W characters
void disp_lat_lon(double lat, double lon, int x, int y, int col);

// ****** flash functions ***********
void Flash::override_default_cs()
{
    // comment the line below if the default SPI CS pin is unused by any device
    //pinMode(DEFAULT_CS, INPUT_PULLUP);
}

void Flash::release_default_cs()
{
    // comment the line below if the default SPI CS pin is unused by any device
    //pinMode(DEFAULT_CS, OUTPUT);
}

// ******* plot_pixel function ******
void SM::draw_line(int x0, int y0, int x1, int y1, uint16_t color) {
  // sanity check
  if (x0<0) x0=0;
  if (x1<0) x1=0;
  if (y0<0) y0=0;
  if (y1<0) y1=0;
  if (x0>=TFT_W) x0=TFT_W-1;
  if (x1>=TFT_W) x1=TFT_W-1;
  if (y0>=TFT_H) y0=TFT_H-1;
  if (y1>=TFT_H) y1=TFT_H-1;
  // handle your TFT here
  gfx->drawLine(x0, y0, x1, y1, color);
}

void SM::plot_pixel(uint16_t color, int x, int y) {
    // sanity check
    if (x<0) x=0;
    if (y<0) y=0;
    if (x>=TFT_W) x=TFT_W-1;
    if (y>=TFT_H) y=TFT_H-1;
    // handle your TFT here
    gfx->drawPixel(x, y, color);
}

// ******* storage_read function *******
int SM::storage_read(uint32_t addr, char* data, uint16_t len){
    // read from Flash memory
    flash.flash_read(addr, data, len);
    return 0;
}

// ******** setup() function ********
void setup() {
  // put your setup code here, to run once:
  Serial.begin(115200);
  // TFT setup
#ifdef GFX_EXTRA_PRE_INIT
  GFX_EXTRA_PRE_INIT();
#endif
  Wire.setClock(1000000); // speed up I2C 
  if (!gfx->begin()) {
    Serial.println("gfx->begin() failed!");
  }
  gfx->fillScreen(BLACK);
  // GPS setup
  #ifndef NO_GPS
  gnss.init(GNSS_BAUDRATE);
  #endif
  // starmap setup
  mag = 5; // render magnitude level
  screen_rect.left = 0;
  screen_rect.right = TFT_W;
  screen_rect.top = 0;
  screen_rect.bottom = TFT_H;
  starmap.siteLat = DEFAULT_LAT;
  starmap.siteLon = DEFAULT_LON;
  // colors
  starmap.col_coord_grid = 0x4a49; // dark gray
  starmap.col_ecliptic = 0xab91; // dark pink-red
  starmap.col_constel = 0x326b; // dark aqua
  starmap.col_stardim = 0xa520; // dim yellow
  starmap.col_starbright = 0xffe0; // bright yellow
  starmap.col_startext = 0x001f; // pure blue
  starmap.col_moon_bright = 0xffff; // white
  starmap.col_moon_dim = 0xe71c; // light gray
  starmap.col_moon_dark = 0xce79; // med gray
  starmap.col_moon_phtext = 0x0000; // black
  starmap.col_bright_st_text = 0x001f; // pure blue
  starmap.col_ecliptic_text = 0xab91; // dark pink-red
  starmap.col_celest_eq_text = 0x4a49; // dark gray
  starmap.col_constel_text = 0x0030; // dark blue

  // check if flash is available and valid
  if (flash.flash_read_id() == 0x8e) {
    // Flash chip is installed, but check if it contains the magic string:
    flash_data[0] = 0;
    flash.flash_read(0x00000, flash_data, 12);
    if (strcmp(flash_data, "starmap v01")==0) {
      flash_present=1; // magic string is present.
    }
  }

  DELAY_MS(3000);
  Serial.print("setup() complete\n\r");
  if(flash_present) {
    Serial.print("flash present and valid\n\r");
  }
}

// ******** loop() function ********
void loop() {
  double lat, lon;
  int valid=0;
  int i;
  int x, y;
  uint32_t ts; // timestamp in milliseconds
  int hr, min, sec; // used to store the displayed time
  int old_min;
  char disp_lat_lon_complete;
  char text_string[24];

  if (!first_fix_done) {
    // wait for GNSS to start providing data
    gfx->fillScreen(BLACK); // clear the TFT screen
    gfx->setCursor(40, 240);
    gfx->setTextColor(WHITE);
    gfx->setTextSize(1);
    gfx->println("Acquiring Satellite Fix..");
    yield();
    get_first_fix(); // executes gnss.get_data() repeatedly until we have a valid fix
    first_fix_done = 1;
  } else {
#ifndef NO_GPS
    gnss.flush_buffer(); // clear out the buffer
    gnss.get_data(PRINT_ENABLE); // get GNSS data and print the raw sentence to Serial
#endif
  }

#ifdef NO_GPS
  tm.tm_sec = 0; // seconds 0-59
  tm.tm_min = 15; // minutes 0-59
  tm.tm_hour = 02; // hour 0-23
  tm.tm_mday = 11; // date 1-31
  tm.tm_mon = 7 - 1; // month 0-11
  tm.tm_year = 2024 - 1900; // year since 1900. Example: 100 means year 2000
  starmap.jdtime = starmap.jtime(&tm); // store the time in Julian days
#else
  // fill starmap object with the GNSS time and location
  copy_gnss_time_to_starmap();
  copy_gnss_loc_to_starmap();
#endif

  gfx->fillScreen(BLACK); // clear the TFT screen
  yield();
  Serial.print("executing paintSky..\n\r");
  starmap.paintSky(mag, &screen_rect); // paint the sky!
  yield();
  Serial.print("done executing paintSky.\n\r");

  // print lat and lon on tft
    //lat = starmap.siteLat;
    //lon = starmap.siteLon;
    //disp_lat_lon(lat, lon, 40, 80, WHITE);


  invalidate_displayed_sat_list(); // clear the list of displayed satellites
  old_min = -1;
  disp_lat_lon_complete = 0;
  // loop until the starmap needs to be redrawn
  ts = millis(); // get current Arduino timestamp
  while (millis() - ts < starmap_update_period) {

    //if (disp_lat_lon_complete==0) {
    //    if ((millis() - ts > 10000)) { // display the lat and lon for 10 seconds
    //        gfx->fillScreen(BLACK); // clear the TFT screen
    //        starmap.paintSky(mag, &screen_rect);
    //        yield();
    //        disp_lat_lon_complete = 1;
    //    }
    //}

    DELAY_MS(1100);
#ifdef NO_GPS
#else
    gnss.get_data(PRINT_DISABLE);
#endif
    // store the time so we can display it
#ifdef NO_GPS
#else
    hr = gnss.rmc.hour + DISPLAYED_TIME_OFFSET;
    if (hr > 23) hr -= 24;
    min = gnss.rmc.min;
    sec = gnss.rmc.sec;
    char sat_drawn = 0;
    for (i=0; i<gnss.satnum; i++) {
        //Serial.print("sat_id[0] = ");
        //Serial.print((uint8_t)sat_id[0]);
        //Serial.print("\n\r");
        if (check_and_add_sat((uint8_t)gnss.gsv[i].prn) < 1) {
            // satellite already in the list, no need to redraw it.
            continue;
        }
        ccw_azimuth_elevation_to_xy(gnss.gsv[i].azim, gnss.gsv[i].elev, &x, &y);
        switch(gnss.gsv[i].source) {
            case SOURCE_GPS:
                draw_flag(x, y, FLAG_USA);
                sat_drawn = 1;
                break;
            case SOURCE_GALILEO:
                draw_flag(x, y, FLAG_EU);
                sat_drawn = 1;
                break;
            case SOURCE_BEIDOU:
                draw_flag(x, y, FLAG_CHINA);
                sat_drawn = 1;
                break;
            case SOURCE_GLONASS:
                draw_flag(x, y, FLAG_RUSSIA);
                sat_drawn = 1;
                break;
            case SOURCE_UNK:
                // should draw some default flag
                sat_drawn = 1;
                Serial.print("warning! unknown satellite source.\n");
                break;
            default:
                Serial.print("error! unexpected satellite source.\n");
                break;
        }
        if (sat_drawn) { // don't want satellite icons to overlap the text
            // overlay N, E, S, W characters onto the screen
            plot_char_10_12(NORTH_SYMBOL, 115, 18, NESW_COLOR);
            plot_char_10_12(EAST_SYMBOL, 3, 125, NESW_COLOR);
            plot_char_10_12(SOUTH_SYMBOL, 115, 234, NESW_COLOR);
            plot_char_10_12(WEST_SYMBOL, 225, 125, NESW_COLOR);
            // display lat and lon
            lat = starmap.siteLat;
            lon = starmap.siteLon;
            disp_lat_lon(lat, lon, 52, 190, WHITE);
        }
    }
#endif
    // display the time if it has changed
    if (min != old_min) {
        old_min = min;
        sprintf(text_string, "%02d:%02d", hr, min);
        gfx->fillRect(91, 199, 60, 16, BLACK);
        gfx->setCursor(92, 200);
        gfx->setTextColor(WHITE);
        gfx->setTextSize(2);
        gfx->println(text_string);
        yield();
    }
  } // end while loop for starmap_update_period


  // starmap has been displayed for starmap_update_period
  // now loop back to the beginning, so that the starmap is completely redrawn
}

// ************ other functions *****************
// plot_char_10_12 only supports N, E, S, W characters
void plot_char_10_12(char c, int x, int y, int color) {
    int i, j;
    // check that the character is in the font
    // only 0-3 are valid (N, E, S, W)
    if (c > 3 || c < 0)
    {
        return;
    }

    for (i = 0; i < 10; i++)
    {
        for (j = 0; j < 12; j++)
        {
            if (font10_12[c][i] & (1 << j))
            {
                gfx->drawPixel(x + i, y - j, color);
            } else {
                gfx->drawPixel(x + i, y - j, BLACK);
            }
        }
    }
}

void draw_flag(int x, int y, int flagtype) {
    int i, j;
    for (i=0; i<7; i++) {
        for (j=0; j<7; j++) {
            switch (flagtype) {
                case FLAG_USA:
                    gfx->drawPixel(x+i-3, y-j+3, flag_usa[i][j]);
                    //starmap.plot_pixel(flag_usa[i][j], x+i-3, y-j+3);
                    break;
                case FLAG_EU:
                    gfx->drawPixel(x+i-3, y-j+3, flag_eu[i][j]);
                    //starmap.plot_pixel(flag_eu[i][j], x+i-3, y-j+3);
                    break;
                case FLAG_CHINA:
                    gfx->drawPixel(x+i-3, y-j+3, flag_china[i][j]);
                    //starmap.plot_pixel(flag_china[i][j], x+i-3, y-j+3);
                    break;
                case FLAG_RUSSIA:
                    if (j<6) {
                        gfx->drawPixel(x+i-3, y-j+3, flag_russia[i][j]);
                        //starmap.plot_pixel(flag_russia[i][j], x+i-3, y-j+3);
                    }
                    break;
            }
        }
    }
    // draw a  border outside each flag, using draw_line
    gfx->drawLine(x-4, y-4, x+4, y-4, FLAG_BORDER);
    gfx->drawLine(x-4, y+4, x-4, y-4, FLAG_BORDER);
    gfx->drawLine(x+4, y+4, x+4, y-4, FLAG_BORDER);
    //starmap.draw_line(x-4, y+4, x+4, y+4, 0x0000);
    //starmap.draw_line(x-4, y+4, x-4, y-4, 0x0000);
    //starmap.draw_line(x+4, y+4, x+4, y-4, 0x0000);
    if (flagtype != FLAG_RUSSIA) {
        gfx->drawLine(x-4, y-4, x+4, y-4, FLAG_BORDER);
        //starmap.draw_line(x-4, y-4, x+4, y-4, 0x0000);
    } else {
        gfx->drawLine(x-4, y-3, x+4, y-3, FLAG_BORDER);
        //starmap.draw_line(x-4, y-3, x+4, y-3, 0x0000);
    }
    yield();
}

// function to convert elevation and azimuth to x,y coordinates
// with north-up, but in a counter-clockwise fashion, since we are 
// plotting an overhead view of the sky
void ccw_azimuth_elevation_to_xy(double az, double el, int *x, int *y) {
  // convert to radians
  az = az * PI / 180.0;
  el = el * PI / 180.0;
  // calculate x and y
  *x = (int)(TFT_W/2 + (TFT_W/2-1) * cos(el) * sin(az));
  *y = (int)(TFT_H/2 - (TFT_H/2-1) * cos(el) * cos(az));
  // we want anticlockwise rotation
  *x = TFT_W - *x;
}

#ifdef NO_GPS
#else
// copy GNSS data to Starmap object
int copy_gnss_time_to_starmap(void) {
  int valid = 1;
  if ((gnss.rmc.sec<0) || (gnss.rmc.sec>59)) valid = 0;
  if ((gnss.rmc.min<0) || (gnss.rmc.min>59)) valid = 0;
  if ((gnss.rmc.hour<0) || (gnss.rmc.hour>23)) valid = 0;
  if ((gnss.rmc.date<1) || (gnss.rmc.date>31)) valid = 0;
  if ((gnss.rmc.month<1) || (gnss.rmc.month>12)) valid = 0;
  if ((gnss.rmc.year<2000) || (gnss.rmc.year>2199)) valid = 0;

  if (!valid) return -1;

  // populate the starmap object with the time
  tm.tm_sec = gnss.rmc.sec; // seconds 0-59
  tm.tm_min = gnss.rmc.min; // minutes 0-59
  tm.tm_hour = gnss.rmc.hour; // hour 0-23
  tm.tm_mday = gnss.rmc.date; // date 1-31
  tm.tm_mon = gnss.rmc.month - 1; // month 0-11
  tm.tm_year = gnss.rmc.year - 1900; // year since 1900. Example: 100 means year 2000
  starmap.jdtime = starmap.jtime(&tm); // store the time in Julian days
  return 1;
}

int copy_gnss_loc_to_starmap(void) {
  int valid = 1;
  if ((gnss.rmc.lat<-90) || (gnss.rmc.lat>90)) valid = 0;
  if ((gnss.rmc.lon<-180) || (gnss.rmc.lon>180)) valid = 0;

  if (!valid) {
    Serial.print("error! invalid GNSS location: lat=");
    Serial.print(gnss.rmc.lat);
    Serial.print(", lon=");
    Serial.print(gnss.rmc.lon);
    Serial.print("\n\r");
    return -1;
  } 

  // populate the starmap object with the location
  starmap.siteLat = gnss.rmc.lat;
  // GNSS longitude is in degrees East, but starmap object requires degrees West
    starmap.siteLon = 0.0 - gnss.rmc.lon;
  return 1;
}
#endif

#ifdef NO_GPS
void get_first_fix(void) {
  //
}
#else
void get_first_fix(void) {
  double lat, lon;
  int valid=0;
  // wait for GNSS to start providing data
  // set to invalid values
  gnss.rmc.lat = 9999;
  gnss.rmc.lon = 9999;
  gnss.rmc.hour = 9999;
  while (!valid) {
    DELAY_MS(1100);
    gnss.flush_buffer(); // clear out the buffer
    gnss.get_data(PRINT_ENABLE); // get GNSS data and print the raw sentence to Serial
    if ((gnss.rmc.lat < 999) && (gnss.rmc.lon < 999) && (gnss.rmc.hour < 999)) {
      // ok we seem to have received at least lat/lon and time
      valid = 1;
    }
  }
}
#endif

void invalidate_displayed_sat_list(void) {
  int i;
  for (i=0; i<24; i++) {
    sat_id[i] = 255;
  }
}

#ifdef NO_GPS
int check_and_add_sat(uint8_t id) {
  return(0);
}
#else
int check_and_add_sat(uint8_t id) {
  int i;
  // check if the id is valid
  if (id == SOURCE_INVALID) {
    Serial.print("error! invalid satellite ID.\n\r");
    return(0);
  }
  // check if the satellite is already in the list
  for (i=0; i<24; i++) {
    if (sat_id[i] == 255) {
      // end of list
      sat_id[i] = id;
      return 1;
    }
    if (sat_id[i] == id) {
      // already in the list
      return 0;
    }
  }
  // list is full
  Serial.print("error! satellite list is full.\n\r");
  return(0);
}
#endif

void disp_lat_lon(double lat, double lon, int x, int y, int col) {
    char text_string[32];
    char neg_lat=0;
    char neg_lon=0;
    int width;
    if (lat < 0) {
        neg_lat = 1;
        lat = 0 - lat;
    }
    if (lon < 0) {
        neg_lon = 1;
        lon = 0 - lon;
    }
    // build text string; there must be an easier way to do this!
    if (neg_lat) {
        if (neg_lon) {
            sprintf(text_string, "LAT:-%02d.%03dN LON:-%02d.%03dW", (int)lat, (int)((lat - (int)lat) * 1000), (int)lon, (int)((lon - (int)lon) * 1000));
        } else {
            sprintf(text_string, "LAT:-%02d.%03dN LON:%02d.%03dW", (int)lat, (int)((lat - (int)lat) * 1000), (int)lon, (int)((lon - (int)lon) * 1000));
        }
    } else {
        if (neg_lon) {
            sprintf(text_string, "LAT:%02d.%03dN LON:-%02d.%03dW", (int)lat, (int)((lat - (int)lat) * 1000), (int)lon, (int)((lon - (int)lon) * 1000));
        } else {
            sprintf(text_string, "LAT:%02d.%03dN LON:%02d.%03dW", (int)lat, (int)((lat - (int)lat) * 1000), (int)lon, (int)((lon - (int)lon) * 1000));
        }
    }

    width = 158;
    if (neg_lat) width += 8;
    if (neg_lon) width += 8;
    gfx->fillRect(x-1, y-1, width, 10, BLACK);
    gfx->setCursor(x, y);
    gfx->setTextColor(col);
    gfx->setTextSize(1);
    gfx->println(text_string);
}
