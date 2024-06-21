// StarmapTest.ino
// rev 1 - May 2024 - shabaz
// This code demonstrates the Starmap library

// ******** includes ********
#include <Starmap.h>
#include <SPI.h>
#include <Wire.h>

#include <HardwareSerial.h>


// ******** defines ********
// defines
#define FOREVER 1
#define DELAY_MS delay
// screen dimensions
#define TFT_W 240
#define TFT_H 240
// default co-ordinates, lat: deg N, lon: deg W
#define DEFAULT_LAT 47.0
#define DEFAULT_LON 122.0


// ***** class based on Starmap ******
class SM : public Starmap {
  // you need to implement either plot_pixel or draw_line
  void plot_pixel(uint16_t color, int x, int y);
  // void draw_line(int x0, int y0, int x1, int y1, uint16_t color);
  // optionally you can also implement text_out
  // void text_out(int x, int y, char* lab, unsigned char len);
};

// ********global variables ********
// create starmap object
SM starmap;
double mag;
rect_s screen_rect;
tm_t tm;


// ******* plot_pixel function ******
void SM::plot_pixel(uint16_t color, int x, int y) {
  if ((x>=TFT_W) || (x<0) || (y>=TFT_H) || (y<0)) return;
  // handle your TFT here
}

// ******** setup() function ********
void setup() {
  // put your setup code here, to run once:
  Serial.begin(115200);
  // starmap setup
  mag = 5; // render magnitude level
  screen_rect.left = 0;
  screen_rect.right = TFT_W;
  screen_rect.top = 0;
  screen_rect.bottom = TFT_H;
  starmap.siteLat = DEFAULT_LAT;
  starmap.siteLon = DEFAULT_LON;
}

// ******** loop() function ********
void loop() {
  // put your main code here, to run repeatedly:
  tm.tm_sec = 0; // seconds 0-59
  tm.tm_min = 18; // minutes 0-59
  tm.tm_hour = 23; // hour 0-23
  tm.tm_mday = 16; // date 1-31
  tm.tm_mon = 7; // month 0-11
  tm.tm_year = 2004 - 1900; // year since 1900. Example: 100 means year 2000
  starmap.jdtime = starmap.jtime(&tm); // store the time in Julian days

  starmap.paintSky(mag, &screen_rect); // paint the sky!
  Serial.print("done painting sky\n");
  
  
  while(FOREVER) {
    DELAY_MS(1000);
  }

}

