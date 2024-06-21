# starmap

The code in this repository is a star/planets engine in an Arduino library.

The ready-made Arduino library zip file can be found by clicking on **Releases**.

If you want to build the library from source, you can clone the code into a folder onto your Linux machine, then go into the starmap folder, and type:

```
chmod 755 make_libraries.py
./make_libraries.py
```

The output will be two zip files, one folder up from the starmap folder. The zip file titled **starmap_library-1.0.0.zip** will be the bundled library that you can add into the Arduino development environment (using **Sketch->Include Library->Add .ZIP file**).

The other zip file that will be generated, called **starmap_lib.zip** is useful for CMake projects for non-Arduino platforms, such as Pi Pico C/C++ SDK. It can be ignored if you're using Arduino.

If you wish to test the code on Linux, prior to uploading to any microcontroller, then you can do that by typing:

```
cd starmap/src
mkdir -p build
cd build
cmake ..
make
```

By doing that, a **starmap** Linux executable will be built in the build folder. You can run it by typing:

```
./starmap
```

It will generate a starmap PNG file (called out.png) which can be inspected.


## Acknowlegdement

Note: This code uses content from https://www.fourmilab.ch/homeplanet/ (public domain code).

See the screenshot here:

<img width="100%" align="left" src="screenshot-public-domain-code.jpg">


