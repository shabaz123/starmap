#!/usr/bin/env python3
import sys
import os

# this function reads the library.properties file
def get_library_properties():
    properties = {}
    with open('library.properties') as f:
        for line in f:
            if '=' in line:
                key, value = line.split('=', 1)
                properties[key] = value.strip()
    return properties

properties = get_library_properties()
name = properties['name'].split(' ')[0].lower()
name += '_library-' + properties['version']

if len(name) < 3:
  print('Error, library name is not present or is too short in library.properties')
  sys.exit(1)

print(f'This script will build an Arduino library file called ../{name}.zip')
print('If this is incorrect, you should abort and check your library.properties file')
response = input('Do you wish to proceed? [y/n]: ')
if response != 'y':
  print('Aborted.')
  sys.exit(1)

print(f'creating folder ../{name}')
if os.path.exists('../' + name):
  os.system('rm -rf ../' + name)
os.makedirs('../' + name)

print('copying files into the folder..')
os.system('cp -r * ../' + name)

print('removing any main.cpp and CMakeLists.txt file..')
if os.path.exists('../' + name + '/src/main.cpp'):
  os.system('rm ../' + name + '/src/main.cpp')
if os.path.exists('../' + name + '/src/CMakeLists.txt'):
  os.system('rm ../' + name + '/src/CMakeLists.txt')

print('removing any dummy Arduino .cpp/.h file..')
if os.path.exists('../' + name + '/src/Arduino.h'):
  os.system('rm ../' + name + '/src/Arduino.h')
if os.path.exists('../' + name + '/src/Arduino.cpp'):
  os.system('rm ../' + name + '/src/Arduino.cpp')

print('removing any build folder..')
if os.path.exists('../' + name + '/src/build'):
  os.system('rm -rf ../' + name + '/src/build')

print('removing python script..')
script_name = os.path.basename(__file__)
os.system('rm ../' + name + '/' + script_name)

print('removing any existing zip file of the target name..')
if os.path.exists('../' + name + '.zip'):
      os.system('rm ../' + name + '.zip')

print('creating zip file..')
os.system('cd ../ && zip -r ' + name + '.zip ' + name)

print('done! Now you can use the Arduino IDE and click on')
print('  Sketch->Include Library->Add .ZIP file')

