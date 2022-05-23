# This is version Python 3 compatible with parallel GPS logging and for the 
# modern Edgetech deck box serial output.
# Provided courtesy of James R McVey (PNNL) on 5/23/2022
import os
import sys
import time
import serial
import pynmea2
import threading
from glob import glob
from datetime import datetime


def getSerialDevice(prompt, default=1):
    """Look through the /dev directory looking for serial devices on iOS/linux

    Args:
        prompt (str): user input display prompt
        default (int, optional): index of serial device. Defaults to 1.

    Returns:
        str: The path to device or None
    """

    # Searches for Return a list of all available keyspan serial devices on the Mac
    deviceList = glob('/dev/tty.USA*')
    if len (deviceList) == 0:
        return None
    
    # Print input prompt
    print(prompt)

    for i, dev in enumerate(deviceList):
        print("[%d] %s" % (i+1, dev))

    try:
        value = int(input("\nChoose [default=%d]: " % default))
    except:
        value = None
        
    if not value or value > len(deviceList):
        value = default
        print("Using Default (%s)" % deviceList[value-1])
    else:
        print("Using %s for your serial device" % deviceList[value-1])
        
    return deviceList[value-1]


def gps_thread(ser, out):
    """Reads GPGGA strings incoming from GPS serial output

    Args:
        ser (serial object): open GPS serial port

    Returns:
        pnynmea object: parsed GPGGA NMEA string
    """
    # Make sure to start reading buffer at new NMEA string
    ser.read_until('\n'.encode())

    # Read NMEA strings and pull out GGA
    while True:
        newline=ser.readline().decode('utf-8')

        if 'GGA' in newline:
            gga=pynmea2.parse(newline,check=True)

            out['nmea'] = gga


def main():
    # Determine serial ports to use for the devices
    if 'win' in sys.platform:
        gpsDev = input("Enter COM port of GPS serial device: ")
        edgeDev = input("Enter COM port of EdgeTech Ranging System: ")
    else:
        gpsDev = getSerialDevice("Enter id of GPS serial device: ", 1)
        edgeDev = getSerialDevice("Enter id of EdgeTech Ranging System: ", 2)

    if (gpsDev == None or edgeDev == None):
        print("Error: Cannot find a serial device.\nIs the Keyspan plugged into the USB port?\n")
        sys.exit(-1)

    # Got something...
    print("GPS:      %s" % gpsDev)
    print("EdgeTech: %s" % edgeDev)
    ## For 8011M Acoustic Deck Box
    # print("\n\n")
    # print("=" * 70)
    # print("Please configure the EdgeTech to the following setup:")
    # print("1) Set the units for ranging to be in msec (instead of meters).")
    # print("   Keypad Commands: <MENU> -> 1 (Range Setup) -> 5 (Units) -> msec")
    # print("2) Configure the EdgeTech to send ranging information to the serial port.")
    # print("   Keypad Commands: <MENU> -> 3 (RS232 Setup) -> 3 (Logging) -> Ranges")
    # print("3) Configure the repetion rate for ranging (set to 20 seconds).")
    # print("   Keypad Commands: <R> -> <MAN/REP> twice -> type 20.0 -> <ENT>")
    # print("4) Start the Ranging.")
    # print("   Keypad Commands: <R> -> <INT>")
    # print("=" * 70)

    try:
        gps = serial.Serial(gpsDev, 9600, timeout=1)
        edgeTech = serial.Serial(edgeDev, 9600, timeout=None)
    except:
        print("\nError: Cannot open communications with one of the serial devices.")
        sys.exit(-1)

    # Flush the buffers --- had a problem with residual data
    gps.flush()
    edgeTech.flush()

    # Ask for a filename
    filename = input("File name (without extension): ")
    cruiseID = input("Cruise ID: ")
    siteID       = input("Site Identifier: ")
    instrumentID = input("Instrument Identifier:  ")
    latDropPoint = input("Drop point (Latitude):  ")
    lonDropPoint = input("Drop point (Longitude): ")
    depthPoint   = input("Initial Depth (meters): ")
    comment      = input("Comment: ")

    # Create a distinctive filename (based on the current time) and open it for writing.
    filename = os.path.expanduser("~/ranging-" + str(time.time()).split('.')[0] )
    print("Output File:   %s" % filename)

    with open(filename, "w") as fp:
        fp.write("Ranging data taken on:  %s\n" % str(datetime.now()))
        fp.write("Cruise:                 %s\n" % cruiseID)
        fp.write("Site:                   %s\n" % siteID)
        fp.write("Instrument:             %s\n" % instrumentID)
        fp.write("Drop Point (Latitude):  %s\n" % latDropPoint)
        fp.write("Drop Point (Longitude): %s\n" % lonDropPoint)
        fp.write("Depth (meters):         %s\n" % depthPoint)
        fp.write("Comment:                %s\n" % comment)
        fp.write("="*50 + "\n\n")
        out = dict.fromkeys(['nmea']) # gps output saved here
        print("\n\nWaiting for ranging data from EdgeTech...")

        # Start thread to read GPS input
        x = threading.Thread(target=gps_thread, args=(gps,out), daemon=True)
        x.start()

        # Start main loop
        # Files will automatically close when the program ends.
        while True:
            # Wait for EdgeTech serial data
            # P.A.C.S Acoustic Deck Box: b'TX = 11.0 RX = 12.0 time = 00.090 Sec.'
            range = edgeTech.readline().decode().split(' ')

            # Pull time and GPS once Edgetech data is recorded
            dt = datetime.now().timetuple()
            currentTime = str(dt.tm_year)+':'+str(dt.tm_yday)+':'+str(dt.tm_hour)+':'+str(dt.tm_min)+':'+str(dt.tm_sec)

            # Read current GPS point
            pos = out['nmea']

            try:
                distance = float(range[-2])*1000 # convert s to ms
                units = 'ms'

                s = "%5d %s Lat: %s %s %s  Lon: %s %s %s  Alt: %s Time(UTC): %s" % (
                    int(distance),
                    units,
                    int(pos.lat[0:2]), #latDeg,
                    float(pos.lat[2:10]), #latMin,
                    pos.lat_dir, #latHem,
                    int(pos.lon[0:3]), #longDeg,
                    float(pos.lon[3:10]), #longMin,
                    pos.lon_dir, #longHem,
                    pos.altitude, #alt,
                    currentTime
                    )
            except:
                s = "Event skipped - Timeout or Badly formatted data was received"
                
            print(s)
            fp.write(s + '\n')


if __name__ ==  '__main__':
    main()
