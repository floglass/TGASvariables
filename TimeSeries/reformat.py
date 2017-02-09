#! usr/bin/env python

"""
Florian Glass, 2017
Reformat database output file of GAIA photometric time series for local use.
"""

import numpy as np
import json


def reformat(line):  # band as in measurement band
    """
    reformat (long string) line into a list of appropriate format (strings and floats).
    get the four lines separately:
     line 0 - source_id, band_name, parallax, parallax_error
     line 1 - times
     line 2 - mag
     line 3 - mag_error

    :param line: the string to be reformatted.
    :returns: list of 4 lists
    """
    header = line[0]
    time_series = line[1]  # data
    
    ts = time_series.split('{')

    ts[3] = ts[3].rstrip("}(\"\")\"\"\n")
    ts[2] = ts[2].rstrip("}(\"\")\"\",\n")
    ts[1] = ts[1].rstrip("}(\"\")\"\",\n")
    ts[0] = header.rstrip("}(\"\")\"\",\n")

    ts[0] = [float(x) if i > 1 else x for i, x in enumerate(ts[0].split(','))]
    # last two elements are parallaxe and parallaxe_error: should be floats
    for i in range(1, 4):
        ts[i] = [float(x) for x in ts[i].split(',')]
    return ts


def checkdata(band, nodatacount):

    line = band.split('(')  # two lines - header and data
    if line[1][1] == ',':
        print "---------- NO DATA ----------"
        nodatacount += 1
        line = [line[0].split(','), [np.nan], [np.nan], [np.nan]]  # header, times, mag, mag_e
    else:
        line = reformat(line)
    return line, nodatacount


def getdata(target=None):
    """
    Load and read 'target', and call reformat on each line.
    Will output a dictionary (of dictionary of dictionary) of the form:
                                           +---------+
                                           |source_id|
                                           +----+----+
                                                |
              +----------------------+----------+-----------+---------------------+
              |                      |                      |                     |
              |                      |                      |                     |
              |                      |                      |                     |
              |                      |                      |                     |
           +--+--+                +--+--+                +--+--+              +---+----+
       +---+Band1+---+        +---+Band2+---+        +---+Band3+---+        +-+Parallax+-+
       |   +--+--+   |        |   +--+--+   |        |   +--+--+   |        | +--------+ |
       |      |      |        |      |      |        |      |      |        |            |
       |      |      |        |      |      |        |      |      |        |            |
    +--+--+ +-+-+ +--+--+  +--+--+ +-+-+ +--+--+  +--+--+ +-+-+ +--+--+  +--+-----+ +----+-----+
    |times| |mag| |mag_e|  |times| |mag| |mag_e|  |times| |mag| |mag_e|  |parallax| |parallax_e|
    +-----+ +---+ +-----+  +-----+ +---+ +-----+  +-----+ +---+ +-----+  +--------+ +----------+
    Flow chart from asciiflow.com

    :param target: the -unformatted- data file from GAIA DB. Each row is a specific band
    :returns: A dictionary containing the formatted data.
    """
    if target is None:
        target = 'source_2085235243870073344.csv'
        
    with open(target) as datafile:
        data_original = datafile.readlines()

    data_reformatted = dict()
    source_id_data = dict()
    no_data_count = 0
    for i in range(1, len(data_original)):
        print "\nProcessing row %s" % i
        row, no_data_count = checkdata(data_original[i], no_data_count)

        header = row[0]
        source_id = header[0]
        phot_band = header[1]
        print "target ID %s, photometric band %s" % (source_id, phot_band)
        parallax = header[2]
        parallax_e = header[3]
        
        times = row[1]
        magnitudes = row[2]
        magnitudes_errors = row[3]

        parallax_dict = dict([("parallax", parallax), ("parallax_e", parallax_e)])
        phot_band_dict = dict([("times", times), ("magnitudes", magnitudes),
                               ("magnitudes_errors", magnitudes_errors)])
        source_id_data[phot_band] = phot_band_dict
        if i % 3 == 0:
            source_id_data["parallax"] = parallax_dict
            data_reformatted[source_id] = source_id_data
            source_id_data = dict()

    print "\n%s bands had no data" % no_data_count
    return data_reformatted

    
if __name__ == "__main__":
    data = getdata('targets.csv')
    json.dump(data, open("data.txt", 'w'))

    # Load with json.load(open("data.txt"))
    # strings are going to be represented as unicode.
