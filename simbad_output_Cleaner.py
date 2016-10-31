#! /usr/bin/python

import os.path
import sys
import numpy as np
import pandas as pd
import cPickle


def simbad_makescript(df_variable_stars_list):
    """ Creates a SIMBAD script to get photometric info -mainly errors- on a subset of stars (the variables) """
    print "Starting process.."
    with open('SIMBAD_script.txt', 'w') as sc:
        print "Opening 'SIMBAD_script.txt' and dumping variable stars names:"
        sc.write("format object form1 \"%IDLIST(hip,tyc),%FLUXLIST(B,V,J,K;N F,E,)\"\n")
        hip_numbers = df_variable_stars_list['hip'].loc[np.isnan(df_variable_stars_list.loc[:, 'hip'])
                                                        is False].tolist()
        print "Hip identification.."
        for hip_number in hip_numbers:
            sc.write("query id hip %s\n" % (int(hip_number)))
        print "..Done"
        tycho_numbers = df_variable_stars_list['tycho2_id'].loc[pd.isnull(df_variable_stars_list.loc[:, 'tycho2_id'])
                                                                is False].tolist()
        print "Tycho2 identification.."
        for tycho_number in tycho_numbers:
            sc.write("query id tyc %s\n" % tycho_number)
        print "..Done"
        print "'SIMBAD_script.txt' created\n----------"
    return


def simbad_outputcolumn_clean(output_file='simbad_output_0.1.txt'):
    """ Go through a SIMBAD output file and checks if a column is missing, by splitting every row into a list and
    adding missing data as NaN as well as replacing '~' with NaN """
    print "Opening %s.." % output_file
    with open(output_file) as dataFile:
        data = dataFile.readlines()
    print "..Done"

    new_output = list()
    print "Processing missing data and cleaning.."
    for i in range(len(data)):
        line = data[i]
        line_list = line.split(',')
        # handling '~' : meaning an error measurement is missing
        if '~' in line:
            line = line.replace('~', str(np.nan))
            line_list = line.split(',')
        if len(line_list) < 11:
            # handling missing columns : B, V, J, K are actually two columns
            for count, string in enumerate(['HIP', 'TYC', 'B', 'V', 'J', 'K']):
                if string not in line:
                    # print "%s not in %s"%(string, line)
                    line_list.insert(count, np.nan)
                    if string in ['B', 'V', 'J', 'K']:
                        line_list.insert(count+1, np.nan)
        if i > 0:
            line_list = simbad_outputlineclean(line_list)
        new_output.append(line_list)
    print "..Done\n----------"
    return new_output, data


def simbad_outputlineclean(line_list):
    """ go through the line, and remove the leading string as well as white spaces """
    for i in range(len(line_list)):
        for string in ['HIP', 'TYC', 'B', 'V', 'J', 'K']:
            # print i, string
            line_list[i] = str(line_list[i]).lstrip(string)
            line_list[i] = str(line_list[i]).lstrip()
    return line_list


def simbad_create_data_frame(cleaned_output):
    """ uses the output from 'simbad_outputcolumn_clean' to create a pandas.DataFrame """
    simbad_df = pd.DataFrame(cleaned_output[1:], columns=cleaned_output[0])
    simbad_df = simbad_df.drop('\n', axis=1)
    return simbad_df


def simbad_output():
    """ process a simbad output file into a workable pickled pandas.DataFrame """
    new_out, data = simbad_outputcolumn_clean()
    dataframe = simbad_create_data_frame(new_out)
    pickle_it(dataframe)
    return


def pickle_it(df_to_pickle, target_file='simbad_mag_errors.pkl'):
    """ pickle a DataFrame 'df_to_pickle' into 'target_file' """
    print "Starting pickling of %s.." % target_file
    with open(target_file, 'wb') as open_file:
        cPickle.dump(df_to_pickle, open_file)
    print "%s created\n----------" % target_file
    return


def un_pickle(target_file):
    """ unpickle a DataFrame from 'target_file' """
    print "Opening pickle file '%s'.." % target_file
    with open(target_file, 'rb') as open_file:
        df_fresh = cPickle.load(open_file)
    print "..Done"
    return df_fresh


#==================================================#
#==================================================#

if __name__=='__main__':
    new_out, data = columnCheck()
    df = createDataFrame(new_out)
    pickleIt(df)
