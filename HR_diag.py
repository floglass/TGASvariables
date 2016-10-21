#! /usr/bin/python

import pandas as pd
import numpy as np
import notoriousTEX as nTEX
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import os.path
import argparse
import cPickle


def import_data(catalog='xmatch_TGAS_Simbad.csv', params=None, nrows=None, delimiter=','):
    """ imports 'catalog', and creates a DataFrame containing the columns specified in 'params'.
    'catalog' is expected to be in the .csv format. """
    print "Loading %s and creating DataFrame.." % catalog
    df = pd.read_csv(catalog, delimiter=delimiter, header=0, usecols=params, nrows=nrows)
    print "..Done\n----------"
    return df


def create_tycho_id(df):
    df['tycho2_id'] = df.TYC1.astype(str).str.cat(df.TYC2.astype(str), sep='-').str.cat(df.TYC3.astype(str), sep='-')
    df = df.rename(columns={'HIP': 'hip'})
    return df


def reindex(df):
    df = df.set_index(['hip', 'tycho2_id'])
    return df


def data_process(df=None, cutoff=0.2, catalog='xmatch_TGAS_Simbad.csv'):
    """ select the data, add columns for plotting """
    print "Selecting objects.."
    print "Cutoff at relative parallax error of %s\n----------" % cutoff
    df['sigma_pi/pi'] = df.loc[:, 'parallax_error']/df.loc[:, 'parallax']

    # only take objects with relative parallax error < cutoff
    df = df.loc[df.loc[:, 'parallax']/df.loc[:, 'parallax_error'] > 1./cutoff]

    # add columns 'sigma_pi/pi (for cutoff), 'B_V' (B-V) and 'M_G' (absolute mag)
    print catalog
    if catalog == 'xmatch_TGAS_Simbad.csv':
        df = df.loc[(df['J'] < 11.) & (df['K'] < 11.)]
        print "min in J: %s" % np.max(df['J'])
        print "max in J: %s" % np.min(df['J'])
        df.insert(10, 'B_V', df.loc[:, 'B']-df.loc[:, 'V'])
        df.insert(10, 'J_K', df.loc[:, 'J']-df.loc[:, 'K'])
        df.insert(10, 'M_G', df.loc[:, 'phot_g_mean_mag']-5.*(np.log10(1000./df.loc[:, 'parallax'])-1.))
        df.insert(10, 'M_J', df.loc[:, 'J']-5.*(np.log10(1000./df.loc[:, 'parallax'])-1.))
        df.insert(10, 'M_K', df.loc[:, 'K']-5.*(np.log10(1000./df.loc[:, 'parallax'])-1.))

    if catalog == 'xmatch_TGAS_VSX.csv':
        df.insert(8, 'Var', 'var')
    print "%s objects selected" % len(df)
    print "..Done\n----------"
    return df


def get_pickle(get_pkl_file='simbad_mag_errors.pkl'):
    """ load and unpickle a pickled pandas.DataFrame 'pkl_file' """
    print "Opening %s and unpickling the DataFrame.." % get_pkl_file
    with open(get_pkl_file, 'r') as pkl_file:
        df = cPickle.load(pkl_file)
    print "..Done"
    return df


def merge_df(df, df_error, merge_on=None):  # avoid mutable default argument with if merge_on is None: default
    """ add two pandas.DataFrames together """
    if merge_on is None:
        merge_on = ['hip', 'tycho2_id']
    df = df.merge(df, df_error, on=merge_on)
    return df


def plot_full(df, list_variables, variable_types, x='J_K', y='M_J'):
    plt.ion()
    print "cutoff at %s" % args.cutoff
    plot_hr_diag(df, x=x, y=y)
    plt.colorbar()
    plot_variable_stars(list_variables, variable_types, x=x, y=y)
    print "Plotting '%s' vs. '%s'" % y, x
    nTEX.plotSettings()
    print "----------"
    plt.show()
    return


def plot_hr_diag(df, x='J_K', y='M_J', marker='.', color='b'):
    """ plotting of the background stars, making the actual HR diagram.
    The plot is a 2d histogram, for better readability. Only bins with at least 10 stars a shown. """
    plt.figure()
    print "Plotting background stars.."
    plt.hist2d(df[x].tolist(), df[y].tolist(), (200, 200), norm=LogNorm(), cmin=10, alpha=.5, marker=marker, c=color)
    plt.axis([-0.5, 1.5, -3., 10])
    plt.gca().invert_yaxis()
    plt.xlabel(r'$J-K$')
    plt.ylabel(r'$%s$' % y)
    print "..Done\n----------"
    return


def plot_variable_stars(variablesdf, variabletype=None, x='J_K', y='M_G', size=40):
    """ Parent function of get_variable_stars. Sequencially select 'variableTypes' variable stars and
    plot them on the HR diagram. """
    if variabletype is None:
        variabletype = ['CEP', 'BCEP', 'BCEPS', 'DSCT', 'SR', 'SRA', 'SRB', 'SRC', 'SRD', 'RR', 'RRAB', 'RRC', 'GDOR',
                        'SPB', 'M']
    markers = ['^', '^', '^', 'v', '<', '<', '<', '<', '<', '>', '>', '>', '*', 's', 'p']
    colors = ['b', 'b', 'b', 'y', 'r', 'r', 'r', 'r', 'r', 'c', 'c', 'c', 'm', 'g', 'w']
    for i in range(len(variabletype)):
        plt.scatter(variablesdf[x].loc[variablesdf.loc[:, 'Type'] == variabletype[i]], variablesdf[y].loc[variablesdf
                    .loc[:, 'Type'] == variabletype[i]], facecolor=colors[i], marker=markers[i], s=size)
        print "plotting %s as %s%s" % (variabletype[i], colors[i], markers[i])
    return


def get_variable_stars(df, df2, variabletype=None):
    """ Child function fo plot_variable_stars. Process the DataFrame to select only stars marked as
    'var_type' variable stars """
    if variabletype is None:
        variabletype = ['CEP', 'BCEP', 'BCEPS', 'DSCT', 'SR', 'SRA', 'SRB', 'SRC', 'SRD', 'RR', 'RRAB', 'RRC',
                        'GDOR', 'SPB', 'M']
    are_variables = df2[df2.loc[:, 'Type'].isin(variabletype)]
    variable_df = pd.merge(df, are_variables, how='inner', on=['hip', 'tycho2_id'])
    variable_df = variable_df[np.isnan(variable_df.loc[:, 'K']) is False]  # only get stars with K measurement
    return variable_df


def plot_errors():

    return


def simbad_makescript(df2):
    """ Creates a SIMBAD script to get photometric info -mainly errors- on a subset of stars (the variables) """
    print "Starting process.."
    with open('SIMBAD_script.txt', 'w') as sc:
        print "Opening 'SIMBAD_script.txt' and dumping variable stars names:"
        sc.write("format object form1 \"%IDLIST(hip,tyc),%FLUXLIST(B,V,J,K;N F,E,)\"\n")
        hip_numbers = df2['hip'].loc[np.isnan(df2.loc[:, 'hip']) is False].tolist()
        print "Hip identification.."
        for hip_number in hip_numbers:
            sc.write("query id hip %s\n" % (int(hip_number)))
        print "..Done"
        tycho_numbers = df2['tycho2_id'].loc[pd.isnull(df2.loc[:, 'tycho2_id']) is False].tolist()
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
            for count,string in enumerate(['HIP', 'TYC', 'B', 'V', 'J', 'K']):
                if string not in line:
                    # print "%s not in %s"%(string, line)
                    line_list.insert(count,np.nan)
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
    df = pd.DataFrame(cleaned_output[1:], columns=cleaned_output[0])
    df = df.drop('\n', axis=1)
    return df


def simbad_output():
    new_out, data = simbad_outputcolumn_clean()
    df = simbad_create_data_frame(new_out)
    pickle_it(df)
    return


def pickle_it(df, target_file='simbad_mag_errors.pkl'):
    """ pickle a DataFrame 'df' into 'target_file' """
    print "Starting pickling of %s.." % target_file
    with open(target_file, 'wb') as pkl_file:
        cPickle.dump(df, pkl_file)
    print "%s created\n----------" % target_file
    return


def un_pickle(target_file):
    """ unpickle a DataFrame from 'target_file' """
    print "Opening pickle file '%s'.." % target_file
    with open(target_file, 'rb') as pkl_file:
        df = cPickle.load(pkl_file)
    print "..Done"
    return df
# ================================================== #
# ================================================== #


if __name__ == "__main__":
    #################################
    # Parser options and arguments: #
    #################################
    parser = argparse.ArgumentParser(description='Creates a HR diagram with variables stars')
    parser.add_argument('-p', '--pickle', dest='pickle', action='store_false', default=True,
                        help="Switch for building pickle files of stars and variable stars. Default=%default.")
    parser.add_argument('-c', '--cutoff', dest='cutoff', action='store', default=0.2, type=float,
                        help="At which percentage in relative parallax error should the data stop being read."
                             "Default=%(default).")
    parser.add_argument('-r', '--nrows', dest='nrows', action='store', type=int, default=None,
                        help="The number of line contained in the initial DataFrame. Default=%default")
    parser.add_argument('-d', '--draw', dest='draw', action='store_false', default=True,
                        help="Switch the plotting off. Default=%default.")
    args = parser.parse_args()

    #########################################################
    # Beginning of the script -- data import and selection: #
    #########################################################
    target = 'df_%s.pkl' % args.cutoff
    target2 = 'df2_%s.pkl' % args.cutoff
    if os.path.isfile(target) and os.path.isfile(target2):
        print "Opening pickle file '%s'.." % target
        with open(target, 'rb') as pkl_file:
            df = cPickle.load(pkl_file)
        print "..Done"
        print "Opening pickle file '%s'.." % target2
        with open(target2) as pkl_file:
            df2 = cPickle.load(pkl_file)
        print "..Done"
    else:
        # Use of 'xmatch_TGAS_Tycho2.csv' instead of 'xmatch_TGAS_Simbad.csv'
        # want to use tycho2's B and V mags
        df = import_data(catalog='xmatch_TGAS_Simbad.csv', params=('hip', 'tycho2_id', 'parallax', 'parallax_error',
                                                                   'phot_g_mean_mag', 'B', 'V', 'J', 'K'),
                         nrows=args.nrows)
        df = data_process(df, catalog='xmatch_TGAS_Simbad.csv', cutoff=args.cutoff)

        df2 = import_data(catalog='xmatch_TGAS_VSX.csv', params=('hip', 'tycho2_id', 'parallax', 'parallax_error',
                                                                 'Name', 'V', 'Type'), nrows=args.nrows)
        df2 = data_process(df2, catalog='xmatch_TGAS_VSX.csv', cutoff=args.cutoff)

    variable_types = ['CEP', 'BCEP', 'BCEPS', 'DSCT', 'SR', 'SRA', 'SRB', 'SRC', 'SRD', 'RR', 'RRAB', 'RRC',
                      'GDOR', 'SPB', 'M']
    list_variables = get_variable_stars(df, df2, variable_types)
    ########################################################
    # Plotting of HR diag and variables stars with errors: #
    ########################################################
    if args.draw is True:
        plot_full(df, list_variables, variable_types, x='J_K', y='M_J')

    #####################################################
    # Creation of the pickle file                       #
    # Check that the pickle file isn't already present: #
    #####################################################
    if args.pickle is True:
        print "Pickle switch is ON:"
        print "Checking if pickle files are present.."
        if os.path.isfile(target) and os.path.isfile(target2):
            print ".. Target files already here!"
            print "File creation aborted\n----------"
        else:
            print ".. Files not found"
            print "Creating new pickle files '%s' and '%s'" % (target, target2)
            pickle_it(df, target)
            pickle_it(df2, target2)
            print "..Pickle creation done\n----------"
