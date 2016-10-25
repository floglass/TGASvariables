#! /usr/bin/python

import pandas as pd
import numpy as np
import notoriousTEX as nTEX
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import os.path
import argparse
import cPickle  # use JSON instead ?

""" This code's goal is to add some photometric measurements from the Tycho2 catalog to the recently released Tycho-Gaia
DR1 catalog. The latter not including the photometric data.
The Tycho2 catalog contains 2539913 stars with BT and VT magnitudes, as well as other astrometric measurements.
The TGAS catalog contains 2057050 stars with very accurate astrometric measurements, as well as a g band mean mag.
We are trying to make a M_V versus B-V HR-diagram, to assess the accuracy of the catalog on the variable stars sample.
"""


def import_data(catalog='xmatch_TGAS_Simbad.csv', params=None, nrows=None, delimiter=','):
    """ imports 'catalog', and creates a pandas.DataFrame containing the columns specified in 'params'.
    'catalog' is expected to be in the .csv format. """
    print "Loading %s and creating DataFrame.." % catalog
    df_imported = pd.read_csv(catalog, delimiter=delimiter, header=0, usecols=params, nrows=nrows)
    print "..Done\n----------"
    return df_imported


def create_tycho_id(tycho2df):
    """ creates a new column 'tycho2_id' in the tycho2 catalog. This is for comparison with the TGAS catalog. """
    tycho2df['tycho2_id'] = tycho2df.TYC1.astype(str).str.cat(tycho2df.TYC2.astype(str), sep='-')\
        .str.cat(tycho2df.TYC3.astype(str), sep='-')
    tycho2df = tycho2df.rename(columns={'HIP': 'hip'})
    return tycho2df


def reindex(tycho2df_toindex):
    """ passes the two columns 'hip' and 'tycho2_id' as indexes of the DataFrame """
    tycho2df_toindex = tycho2df_toindex.set_index(['hip', 'tycho2_id'])
    return tycho2df_toindex


def data_process(df_toprocess=None, cutoff=0.2, bv_cutoff=0.15, catalog=None):
    """ select data with relative parallax error less than 'cutoff', add absolute magnitude columns for plotting """

    print "Selecting objects.."
    df_toprocess['sigma_pi/pi'] = df_toprocess.loc[:, 'parallax_error'].astype(float) / df_toprocess.loc[:, 'parallax']\
        .astype(float)
    print "..Done\nCutoff at relative parallax error of %s\n----------" % cutoff

    # only take objects with relative parallax error < cutoff
    df_toprocess = df_toprocess.loc[df_toprocess.loc[:, 'parallax'] /
                                    df_toprocess.loc[:, 'parallax_error'] > 1. / cutoff]

    print catalog
    if catalog is None:
        print "Replacing whitespace by nan"
        df_toprocess = df_toprocess.replace('      ', np.nan)  # some cells are '      ' instead of nan

        print "Converting BTmag and VTmag to floats.."
        df_toprocess.BTmag = df_toprocess.BTmag.astype(float)
        df_toprocess.VTmag = df_toprocess.VTmag.astype(float)
        # Some values are NaN:
        print "Removing objects with missing BT or VT measurements.."
        df_toprocess = df_toprocess[df_toprocess.BTmag.notnull()]
        df_toprocess = df_toprocess[df_toprocess.VTmag.notnull()]

        print "Computing B-V and M_V.."
        df_toprocess['B_V'] = df_toprocess.BTmag - df_toprocess.VTmag
        df_toprocess['M_V'] = df_toprocess.VTmag - 5. * (np.log10(1000. / df_toprocess.parallax) - 1.)

        print "Converting sigma BT and sigma VT to float.."
        df_toprocess.e_BTmag = df_toprocess.e_BTmag.astype(float)
        df_toprocess.e_VTmag = df_toprocess.e_VTmag.astype(float)

        print "Computing sigma B-V.."
        df_toprocess['e_B_V'] = np.sqrt(df_toprocess.e_BTmag.pow(2)+df_toprocess.e_VTmag.pow(2))

        print "Applying selection on sigma BT-VT < %s.." % bv_cutoff
        df_toprocess = df_toprocess[df_toprocess.e_B_V < bv_cutoff]

    if catalog == 'xmatch_TGAS_Simbad.csv':
        df_toprocess = df_toprocess.loc[(df_toprocess['J'] < 11.) & (df_toprocess['K'] < 11.)]
        print "min in J: %s" % np.max(df_toprocess['J'])
        print "max in J: %s" % np.min(df_toprocess['J'])
        df_toprocess.insert(10, 'B_V', df_toprocess.loc[:, 'B'] - df_toprocess.loc[:, 'V'])

        df_toprocess.insert(10, 'J_K', df_toprocess.loc[:, 'J'] - df_toprocess.loc[:, 'K'])
        df_toprocess.insert(10, 'M_G', df_toprocess.loc[:, 'phot_g_mean_mag'] - 5. *
                            (np.log10(1000. / df_toprocess.loc[:, 'parallax']) - 1.))
        df_toprocess.insert(10, 'M_J', df_toprocess.loc[:, 'J'] - 5. *
                            (np.log10(1000. / df_toprocess.loc[:, 'parallax']) - 1.))
        df_toprocess.insert(10, 'M_K', df_toprocess.loc[:, 'K'] - 5. *
                            (np.log10(1000. / df_toprocess.loc[:, 'parallax']) - 1.))

    if catalog == 'xmatch_TGAS_VSX.csv':
        df_toprocess.insert(8, 'Var', 'var')
    print "%s objects selected" % len(df_toprocess)
    print "..Done\n----------"
    return df_toprocess


def get_pickle(get_file='simbad_mag_errors.pkl'):
    """ load and unpickle a pickled pandas.DataFrame 'pkl_file' """
    print "Opening %s and unpickling the DataFrame.." % get_file
    with open(get_file, 'r') as opened_file:
        df_unpickled = cPickle.load(opened_file)
    print "..Done"
    return df_unpickled


def merge_df(merge_on_df, merge_with_df, merge_column=None):
    """ add two pandas.DataFrames together on columns 'hip' and 'tycho2_id' columns.
    Avoid mutable default argument in 'merge_column' (if merge_column is None: default)"""
    if merge_column is None:
        merge_column = ['hip', 'tycho2_id']
    merge_on_df = merge_on_df.merge(merge_on_df, merge_with_df, on=merge_column)
    return merge_on_df


def plot_full(plot_df, list_variable_stars, variable_stars_types=None, x='B_V', y='M_V', cutoff=0.2, bvcutoff=0.05):
    """
    plots the full thing : HR diagram with TGAS targets, and highlighted variable stars

    :param plot_df: pandas DataFrame to plot
    :param list_variable_stars: pandas DataFrame of variable stars
    :param variable_stars_types: list of variable stars types
    :param x: abscissa of the graph
    :param y: ordinate of the graph
    :return:
    """
    if variable_stars_types is None:
        variable_stars_types = ['CEP', 'BCEP', 'BCEPS', 'DSCT', 'SR', 'SRA', 'SRB', 'SRC', 'SRD', 'RR', 'RRAB', 'RRC',
                                'GDOR', 'SPB', 'M']
    plt.ion()
    print "cutoff at %s" % cutoff
    print "Plotting '%s' vs. '%s'" % (y, x)
    plot_hr_diag(plot_df, x=x, y=y, cutoff=cutoff, bvcutoff=bvcutoff)
    plt.colorbar()
    plot_variable_stars(list_variable_stars, variable_stars_types, x=x, y=y)
    nTEX.plotSettings()
    print "----------"
    plt.show()
    return


def plot_hr_diag(hr_df, x='B_V', y='M_V', cutoff=0.2, bvcutoff=0.05):
    """ plotting of the background stars (HR diagram).
    The plot is a 2d histogram, for better readability. Only bins with at least 10 stars a shown. """
    plt.figure()
    print "Plotting background stars.."
    plt.hist2d(hr_df[x].tolist(), hr_df[y].tolist(), (200, 200), norm=LogNorm(), cmin=10, alpha=.5)
    plt.axis([-0.5, 2., -2., 8.])
    plt.gca().invert_yaxis()
    plt.xlabel(r'$B-V$')
    plt.ylabel(r'$%s$' % y)
    plt.suptitle(r'cutoff = %s, $\sigma_{B-V}< %s$' % (cutoff, bvcutoff))
    print "..Done\n----------"
    return


def plot_variable_stars(variablesdf, variabletype=None, x='B_V', y='M_V', size=40):
    """
    Parent function of get_variable_stars. Sequencially select 'variableTypes' variable stars and
    plot them on the HR diagram.

    :param variablesdf: pandas DataFrame containing the variable stars
    :param variabletype: list of names of the variable stars' types
    :param x: abscissa of the graph
    :param y: ordinate of the graph
    :param size: size of markers on graph
    :return:
    """
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


def get_variable_stars(df_data, df_variables_names, variabletype=None):
    """ Child function fo plot_variable_stars. Process the DataFrame to select only stars marked as
    'var_type' variable stars. """
    if variabletype is None:
        variabletype = ['CEP', 'BCEP', 'BCEPS', 'DSCT', 'SR', 'SRA', 'SRB', 'SRC', 'SRD', 'RR', 'RRAB', 'RRC',
                        'GDOR', 'SPB', 'M']

    print "Selecting variable stars.."
    are_variables = df_variables_names[df_variables_names.loc[:, 'Type'].isin(variabletype)]
    typesdf = are_variables[['hip', 'tycho2_id', 'Type']]
    print "..Done"
    print "Preparing subselection of initial DataFrame.."
    print "..Making Hipparcos list.."
    hip_list = are_variables.hip.tolist()
    hip_list = np.array(hip_list)
    hip_list = hip_list[~np.isnan(hip_list)]  # remove the nans
    hip_list = list(hip_list)
    print "..Making Tycho2 list.."
    tycho2_list = are_variables.tycho2_id.tolist()
    tycho2_list = np.array(tycho2_list)
    tycho2_list = tycho2_list[tycho2_list != 'nan']  # tycho2 is str
    tycho2_list = list(tycho2_list)
    print "..Done\n----------"

    print "Getting Hipparcos and Tycho variable objects.."
    hip_objects = df_data[df_data.hip.isin(hip_list)]
    hip_objects = pd.merge(hip_objects, typesdf, on='hip', how='inner')
    hip_objects = hip_objects.drop('tycho2_id_y', axis=1)
    hip_objects = hip_objects.rename(columns={'hip_x': 'hip', 'tycho2_id_x': 'tycho2_id'})

    tycho_objects = df_data[df_data.tycho2_id.isin(tycho2_list)]
    tycho_objects = pd.merge(tycho_objects, typesdf, on='tycho2_id', how='inner')
    tycho_objects = tycho_objects.drop('hip_y', axis=1)
    tycho_objects = tycho_objects.rename(columns={'hip_x': 'hip', 'tycho2_id_x': 'tycho2_id'})
    print "..Done\n----------"

    variable_df = pd.concat([hip_objects, tycho_objects], axis=0, ignore_index=True)
    return variable_df


def plot_errors():
    """ plots the error bars on the variable stars -- work in progress """
    return


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


def get_hip_sources(df_main_catalog=None, variable_class=None):
    if variable_class is None:
        variable_class = ['CEP', 'BCEP', 'BCEPS', 'DSCT', 'SR', 'SRA', 'SRB', 'SRC', 'SRD', 'RR', 'RRAB', 'RRC',
                          'GDOR', 'SPB', 'M']
    hip_sources = pd.read_csv('Hip/variables.tsv', delimiter=';', header=0)
    variable_stars = hip_sources[hip_sources.loc[:, 'VarType'].isin(variable_class)]
    variable_stars.hip = variable_stars.hip.astype(float)

    hip_objects = df_main_catalog[df_main_catalog.hip.isin(variable_stars.hip)]
    hip_objects = pd.merge(hip_objects, variable_stars, on='hip', how='inner')
    hip_objects = hip_objects.rename(columns={'VarType': 'Type'})

    return hip_objects
# ================================================== #
# ================================================== #


if __name__ == "__main__":
    #################################
    # Parser options and arguments: #
    #################################
    parser = argparse.ArgumentParser(description='Creates a HR diagram with variables stars',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-p', '--pickle', dest='pickle', action='store_false', default=True,
                        help="Switch for building pickle files of stars and variable stars. Default=default.")
    parser.add_argument('-c', '--cutoff', dest='cutoff', action='store', default=0.2, type=float,
                        help="At which percentage in relative parallax error should the data stop being read."
                             "Default=default.")
    parser.add_argument('-r', '--nrows', dest='nrows', action='store', type=int, default=None,
                        help="The number of line contained in the initial DataFrame. Default=default")
    parser.add_argument('-d', '--draw', dest='draw', action='store_false', default=True,
                        help="Switch the plotting off. Default=default.")
    parser.add_argument('-m', '--magnitude', dest='bvcutoff', action='store', default=0.15, type=float,
                        help="cutoff in 'B-V' magnitude. Default=default.")
    args = parser.parse_args()

    #########################################################
    # Beginning of the script -- data import and selection: #
    #########################################################
    # target = 'df_%s.pkl' % args.cutoff
    target = 'xmatch_TGAS_Tycho2_ByHand.pkl'  # TGAS + photometric measurements from Tycho2
    # target2 = 'df2_%s.pkl' % args.cutoff  # comes from xmatch TGAS and VSX

    print "Opening pickle file '%s'.." % target
    df = pd.read_pickle(target)
    print "..Done"
    # print "Opening pickle file '%s'.." % target2
    # df2 = pd.read_pickle(target2)
    # print "..Done"
    df = data_process(df, cutoff=args.cutoff, bv_cutoff=args.bvcutoff)
    """
    else:
        # Use of 'xmatch_TGAS_Tycho2.csv' instead of 'xmatch_TGAS_Simbad.csv'
        # want to use tycho2's B and V mags

        df = import_data(catalog='xmatch_TGAS_Simbad.csv', params=('hip', 'tycho2_id', 'parallax', 'parallax_error',
                                                                   'phot_g_mean_mag', 'B', 'V', 'J', 'K'),
                         nrows=args.nrows)
        df = data_process(df, catalog='xmatch_TGAS_Simbad.csv', cutoff=args.cutoff)
    """
    df2 = import_data(catalog='xmatch_TGAS_VSX.csv', params=('hip', 'tycho2_id', 'parallax', 'parallax_error',
                                                             'Name', 'V', 'Type'), nrows=args.nrows)
    df2 = data_process(df2, catalog='xmatch_TGAS_VSX.csv', cutoff=args.cutoff)

    variable_types = ['CEP', 'BCEP', 'BCEPS', 'DSCT', 'SR', 'SRA', 'SRB', 'SRC', 'SRD', 'RR', 'RRAB', 'RRC',
                      'GDOR', 'SPB', 'M']

    list_variables = get_variable_stars(df, df2, variable_types)
    ########################################################
    # Plotting of HR diag and variables stars (with errors): #
    ########################################################
    if args.draw is True:
        plot_full(df, list_variables, variable_types, x='B_V', y='M_V', cutoff=args.cutoff, bvcutoff=args.bvcutoff)
