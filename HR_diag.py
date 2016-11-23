#! /usr/bin/python

import pandas as pd
import numpy as np
import notoriousTEX as nTEX
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import argparse
import cPickle  # use JSON instead ?

""" This code's goal is to add some photometric measurements from the Tycho2 catalog to the recently released Tycho-Gaia
DR1 catalog. The latter not including the photometric data.
The Tycho2 catalog contains 2539913 stars with BT and VT magnitudes, as well as other astrometric measurements.
The TGAS catalog contains 2057050 stars with very accurate astrometric measurements, as well as a g band mean mag.
We are trying to make a M_VT versus BT-VT Hertzprung-Russell diagram, to check variable stars in the catalog sample.
The photometry used is Tycho2's, hence the BT, VT and M_VT magnitudes (instead of the more traditional Johnson's B, V)
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
    """ passes the two columns 'hip' and 'tycho2_id' as indexes of the DataFrame
    -- DEPRECATED --"""
    tycho2df_toindex = tycho2df_toindex.set_index(['hip', 'tycho2_id'])
    return tycho2df_toindex


def data_process(df_toprocess=None, cutoff=0.2, bv_cutoff=0.15, catalog=None):
    """
    select data with relative parallax error less than 'cutoff', add absolute magnitude columns for plotting.
    If catalog is not None, the cutoff on B-V will not be applied (ensures initial variable stars DataFrame is not
    constrained in magnitudes)

    :param df_toprocess: pandas.DataFrame
    :param cutoff: the maximum relative error on parallax allowed
    :param bv_cutoff: the maximum error on temperature (B-V mag)
    :param catalog: catalog name from which the stars are taken, in a .csv form
    :return: processed DataFrame
    """

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
        df_toprocess = df_toprocess[df_toprocess.V == 0]
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
    :type plot_df: pandas.core.frame.DataFrame
    :param list_variable_stars: pandas DataFrame of variable stars
    :type list_variable_stars: pandas.core.frame.DataFrame
    :param variable_stars_types: list of variable stars types
    :param x: abscissa of the graph
    :param y: ordinate of the graph
    :param cutoff: max sigma_parallax / parallax allowed on stars
    :param bvcutoff: max sigma_B-V allowed on stars
    :return:
    """
    if variable_stars_types is None:
        variable_stars_types = ['CEP', 'BCEP', 'BCEPS', 'DSCT', 'SR', 'SRA', 'SRB', 'SRC', 'SRD', 'RR', 'RRAB', 'RRC',
                                'GDOR', 'SPB', 'M', 'LPV']
    plt.ion()
    print "Plot started.."
    print "Cutoff at %s" % cutoff
    print "Plotting '%s' vs. '%s'" % (y, x)
    plot_hr_diag(plot_df, x=x, y=y, cutoff=cutoff, bvcutoff=bvcutoff)
    cb = plt.colorbar()
    cb.set_label("number of TGAS sources per color-mag bin")
    plot_variable_stars(list_variable_stars, variable_stars_types, x=x, y=y)
    nTEX.plotSettings(width=11., height=10.)
    print "..Done\n----------"
    plt.show()
    return


def plot_hr_diag(hr_df, x='B_V', y='M_V', cutoff=0.2, bvcutoff=0.05):
    """ plot the background stars (HR diagram).
    The plot is a 2d histogram, for better readability. Only bins with at least 10 stars a shown.

    :param hr_df: pandas.DataFrame of the background stars
    :param x: hr_df column name to put on the x-axis
    :param y: hr_df column name for y-axis
    :param cutoff: the relative parallax error cutoff
    :param bvcutoff: the B-V magnitude cutoff
    :return:
    """
    plt.figure(figsize=(11., 10.))
    print "Plotting background stars.."
    plt.set_cmap('gray_r')
    plt.hist2d(hr_df[x].tolist(), hr_df[y].tolist(), (200, 200), norm=LogNorm(), cmin=10)
    plt.axis([-0.2, 2.35, -3., 7.])
    plt.gca().invert_yaxis()
    plt.xlabel(r'$BT-VT$ (mag)')
    plt.ylabel(r'$M_{VT}$ (mag)')  # Plotting M_{VT}
    plt.title(r'$\sigma_\pi / \pi < %s, \sigma_{BT-VT}< %s$ mag' % (cutoff, bvcutoff))
    print "..Done"
    return


def plot_variable_stars(variablesdf, variabletype=None, x='B_V', y='M_V'):
    """
    Parent function of get_variable_stars. Sequencially select 'variableTypes' variable stars and
    plot them on the HR diagram.

    :param variablesdf: pandas DataFrame containing the variable stars
    :param variabletype: list of names of the variable stars' types
    :param x: abscissa of the graph
    :param y: ordinate of the graph
    :return:
    """
    if variabletype is None:
        variabletype = ['CEP', 'BCEP', 'BCEPS', 'DSCT', 'SR', 'SRA', 'SRB', 'SRC', 'SRD', 'RR', 'RRAB', 'RRC', 'GDOR',
                        'SPB', 'M', 'LPV']
    markers = ['^', 'D', 'D', 'v', 's', 'D', 'D', 'D', 'D', 's', 'D', 'D', 'D', 'o', 'p', 'o']
    colors = ['k', 'k', 'k', '#00c000', 'r', 'r', 'r', 'r', 'r', 'm', 'm', 'm', '#00c0ff', (1, .7, 0), 'w', 'w']
    sizes = [50,  40,  40,  40,  50,  40,  40,  40,  40,  50,  50,  50,  40,  40,  45, 40]
    labels = ['', "BCEP, BCEPS", '', 'DSCT', 'SR', "SRA, SRB, SRC, SRD", '', '', '', 'RR', "RRAB, RRC", '', 'GDOR',
              'SPB', '', 'LPV']
    for i in range(len(variabletype)):
        if i in [2, 6, 7, 8, 11]:
            my_label = None
        else:
            my_label = "%s" % labels[i]
        plt.scatter(variablesdf[x].loc[variablesdf.loc[:, 'Type'] == variabletype[i]], variablesdf[y]
                    .loc[variablesdf.loc[:, 'Type'] == variabletype[i]], facecolor=colors[i], marker=markers[i],
                    s=sizes[i], label=my_label)
        print "plotting %s as %s%s" % (variabletype[i], colors[i], markers[i])
    return


def get_variable_stars(df_data, df_variables_names, variabletype=None):
    """ Child function fo plot_variable_stars. Process the DataFrame to select only stars marked as
    'var_type' variable stars.

    :param df_data: pandas.DataFrame with the sample stars (background and variables)
    :param df_variables_names: pandas.DataFrame containing the ID of variables stars
    :param variabletype: types of variables of interest - leave it to None
    :return: a subset of the df_data, containing only the variable stars
    """
    if variabletype is None:
        variabletype = ['CEP', 'BCEP', 'BCEPS', 'DSCT', 'SR', 'SRA', 'SRB', 'SRC', 'SRD', 'RR', 'RRAB', 'RRC',
                        'GDOR', 'SPB', 'M', 'LPV']

    print "Selecting variable stars.."
    are_variables = df_variables_names[df_variables_names.loc[:, 'Type'].isin(variabletype)]
    typesdf = are_variables[['hip', 'tycho2_id', 'source_id', 'Type', 'Name']]
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


def get_hip_sources(df_main_catalog=None, variable_class=None):
    if variable_class is None:
        variable_class = ['CEP', 'BCEP', 'BCEPS', 'DSCT', 'SR', 'SRA', 'SRB', 'SRC', 'SRD', 'RR', 'RRAB', 'RRC',
                          'GDOR', 'SPB', 'M', 'LPV']
    hip_sources = pd.read_csv('Hip/variables.tsv', delimiter=';', header=0)
    variable_stars = hip_sources[hip_sources.loc[:, 'VarType'].isin(variable_class)]
    variable_stars.hip = variable_stars.hip.astype(float)

    hip_objects = df_main_catalog[df_main_catalog.hip.isin(variable_stars.hip)]
    hip_objects = pd.merge(hip_objects, variable_stars, on='hip', how='inner')
    hip_objects = hip_objects.rename(columns={'VarType': 'Type'})

    return hip_objects


def remove_misclassified_objects(data_frame):
    """
    -- added a Variability flag check in data_process, where only CONFIRMED variables are now kept --

    drop misclassified objects, like SR stars on the Main Sequence, optical binaries, etc.
    :param data_frame: pandas.DataFrame containing all the variable stars
    :return: cleaned up DataFrame, with misclassified objects removed
    """

    misclassified_objects = ['7720-1455-1', '7911-499-1',  # Mira
                             '8990-3504-1', '899-471-1',  # SRs with very poor light curves
                             '6430-88-1',  # SV For, obsolete measurement is p=91d, newer is p=16h
                             '9017-396-1', '7676-2953-1', '8296-3303-1',  # SR
                             '2365-2764-1', '4109-638-1', '2058-56-1',  # Cepheids
                             '3642-2459-1', '3999-1391-1', '2607-1448-1',  # Cepheids
                             '3655-469-1', '1476-148-1', '1233-531-1',  # RR Lyrae
                             '3029-738-1', '6863-1255-1', '6954-1236-1', '9380-420-1'  # RR Lyrae
                             '6192-461-1',  # DSCT
                             '2550-686-1', '4992-357-1', '9380-420-1', '8562-728-1', '6567-2007-1', '6040-2003-1',  # RR
                             '2553-1108-1',
                             '4851-2441-1', '8962-577-1',
                             '3135-132-1', '8976-3674-1', '3136-437-1', '1506-618-1', '7046-1715-1', '3140-3046-1',
                             '2000-162-1', '6210-755-1', '3547-1807-1', '8836-935-1', '3033-273-1', '7606-437-1',
                             '3049-180-1', '9198-1862-1', '8192-626-1', '7703-1577-1', '8594-433-1', '6833-280-1'
                             ]
    dsct = data_frame[(data_frame.Type == 'DSCT') & (data_frame.B_V > 0.4) & (data_frame.M_V > 2.5)].tycho2_id.tolist()
    dsct2 = data_frame[(data_frame.Type == 'DSCT') & (data_frame.B_V > 0.25) & (data_frame.M_V > 3.)].tycho2_id.tolist()
    print "Dropping objects DSCT: %s %s" % (dsct, dsct2)
    data_frame = data_frame.drop(data_frame[data_frame.tycho2_id.isin(dsct)].index)
    data_frame = data_frame.drop(data_frame[data_frame.tycho2_id.isin(dsct2)].index)
    print "Dropping objects: %s" % misclassified_objects
    data_frame = data_frame.drop(data_frame[data_frame.tycho2_id.isin(misclassified_objects)].index)
    print "..Done\n----------"
    return data_frame


def deredden_cepheids(df_variables):
    """
    compute the dereddened values of B-V and M_V for the six Cepheids in our sample (parallax cutoff = 0.25 and
    B-V error < 0.1).
    :param df_variables: DataFrame of the variable stars
    :type df_variables: pandas.core.frame.DataFrame
    :return: updated DataFrame of variable stars
    """
    extinction_coefficients = {'2365-2764-1': np.array([0.2622, 0.844]), '4109-638-1': np.array([0.0524, 0.1576]),
                               '2058-56-1': np.array([0.0751, 0.248]), '3642-2459-1': np.array([0.1907, 0.608]),
                               '3999-1391-1': np.array([0.3911, 1.2480]), '2607-1448-1': np.array([0.0430, 0.1310])}
    print "Dereddening Cepheids:"
    for tyc in extinction_coefficients.keys():
        print "%s.." % tyc
        b_minus_v = df_variables[df_variables.tycho2_id == tyc].B_V
        m_v = df_variables[df_variables.tycho2_id == tyc].M_V
        extinc = extinction_coefficients[tyc]
        df_variables.set_value(df_variables.tycho2_id == tyc, 'B_V', b_minus_v - extinc[0])
        df_variables.set_value(df_variables.tycho2_id == tyc, 'M_V', m_v - extinc[1])
    print "..Done\n----------"

    return df_variables


def plot_dereddening():
    """
    plot a Cepheid at its reddened position on the HR diag. (assume that deredden_cepheids() have been used)
    Reddening coefficients are taken from http://irsa.ipac.caltech.edu/applications/DUST/ using object's RA/DEC
    :return:
    """
    extinction_coefficients = {'2365-2764-1': np.array([0.2622, 0.844]), '4109-638-1': np.array([0.0524, 0.1576]),
                               '2058-56-1': np.array([0.0751, 0.248]), '3642-2459-1': np.array([0.1907, 0.608]),
                               '3999-1391-1': np.array([0.3911, 1.2480]), '2607-1448-1': np.array([0.0430, 0.1310])}
    cepheids = {'2365-2764-1': np.array([0.959, 2.09]), '4109-638-1': np.array([0.705, 2.385]), '2058-56-1':
                np.array([1.222, 1.333]), '3642-2459-1': np.array([1.088, 2.0518]), '3999-1391-1':
                np.array([1.360, 1.2567]), '2607-1448-1': np.array([1.484, 0.6963])}
    periods = {'2365-2764-1': 1.61, '4109-638-1': 15.31, '2058-56-1': 63.08, '3642-2459-1': 1.86, '3999-1391-1': 24.98,
               '2607-1448-1': 8.54}
    max_periods = max(periods.values())

    new_positions_bv_mv = []  # in M_V vs B-V space
    colors = []
    theoretical_position = []
    for obj in extinction_coefficients.keys():
        # new_positions_bv_mv.append(cepheids[obj]-extinction_coefficients[obj])
        new_positions_bv_mv.append(cepheids[obj])
        colors.append(periods[obj]/max_periods)
        theoretical_position.append(-2.78*np.log10(periods[obj])-1.35)

    for pos in range(len(new_positions_bv_mv)):
        plt.scatter(new_positions_bv_mv[pos][0], new_positions_bv_mv[pos][1], marker='^', facecolor='w', s=40)
        plt.scatter(new_positions_bv_mv[pos][0], theoretical_position[pos], marker='o', facecolor='r', s=50)
    return new_positions_bv_mv, colors

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

    df = data_process(df, cutoff=args.cutoff, bv_cutoff=args.bvcutoff)

    df2 = import_data(catalog='xmatch_TGAS_VSX.csv', params=('hip', 'tycho2_id', 'source_id', 'parallax',
                                                             'parallax_error', 'Name', 'V', 'Type'), nrows=args.nrows)
    df2 = data_process(df2, catalog='xmatch_TGAS_VSX.csv', cutoff=args.cutoff)

    variable_types = ['CEP', 'BCEP', 'BCEPS', 'DSCT', 'SR', 'SRA', 'SRB', 'SRC', 'SRD', 'RR', 'RRAB', 'RRC',
                      'GDOR', 'SPB', 'M', 'LPV']

    list_variables = get_variable_stars(df, df2, variable_types)

    list_variables = remove_misclassified_objects(list_variables)
    list_variables = deredden_cepheids(list_variables)
    print "Saving 'list_variables.csv' on disk.."
    list_variables.to_csv("list_variables.csv")
    print "..Done\n----------"

    ########################################################
    # Plotting of HR diag and variables stars (with errors): #
    ########################################################
    if args.draw is True:
        plot_full(df, list_variables, variable_types, x='B_V', y='M_V', cutoff=args.cutoff, bvcutoff=args.bvcutoff)
        plt.legend(loc="lower right", fontsize=15, scatterpoints=1)
        fig = plt.figure(num=1)
        fig.set_figheight(10.)
        fig.set_figwidth(11.)
