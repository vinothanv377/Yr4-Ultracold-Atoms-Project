# -*- coding: utf-8 -*-
# pylint: disable=E0012, fixme, invalid-name, no-member, R0913

''' Assignment 3 (working with files; compulsory)

Assignment Tasks: 3

Aim:
    Read experimental data from files and plot data averaged over
    the multiple experimental runs.

Restrictions:
    Do not change anything outside TODO blocks.
    Do not use import, scatter
    Do not add pylint directives.

Guidance:
    Data files are compressed in a single zip file. All data files have
    identical structure (number of lines, etc.).
    Download the zip file, extract the data files into a sub-folder in
    your u:/Python folder.
    Load one of the files into the Spyder editor to understand its structure.

    Task 1:
        Define identifiers used as constants to be used throughout the script,
        including inside functions.
        - Set the data folder name. Ensure cross-platform compatibility
        - Set the list of data column names to reflect the names given in the
          data file you inspected. (You could in principle read these out
          automatically from the file, but it it would serve no purpose to do
          so in this case.)
        - Set number of lines to skip at start of file before reading data.

    Task 2:
        - Ensure that your function will work cross-platform.
        - Use np.loadtxt() with suitable named parameters.
        - Read the first file to determine the data shape, then create an empty
          3d array of suitable size and splice each file's 2d data into it.

    Task 3:
        Generate a fully labelled graph. Use plot(), not scatter().

Author of template:
    Wolfgang Theis
    School of Physics and Astronomy
    University of Birmingham
'''

import os
import numpy as np
import matplotlib.pyplot as plt

# TODO: Assignment Task 1: define the following constants

# path name of folder in which the data files are located
FOLDER = r'//adf/storage/V/V/VXV762/Python/FH 190 Degrees Data'
# list of column names (as in the file)
COL_NAMES = ['Time (s)', 'Channel A (V)', 'Channel B (V)']
# lines to skip when reading data files
LINES_TO_SKIP = 1
# End of Task 1; proceed to Task 2


def read_all_csvfiles_in_folder(folder_name):
    """ Read all CSV files in a folder


    Parameters
    ----------

    folder_name: string
        folder name including path

    Returns
    -------
    data: array of floats
        data read from all the files in the folder (3d array)
        data[:,:,k] is the data from file k
    """

    # TODO: Assignment Task 2: write function body
    list_of_filenames = os.listdir(folder_name)
    for k, fn in enumerate(list_of_filenames):#loops through all the files
        full_path = folder_name + os.sep + fn
        filedata = np.loadtxt(full_path, skiprows=LINES_TO_SKIP, delimiter=',')
        #loads the text
        if k == 0:
            fshape = np.shape(filedata)#shape of the file
            data = np.empty((fshape[0], fshape[1], len(list_of_filenames)), dtype=float)
            #shape of the data
        data[:, :, k] = filedata #format of the file data
    return data





    # End of Task 2;  proceed to Task 3


def generate_plot(ax, data):
    """ Create a fully labelled plot

    Plots data (Channel B versus Channel A) from the first file as
    small blue dots.
    Plots the averaged data (averaged over the set of repeat experiments,
    i.e. the files) as small red dots.
    Provides a legend.
    Use the COL_NAMES defined at the top for the axis labels.

    Uses y=0 as the lower limit for the y axis.
    Restricts the x axis to the data range.
    Shows 'Franck-Hertz Experiment' as the title.
    Provides an annotation text with your name and date (the date is that of
    writing the program, not date of running the program).

    Parameters
    ----------

    ax: matplotlib axes
        object to draw the graph into

    data: array of floats
        3d array [row, column, number of file]

    Returns
    -------
    ax: matplotlib axes
        fully labelled plot
    """
    # TODO: Assignment Task 3: write function body
    ax.plot(data[:, 1, 0], data[:, 2, 0], 'b:', label='Data File 1') #plot of file 1
    avg_data = np.mean(data, axis=2) #average of the data per input along the file list
    ax.plot(avg_data[:, 1], avg_data[:, 2], 'ro', label='Average Data')
    #plot of the averaged data across files
    ax.set_xlabel(COL_NAMES[1])
    ax.set_ylabel(COL_NAMES[2])
    ax.set_ylim((0))
    ax.legend(loc='best')
    ax.set_title("Franck-Hertz Experiment")
    ax.set_xlim(0, np.max(data[:, :, :]))
    ax.text(1, 0.010, 'Vinothan Vaheesan\n26.11.2016', transform=ax.transAxes,
            horizontalalignment='right', fontsize=14)
    return ax




    # End of Task 3;  no further tasks.


def load_and_plot():
    ''' Load the data and plot it'''
    data = read_all_csvfiles_in_folder(FOLDER)
    fig = plt.figure(num=1, figsize=(16, 12))
    ax = fig.add_subplot(1, 1, 1)
    generate_plot(ax, data)
    plt.show()



if __name__ == '__main__':
    # pylint: disable=E0012, wrong-import-position, wrong-import-order
    # Feel free to modify this flag to switch between plotting and unit tests
    plot_it = True

    if plot_it:
        # plot a test figure
        load_and_plot()
    else:
        # do the unit tests
        import unittest

        class load_Tests(unittest.TestCase):
            ''' tests read_all_csvfiles_in_folder() '''

            def test_1(self):
                ''' tests correct array shape '''
                self.assertEqual(read_all_csvfiles_in_folder(FOLDER).shape, (1024, 3, 100))

            def test_2(self):
                ''' tests data not identical '''
                d = read_all_csvfiles_in_folder(FOLDER)
                self.assertFalse(np.all(d[:, :, 0] == d[:, :, 1]))


        unittest.main(exit=False)
        
