# Authors: Nathanyel Schut and Dirk Kuiper
# Studentnumbers: 12907995 and 12416657
# Date: February 2021
# Script for analyzing Single Dot Spectroscopy measurement data

import pandas as pd
import numpy as np
from statistics import mean, stdev
import matplotlib.pyplot as plt
from lmfit import models
from PyQt5.QtCore import (Qt, pyqtSignal)

# Variables to adjust threshold to filter out solar flares
raw_data_threshold = 2000                       # Raw data
raw_data_no_background_threshold = 700          # Raw data but with background reduction
corrected_data_threshold = 10                   # Corrected data
corrected_data_no_background_threshold = 3      # Corrected data with background reduction

# All the constants we're using
v_light = 299792458
planck_const = 6.62607004 * 10**(-34)
J = 6.24150913 * 10**(18)   # eV

# All the factors used in the correction curve formula, 
# these need to be changed if you have a different correction formula for your lamp
a = 41.6428660875473
b = -4914.66877226876
c = 0.71262780761720
d = 644.241545917904
e = -437088.417274657
f = 100016023.454422
g = -6757961430.13629
h = 0

# Finding the expected intensity of the thungsten lamp
def correction_factors(meas_corr_curve_file_path, CCD_height, CCD_width, lower_wavelength, upper_wavelength):
    """Function for calculating the correction factors with which to multiply the measured intensities of 
    the sample.
    
    Arguments:
        meas_corr_curve_file_path (str or path): path for the file containing the measured intensities of the thungsten lamp.
        CCD_height (int): number of pixels along the vertical axis of the CCD.
        CCD_width (int): number of pixels along the vertical axis of the CCD.
        lower_wavelength (float): the lowest wavelength of the spectrum in nanometers.
        upper_wavelength (float): the highest wavelength of the spectrum in nanometers.
    
    Returns:
        A list of correction factors sorted by corresponding wavelength.
    """

    # Constructs a filepath and a dataframe for the correction curve data
    df = pd.read_csv(meas_corr_curve_file_path)
    df1 = pd.DataFrame({"1": df.Intensity.loc[0:CCD_width - 1]})
    for i in range(2, CCD_height + 1):
        df1[f'{i}'] = df.Intensity.loc[((i - 1)*CCD_width):(i*CCD_width - 1)].reset_index(drop=True)
    
    # Takes the mean intensity value for every wavelength in the correction curve and appends it to a list
    samplewavelength = df.Wavelength.loc[0:CCD_width - 1]
    df1['Wavelength'] = samplewavelength
    df1 = df1.set_index('Wavelength')
    df1['mean'] = df1.mean(axis=1)
    meas_lamp_radiance = df1['mean'].to_numpy()

    # Loops through the mean values and constructs the radiance via the given formula
    lampradiance = []
    for i in samplewavelength:
        # This is the formula for the correction factor, this formula needs to be changed along with the constants
        # if you use a different correction curve
        lampvalue = (i**(-5))*np.exp(a+b/i)*(c+d/i+e/(i**2)+f/(i**3)+g/(i**4)+h/(i**5))
        lampradiance.append(lampvalue)
    
    # Constructs a list with correction factors for every wavelength from the calculated radiances vs the measured radicances
    correctionfactors = []
    for i in range (0,len(samplewavelength)):
        correctionfactor = lampradiance[i]/meas_lamp_radiance[i]
        correctionfactors.append(correctionfactor)
    
    return correctionfactors

class QDot_Spectroscopy():
    """A class to analyze single dot spectroscopy data.
    
        
    """
    def __init__(self, file_path, meas_corr_curve_file_path, lower_wavelength, upper_wavelength, CCD_height, CCD_width, CCD_height_corr, CCD_width_corr):
        """Method for creating the dataframes which hold the data for the 2d maps. 
        
        From these dataframes, more information can be extracted such as the location of
        the quantum dots.
        
        Arguments:
            file_path (string or path): path of the datafile to be analyzed.
            meas_corr_curve_file_path (string or path): path of the correction file for correcting the measurement data.
            lower_wavelength (float): the lowest wavelength of the spectrum in nanometers.
            upper_wavelength (float): the highest wavelength of the spectrum in nanometers.
            CCD_height (int): number of pixels along the vertical axis of the CCD.
            CCD_width (int): number of pixels along the vertical axis of the CCD.
            CCD_height_corr (int): number of pixels along the vertical axis of the CCD for the correction data.
            CCD_width_corr (int): number of pixels along the vertical axis of the CCD for the correction data.
        """
        
        self.file_path = file_path
        self.df = pd.read_csv(self.file_path)
        self.CCD_height = CCD_height
        self.CCD_width = CCD_width
        # self.fig = plt.figure(num=1)
        self.correctionfactors = correction_factors(meas_corr_curve_file_path, CCD_height=CCD_height_corr , CCD_width=CCD_width_corr, lower_wavelength=lower_wavelength, upper_wavelength=upper_wavelength)

        # Determine wavelength and energy values
        self.wavelengths = np.linspace(lower_wavelength, upper_wavelength, self.CCD_width)
        self.energies = (J * planck_const * v_light) / (self.wavelengths * 10**(-9))

        # Bin correctionfactors
        self.correctionfactors_res = []
        for i in range(1, self.CCD_width + 1):
            slice_begin = (i - 1) * int((len(self.correctionfactors) / self.CCD_width))
            slice_end = i * int((len(self.correctionfactors) / self.CCD_width))

            binned_data = self.correctionfactors[slice_begin:slice_end]

            self.correctionfactors_res.append(mean(binned_data))

        # self.correctionfactors_res = pd.DataFrame({"Correction Factors":self.correctionfactors_res})

        # Construct matrix from data. Columns represent different measurements, rows represent measurement data
        #        |   Raw   |  Corr   |
        #--------+---------+---------+
        # Raw    |   df1   |   df3   |
        # No bkg |   df2   |   df4   |

        self.correctionfactors_res.reverse()
        self.df1 = pd.DataFrame({"Correction Factors": self.correctionfactors_res})
        self.df3 = pd.DataFrame({"Correction Factors": self.correctionfactors_res})
        for i in range(1, self.CCD_height + 1):
            self.df1[f'{i}'] = self.df.Intensity.loc[((i - 1)*self.CCD_width):(i*self.CCD_width - 1)].reset_index(drop=True)
            self.df3[f'{i}'] = self.df.Intensity.loc[((i - 1)*self.CCD_width):(i*self.CCD_width - 1)].reset_index(drop=True)
            self.df3[f'{i}'] = self.df3[f'{i}'] * self.df1["Correction Factors"]
        
        self.df1 = self.df1.drop('Correction Factors', axis=1)
        self.df3 = self.df3.drop('Correction Factors', axis=1)

        # Construct differential matrix
        self.df2 = pd.DataFrame({"1": self.df1['3'] - self.df1['2']})
        self.df4 = pd.DataFrame({"1": self.df3['3'] - self.df3['2']})
        
        for i in range(4, self.CCD_height + 1):
            self.df2[f'{i - 2}'] = self.df1[f'{i}'] - self.df1[f'{i - 1}']
            self.df4[f'{i - 2}'] = self.df3[f'{i}'] - self.df3[f'{i - 1}']   

        # Change of index
        self.df1['index'] = self.energies
        self.df2['index'] = self.energies
        self.df3['index'] = self.energies
        self.df4['index'] = self.energies
        self.df1 = self.df1.set_index('index')
        self.df2 = self.df2.set_index('index')
        self.df3 = self.df3.set_index('index')
        self.df4 = self.df4.set_index('index')

        # Sets all the values that are negative to 0 and filters out the solar flares by setting an upper limit for the intensity values
        self.df2[self.df2 < 0] = 0
        self.df4[self.df4 < 0] = 0
    
        self.df1[self.df1 > raw_data_threshold] = 0
        self.df2[self.df2 > raw_data_no_background_threshold] = 0
        self.df3[self.df3 > corrected_data_threshold] = 0
        self.df4[self.df4 > corrected_data_no_background_threshold] = 0

    def matrix_map(self, bkg_reduction=True, data_correction=True):
        """A method for conveniently selecting the required data for the 2d maps.

        This methods allows the user to choose whether or not to return corrected or raw data.

        Arguments:
            bkg_reduction (Bool, optional, default=True): Boolean for enabling or disabling background reduction.
            data_correction (Bool, optional, default=True): Boolean for enabling or disabling correction of the data.

        Returns:
            A data frame containing either raw or corrected data depending on the users choice.
        """

        if bkg_reduction is True:
            if data_correction is True:
                data = self.df4
            
            else:
                data = self.df2

        else:
            if data_correction is True:
                data = self.df3
            
            else:
                data = self.df1

        return data
        
    def QDot_detection(self):
        """Method for detecting the lines containing a posiible quantum dot.

        This method sums up the intensities for each lines and compares the total intensity
        to the 3 sigma threshold of the total intensities. If the total intensity of a line is
        higher than this threshold, the line will be regarded as containing a quantum dot.

        Returns:
            A list of line numbers containing a quantum dot.
        """

        # Creates a list with the total intensities from the lines as to analyze which lines contain quantum dots
        total_intensity_list = []
        for (columnName, columnData) in self.df4.iteritems():
            total_intensity = 0
            for i in columnData.values:
                total_intensity += i
            total_intensity_list.append(total_intensity)
        
        # Construct the 3-sigma threshold
        avg_tot_intensity = mean(total_intensity_list)
        stdev_tot_intensity = stdev(total_intensity_list)

        threshold = 3 * stdev_tot_intensity + avg_tot_intensity

        QDot_slits = [total_intensity_list.index(i) + 1 for i in total_intensity_list if i >= threshold]
        
        # If 2 lines next to each other are labeled as quantum dots, the slit with the lowest total intensity will be discarded
        to_be_deleted_slits = []
        for i in range(0, len(QDot_slits) - 1):
            if QDot_slits[i + 1] - QDot_slits[i] == 1:
                if total_intensity_list[QDot_slits[i + 1] - 1] > total_intensity_list[QDot_slits[i] - 1]:
                    to_be_deleted_slits.append(QDot_slits[i])
                elif total_intensity_list[QDot_slits[i + 1] - 1] < total_intensity_list[QDot_slits[i] - 1]:
                    to_be_deleted_slits.append(QDot_slits[i + 1])
        
        for slit in to_be_deleted_slits:
            QDot_slits.remove(slit)

        # Optional code to plot the total intensities of every slit in the 2d map.
        # -------------------------------------------------------------------------
        # fig = plt.figure(figsize=(10,7))
        # plt.plot(total_intensity_list, label="Total intensity")
        # plt.plot([x - 1 for x in QDot_slits], [total_intensity_list[x - 1] for x in QDot_slits], 'rx', label="SI-NP")
        # plt.hlines(avg_tot_intensity, 0, 200, colors='red', label="Average total intensity")
        # plt.hlines(threshold, 0, 200, colors='green', label='3-sigma threshold')
        # plt.title("Total intensities for a single datafile")
        # plt.xlabel("Position along the slit (pixels)")
        # plt.ylabel("Total intensity (arbitrary units)")
        # plt.xlim(0,200)
        # plt.ylim(0,60)
        # plt.legend()
        # plt.show()

        return QDot_slits

    def gaussian_fit(self):
        """Method for fitting a gaussian curve to the spectrum of possible quantum dots.
        
        Returns:
            A data frame containing the spectrum of a quantum dot as well as the gaussian curve plot data for that quantum dot.
            This method also returns a data frame with the fit statistics of all the fits.
        """

        self.df5 = pd.DataFrame(columns=['Slit Number', 'Centre', 'Centre_err', 'Sigma', 'Sigma_err', 'FWHM', 'FWHM_err', 'Height', 'Height_err'])
        QDot_slits = self.QDot_detection()

        if len(QDot_slits) > 0:    
            self.plot_data = pd.DataFrame(columns=[f"{QDot_slits[0]}"], index=self.energies)
        else:
            self.plot_data = pd.DataFrame(index=self.energies)

        for slit_number in QDot_slits:
            sel = self.df4[f'{slit_number}']
            self.plot_data[f'{slit_number}'] = sel
            
            # Makes a good first guess for the fit values of the gaussian
            max_intensity = max(sel)
            central_energy = sel[sel==max_intensity].index.values
            central_energy = central_energy[0]

            # Fits a gaussian model to the selected data and shows the output
            gauss = models.GaussianModel()
            fit = gauss.fit(sel, x=self.energies, weights=1 / np.sqrt(sel), center = central_energy, amplitude = max_intensity, sigma = 1, nan_policy= 'omit')
            
            self.plot_data[f'{slit_number} best fit'] = fit.best_fit

            # Appends the fit data for the variables to a new dataframe and shows the fit results with errors
            fit_variables = [slit_number]
            for key in fit.params:
                if key in ['center', 'sigma', 'fwhm', 'height']:
                    fit_variables.append(fit.params[key].value)
                    fit_variables.append(fit.params[key].stderr)
            
            self.df5 = self.df5.append({'Slit Number': fit_variables[0], 'Centre': fit_variables[1], 'Centre_err': fit_variables[2], 'Sigma': fit_variables[3], 'Sigma_err': fit_variables[4], 'FWHM': fit_variables[5], 'FWHM_err': fit_variables[6], 'Height': fit_variables[7], 'Height_err': fit_variables[8]}, ignore_index=True)
        
        return self.plot_data, self.df5

def multiple_files_analysis(lower_wavelength, upper_wavelength, CCD_height, CCD_width, CCD_height_corr, CCD_width_corr, file_paths, file_path_corr_data, progress_update):
    """Method for analyzing multiple files consecutively.

    This method takes a list of paths for data files and runs the analysis program for each data file.
    Additionally this method extracts the fit data for all the quantum dots in each file and creates
    lists for the FWHM values and the central energy values of all the fits.

    Arguments:
        lower_wavelength (float): the lowest wavelength of the spectrum in nanometers.
        upper_wavelength (float): the highest wavelength of the spectrum in nanometers.
        CCD_height (int): number of pixels along the vertical axis of the CCD.
        CCD_width (int): number of pixels along the vertical axis of the CCD.
        CCD_height_corr (int): number of pixels along the vertical axis of the CCD for the correction data.
        CCD_width_corr (int): number of pixels along the vertical axis of the CCD for the correction data.
        file_paths (list of strings or paths): paths of the datafiles to be analyzed.
        file_path_corr_data (string or path): path of the correction file for correcting the measurement data.
        progress_update (pyqtSignal): a signal to update.
    
    Returns:
        A 2d array containing the data for all the 2d maps, the plot data of the quantum dots and the fit statistics
        of the quantum dots for every file. The program also returns the FWHM values and the central energy values.  
    """
    
    all_files_data = []
    FWHM_data = []
    central_energy_data = []
    counter = 1

    for file_path in file_paths:
        analysis = QDot_Spectroscopy(file_path=r"{}".format(file_path), meas_corr_curve_file_path=r"{}".format(file_path_corr_data), lower_wavelength=lower_wavelength, upper_wavelength=upper_wavelength, CCD_height=CCD_height, CCD_width=CCD_width, CCD_height_corr=CCD_height_corr , CCD_width_corr=CCD_width_corr)

        twod_map_raw = analysis.matrix_map(bkg_reduction=False, data_correction=False)
        twod_map_no_bkg = analysis.matrix_map(bkg_reduction=True, data_correction=False)
        twod_map_raw_corr = analysis.matrix_map(bkg_reduction=False, data_correction=True)
        twod_map_no_bkg_corr = analysis.matrix_map(bkg_reduction=True, data_correction=True)
        Q_Dot_plot_data, fit_statistics = analysis.gaussian_fit()

        file_analysis = [twod_map_raw, twod_map_no_bkg, twod_map_raw_corr, twod_map_no_bkg_corr, Q_Dot_plot_data, fit_statistics]
        all_files_data.append(file_analysis)

        # Creates a histogram from the collected FWHM and central energy data from all the analyzed datafales containing quantumdots
        for FWHM_value in fit_statistics['FWHM'].to_numpy():
            FWHM_data.append(FWHM_value)
        for CE_value in fit_statistics['Centre'].to_numpy():
            central_energy_data.append(CE_value)

        progress_update.emit(counter * 100/len(file_paths))
        counter += 1
    
    return all_files_data, FWHM_data, central_energy_data

def histogram(hist_data, bins, minimum, maximum):
    """Method for fitting a Gaussian curve to the histogram data.
    
    Arguments:
        hist_data (array like): list or array containing the data for which to make a histogram.
        bins (int): ammount of bins for the histogram.
        minimum (float): the minimum value to contain in the histogram.
        maximum (float): the maximum value to contain in the histogram.

    Returns:
        The data for the histogram plot as well as the Gaussian curve plot.
    """
    hist_data2 = [x for x in hist_data if x <= maximum and x >= minimum]
    y, x = np.histogram(hist_data2, bins=bins)
    x_new = [((x[i] + x[i - 1]) / 2) for i in range(1, len(x))]
    max_freq = max(y)
    mean = x_new[y.argmax()]
    sigma = (x[1] - x[0]) * 4

    gauss = models.GaussianModel()
    fit = gauss.fit(y, x=x_new, center = mean, amplitude = max_freq, sigma = sigma, nan_policy= 'omit')

    return x, y, x_new, fit.best_fit