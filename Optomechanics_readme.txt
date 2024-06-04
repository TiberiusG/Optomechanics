C4--XX--00007.trc is a binary file saved by our LeCroy HDO6104 oscilloscope. The file was recorded when the air around the optically levitated particle had a pressure of 6.7e0 mbar, in a lab whose temperature was 20 degrees Celsius.

hpsd.m is a function that computes the Power Spectral Density of the time trace.
The input parameters are:
A == the amplitude signal of the oscilloscope, or the vertical component of the scope file;
nos == number of segments to which the time trace is divided and then averaged in the PSD estimate algorithm.
the default is nos = 10, which means that a 1 s time trace, for example, is divided into 10 x 100 ms traces, then the PSD is calculated for each of the 10 x 100 ms traces and then averaged.
fs == the sampling frequency used to acquire the data.

calculate_and_save_Hpsd_for_all_files.m
calculates and saves the PSDs for all the .trc files in the folder containing the sequence written on line 16.
One needs the .trc scope files and the ReadLeCroyBinaryWaveform.m program to be in the same folder as calculate_and_save_Hpsd_for_all_files.m.
 
open_trc_hPSD_Hfit_write_abc_peaks.m
The values of the parameters:
f_min, f_max, f0, nos and k must be introduced manually.

This program does the following:
-opens the scope file that contains the reference data, in this case 'C4--XX--00004.trc';
-computes the PSD of the reference time trace;
-displays the PSD between f_min = 102e3 and f_max = 107e3 Hz;
-asigns to this peak the value 3, k = 3 (particle translation perpendicular to the electric field).
This value will then be written in the name of the saved file, to distinguish the particle's translational peaks:
1 == axial, lowest frequency;
2 == parallel to the electric field of the laser in the focal spot;
3 == perpendicular to the electric field of the laser in the focal spot.

[
The values of the parameters:
f_min, f_max, f0, nos and k must be introduced manually.
]

-fits the function:
fit = A1/pi * 2*f0^2*g ./...
    ( (x.^2 - f0^2).^2 + x.^2*g^2 ) + n;
to the experimental PSD between f_min and f_max;

-outputs the found fit parameters:
A1 
f0
g
n
in a .txt document whose name contains information about the scope channel, time trace number and peak number, among other things.


in_A1f0g_out_rmg_210407.m
reads the .txt file containing the fit parameters extracted by open_trc_hPSD_Hfit_write_abc_peaks.m
and computes the particle paramters such as radius, conversion factor from V to m, damping factor and others.
It then saves these extracted parameter into a .txt file.



ReadLeCroyBinaryWaveform.m interprets the .trc file saved by the LeCroy oscilloscope. This file has been written by:
 Original version (c)2001 Hochschule f?r Technik+Architektur Luzern
 Fachstelle Elektronik
 6048 Horw, Switzerland
 Slightly modified by Alan Blankman, LeCroy Corporation, 2006
 July 2007; added decoding of TrigTime_array so that sequence mode waveforms decode correctly
 Slightly modified by Jean-Daniel Desch?nes, Universite Laval, 2010

More information can be found within the file itself.