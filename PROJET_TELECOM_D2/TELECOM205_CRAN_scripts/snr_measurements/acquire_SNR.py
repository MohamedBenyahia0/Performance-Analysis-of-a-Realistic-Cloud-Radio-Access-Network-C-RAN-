# This python script is used to acquire the SNR of the PlutoSDR
# It is based on the example from the Analog Devices examples
# https://raw.githubusercontent.com/analogdevicesinc/pyadi-iio/master/examples/pluto.py


import time

import adi
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from scipy import fft

import utilities

# Ask user the distance for file exports
distance = input("Please enter the distance in meters: ")

# RF parameters
Fc                = 0.6e9;   # Central frequency of the band
BW                = 20e6;    # Considered bandwidth of the system
Fsig_offset       = 1e6;     # Offset of the signal from the center frequency

# Create radio
sdr = adi.Pluto("ip:pluto.local")

# Number of acquisitions
N_acq = 15
# Samples per capture
samples_per_frame = 2**13

# Configure properties
sdr.rx_rf_bandwidth           = int(BW)
sdr.sample_rate               = int(30.72e6)
sdr.rx_lo                     = int(Fc)
sdr.rx_enabled_channels       = [0]
sdr.gain_control_mode_chan0   = "manual"
sdr.rx_buffer_size            = samples_per_frame

# Read properties
print("RX LO %s" % (sdr.rx_lo))

# Define a vector of gains to try
gain_start    = -20.
gain_stop     = 40.
gain_step     = 5.
gain_vec      = np.arange(gain_start, gain_stop, gain_step)

# Define a matrix to store the captured data
data_mat = np.zeros((samples_per_frame, len(gain_vec), N_acq), dtype=np.complex64)

# Collect data
for n_gain in range(len(gain_vec)):
    print("Acquiring data for gain = "+str(gain_vec[n_gain])+"dB")
    print("["+str(n_gain+1)+"/"+str(len(gain_vec))+"]")
    sdr.rx_hardwaregain_chan0   = gain_vec[n_gain]
    for n_acq in range(N_acq):
        print('.', end='', flush=True)
        data_mat[:,n_gain,n_acq] = sdr.rx()
        time.sleep(0.2)

# Compute the PSD for each capture
print("Compute the PSD for each capture")
data_psd_mat = np.zeros((samples_per_frame, len(gain_vec), N_acq), dtype=np.float64)
for n_gain in range(len(gain_vec)):
    for n_acq in range(N_acq):
        f, data_psd_mat[:,n_gain,n_acq] = signal.periodogram(data_mat[:,n_gain,n_acq],BW,window='blackman',return_onesided=False)

# Compute the average PSD for each gain value
print("Compute the average PSD for each gain value")
data_psd_mat_avg = np.mean(data_psd_mat, axis=2)

# Compute sig_power
sig_power_vec = utilities.sig_power(data_psd_mat_avg, Fsig_offset, sdr.sample_rate)
sig_power_vec_dB = 10*np.log10(sig_power_vec)
# print(str(10*np.log10(sig_power)))

#################################################
# Plotting section and file exports
#################################################

# # Debug
# plt.clf()
# plt.semilogy(f, data_psd_mat_avg[:,0])
# # plt.ylim([1e-7, 1e2])
# plt.xlabel("frequency [Hz]")
# plt.ylabel("PSD [V**2/Hz]")
# plt.draw()
# plt.show()

plot_filename = "distance_"+str(distance)+"m_spectrums.pdf"
utilities.plot_spectrums(gain_vec, data_psd_mat_avg, f, plot_filename)

# Repack data for CSV export
data_csv=np.column_stack((gain_vec, sig_power_vec_dB))

# Export signal power to a CSV file
csv_filename = "distance_"+str(distance)+"m_sig_pows.csv"
# https://stackoverflow.com/questions/36210977/python-numpy-savetxt-header-has-extra-character
np.savetxt(csv_filename, data_csv, delimiter=",", header="Gain [dB], Signal Power [dB]", comments="")

# Export average PSDs to a CSV file
csv_filename = "distance_"+str(distance)+"m_psds.csv"
np.savetxt(csv_filename, data_psd_mat_avg, delimiter=",", header="PSD [V**2/Hz]", comments="")
