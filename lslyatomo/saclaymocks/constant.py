# physical constant and parameters
import scipy as sp
from scipy import interpolate
import scipy.constants

c = scipy.constants.speed_of_light / 1000.  # in km/s

n_qso_exp = 100.  # expected number of qso per square degrees
qso_nz_adhoc = 0.213  # adhoc factor to have the right number of QSO per square degrees
rand_qso_nb = 0.006

lya = 1215.67 ## angstrom https://en.wikipedia.org/wiki/Lyman-alpha_line  1215.668 and 1215.674
# lylimit = lya * 3 /4. # 911.75  ok with https://en.wikipedia.org/wiki/Hydrogen_spectral_series
lylimit = 0.  # do not cut pixels bellow lylimit
lyb = lylimit * 9./8. # 1025.72  https://en.wikipedia.org/wiki/Hydrogen_spectral_series : 1025.7
# lambda_min = 3530.  # Lambda min in A (given by Stephen Bailey - 04/05/2018)
lambda_min = 3476.

# Cosmo params given by Helion in a mail (30/03/2018)
# should correspond to Planck 2015
h = 0.6731
omega_M_0 = 0.31457
omega_b_0 = 0.045  # only use in run_camb.py
omega_lambda_0 = 0.68543
omega_k_0 = 0.0
ns = 0.96  # only used in run_camb.py

QSO_bias = 3.7  # should become b(z) in data.py / QSO.py
z_QSO_bias_1 = 1.9  # QSO_bias taken from Helion thesis
z_QSO_bias_2 = 2.75
z_QSO_bias_3 = 3.6
z0 = 1.70975268202
# box_pixel = 1600./1024.   # rather be stored in box.fits

H0 = 100.  # In km/s /(Mpc/h)

deg2rad = sp.pi/180.
rad2deg = 180/sp.pi

# Hardcoded params:
sigma_g = 1.19
rho_sum = 16452460
