### Parameters
import numpy as np
import pylab
import matplotlib as mpl

font_size = 12
mpl.rcParams['axes.titlesize'] = font_size+2
mpl.rcParams['xtick.labelsize'] = font_size
mpl.rcParams['ytick.labelsize'] = font_size
mpl.rcParams['axes.labelsize'] = font_size
mpl.rcParams['legend.fontsize'] = font_size-2
mpl.rcParams['font.size'] = font_size
fl_size = 18

# f-I
a = 270.
b = 108.
d = 0.154
gamma = 0.641
# NMDA
tau_s = 60.
# noise
tau_n = 2.
sigma_n = 0.015
sigma_n_ring = sigma_n
# ring
sigma = 43.2
sigma_s = sigma
n = 256
# visual
Ndir = 8
# connection strength
g_ww = 0.3725
g_wwi = -0.1137
g_st = 0.1
g_tc = 0.09
std_t = 0.
std_sig = 0.
std_vis = 0.
std_I = 0.
std_J_arr = [0.,0.5,1.,1.5,1.75,2.,2.25,2.5,2.75]

# Simulation parameters
cue_on = 0
preoffer = 500
offeron = 1000
offeroff = 1000
targeton = 1000
go = 0
tgt_start = preoffer + offeron + offeroff
dur = tgt_start + targeton + go
ts = 0.5
tp = 1
# The timestep for saving data. The codes for analysis is based on tp = 1

# Clustering
eps = 65
eps_sim = 50
min_samples = 20
epilson = 30

# Color
cA = '#ff8c00'
cB = '#1e90ff'
cTG = '#f08080'
cCT = '#3cb371'
cTS1 = '#8a2be2'
cTS2 = '#b8860b'
cEarly = '#191970'
cMid = '#006400'
cLate = '#a52a2a'
cUN = ['b','g','r','c','y','k','k','k','k','k','k','k']
alpha = 0.2
#foldername = ['jp2jm0_5gtt50','jp2jm0_5gtt90','jp2_32jm0_8gtt90std250']
tgt = 90
opp = 270
alldir = np.arange(0.,360.1,360./n)

### Initialization
jI = np.array([0.03,0.015])
overwrite = 0

# Random seed specified (for OU noisy input)
seednet = 101
seedrun = 51

# Path
matlabpath = 'matlab/'
datapath = 'datafiles/'
figpath = 'figures/'

# Homogeneous network
def para_homogeneous_scene1():
    J_p = 2.
    J_m = -0.35
    g_wt = 0.01
    g_tt = 0.
    wmIno0 = 0.3197
    ctIno0 = wmIno0
    std_J = 0.
    std_t = 0.
    std_I = 0.
    std_sig = 0.
    return J_p,J_m,g_wt,g_tt,wmIno0,ctIno0,std_J,std_t,std_I,std_sig

def para_homogeneous_scene2():
    J_p,J_m,g_wt,g_tt,wmIno0,ctIno0,std_J,std_t,std_I,std_sig = para_homogeneous_scene1()
    g_tt = 1.
    return J_p,J_m,g_wt,g_tt,wmIno0,ctIno0,std_J,std_t,std_I,std_sig

def para_homogeneous_scan():
    return np.array([1.9,1.9,2.02]),np.array([-0.6,-0.35,-0.35])

# Heterogeneous network
def para_heterogeneous():
    J_p = 2.32
    J_m = -0.8
    g_wt = 0.03
    g_tt = 0.9
    wmIno0 = 0.3297
    ctIno0 = wmIno0
    std_J = 2.5
    std_t = 0.
    std_I = 0.
    std_sig = 0.
    return J_p,J_m,g_wt,g_tt,wmIno0,ctIno0,std_J,std_t,std_I,std_sig

# Heterogeneous neurons
def para_heterogeneous_neuron():
    J_p,J_m,g_wt,g_tt,wmIno0,ctIno0,std_J,std_t,std_I,std_sig = para_heterogeneous()
    std_J = 2.
    std_I = 0.025
    std_sig = std_I
    return J_p,J_m,g_wt,g_tt,wmIno0,ctIno0,std_J,std_t,std_I,std_sig
