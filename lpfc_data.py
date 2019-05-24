import numpy as np
import pylab
import matplotlib as mpl
import os.path

# WM
J_p,J_m,g_wt,g_tt,wmIno0,ctIno0,std_J,std_t,std_I,std_sig = para_homogeneous_scene1()
for CJ in [0,1]:
    if not os.path.isfile(datapath+'lpfc_wm_CJ_'+str(CJ)+'_.txt'):
        print 'WM trace'
        wm_trace(CJ)
        
# Homogeneous network
J_p,J_m,g_wt,g_tt,wmIno0,ctIno0,std_J,std_t,std_I,std_sig = para_homogeneous_scene1()
for g_tt in np.arange(0.,1.01,0.1):
    fname = 'jp'+str(J_p)+'jm'+str(-J_m)+'gwt'+str(g_wt)+'gtt'+str(g_tt)+'I'+str(wmIno0)+'std'+str(std_J)+'seednet'+str(seednet)+'seedrun'+str(seedrun)
    if not os.path.exists(datapath+fname) or overwrite:
        JM1,JM2,JMct,JMx1,JMx2 = define_connect()
        if not os.path.exists(datapath+fname):
            os.makedirs(datapath+fname)
        print 'Running simulations for '+fname
        for CJ in range(2):
            for dir in range(0,360,45):
                print 'dir: '+str(dir)
                trial(datapath+fname,dir,CJ)

J_p_array,J_m_array = para_homogeneous_scan()
for [J_p,J_m] in np.array([J_p_array,J_m_array]).T:
    for g_tt in np.arange(0.,1.01,0.1):
        fname = 'jp'+str(J_p)+'jm'+str(-J_m)+'gwt'+str(g_wt)+'gtt'+str(g_tt)+'I'+str(wmIno0)+'std'+str(std_J)+'seednet'+str(seednet)+'seedrun'+str(seedrun)
        if not os.path.exists(datapath+fname) or overwrite:
            JM1,JM2,JMct,JMx1,JMx2 = define_connect()
            if not os.path.exists(datapath+fname):
                os.makedirs(datapath+fname)
            print 'Running simulations for '+fname
            for CJ in range(2):
                for dir in range(0,360,45):
                    print 'dir: '+str(dir)
                    trial(datapath+fname,dir,CJ)
                    #outfigure(datapath+fname,dir,CJ)


# Heterogeneous network
J_p,J_m,g_wt,g_tt,wmIno0,ctIno0,std_J,std_t,std_I,std_sig = para_heterogeneous()
for std_J in std_J_arr[:-1]:
    for g_tt in np.arange(0.4,1.1,0.1):
        fname = 'jp'+str(J_p)+'jm'+str(-J_m)+'gwt'+str(g_wt)+'gtt'+str(g_tt)+'I'+str(wmIno0)+'std'+str(std_J)+'seednet'+str(seednet)+'seedrun'+str(seedrun)
        if not os.path.exists(datapath+fname) or overwrite:
            JM1,JM2,JMct,JMx1,JMx2 = define_connect()
            if not os.path.exists(datapath+fname):
                os.makedirs(datapath+fname)
            print 'Running simulations for '+fname
            for CJ in range(2):
                for dir in range(0,360,45):
                    print 'dir: '+str(dir)
                    trial(datapath+fname,dir,CJ)
                    #outfigure(datapath+fname,dir,CJ)

# Heterogeneous network
J_p,J_m,g_wt,g_tt,wmIno0,ctIno0,std_J,std_t,std_I,std_sig = para_heterogeneous()
for std_J in std_J_arr:
    fname = 'jp'+str(J_p)+'jm'+str(-J_m)+'gwt'+str(g_wt)+'gtt'+str(g_tt)+'I'+str(wmIno0)+'std'+str(std_J)+'seednet'+str(seednet)+'seedrun'+str(seedrun)
    if not os.path.exists(datapath+fname) or overwrite:
        JM1,JM2,JMct,JMx1,JMx2 = define_connect()
        if not os.path.exists(datapath+fname):
            os.makedirs(datapath+fname)
        print 'Running simulations for '+fname
        for CJ in range(2):
            for dir in range(0,360,45):
                print 'dir: '+str(dir)
                trial(datapath+fname,dir,CJ)
                #outfigure(datapath+fname,dir,CJ)

# Heterogeneous neurons
J_p,J_m,g_wt,g_tt,wmIno0,ctIno0,std_J,std_t,std_I,std_sig = para_heterogeneous()
std_J = 2.
for g_tt in np.arange(0.4,1.1,0.1):
    for std_Is in np.arange(0.,0.026,0.005):
        std_I = std_Is
        std_sig = std_Is
        fname = 'jp'+str(J_p)+'jm'+str(-J_m)+'gwt'+str(g_wt)+'gtt'+str(g_tt)+'I'+str(wmIno0)+'std'+str(std_J)+'stdI'+str(std_I)+'stds'+str(std_sig)+'seednet'+str(seednet)+'seedrun'+str(seedrun)
        if not os.path.exists(datapath+fname) or overwrite:
            JM1,JM2,JMct,JMx1,JMx2 = define_connect()
            if not os.path.exists(datapath+fname):
                os.makedirs(datapath+fname)
            print 'Running simulations for '+fname
            for CJ in range(2):
                for dir in range(0,360,45):
                    print 'dir: '+str(dir)
                    trial(datapath+fname,dir,CJ)
                    #outfigure(datapath+fname,dir,CJ)
