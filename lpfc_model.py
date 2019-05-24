import numpy as np
import pylab
from scipy.stats import norm

def f(I):
    return (a*I-b)/(1-np.exp(-d*(a*I-b)))/1000.

def degdiff(n1,n2):
    """degdiff returns the angle in degree between the target orientations n1 and n2"""
    return min(abs(n1-n2),n-abs(n1-n2))*360./n

def anglediff(a1,a2):
    return min(abs(a1-a2),360-abs(a1-a2))

    """
def gauss(x,sigma=sigma):
    return np.exp(-x**2/(2*sigma**2))
"""

def gauss(deg1,deg2=0.,sigma=sigma):
    x = min(abs(deg1-deg2),360-abs(deg1-deg2))
    return np.exp(-x**2/(2*sigma**2))

def visgauss(t1,t2,sigma_s,npeak=1,n=n,bamp=0.):
    vgauss = np.zeros((2,n))
    if npeak == 1:
        for k in range(n):
            vgauss[0,k] = gauss(degdiff(k,t1*n/360.),sigma=sigma_s)
            vgauss[1,k] = gauss(degdiff(k,t2*n/360.),sigma=sigma_s)
    elif npeak == 2:
        for k in range(n):
            vgauss[0,k] = gauss(degdiff(k,t1*n/360.),sigma=sigma_s) + bamp*gauss(degdiff(k,t2*n/360.),sigma=sigma_s)
            vgauss[1,k] = gauss(degdiff(k,t2*n/360.),sigma=sigma_s) + bamp*gauss(degdiff(k,t1*n/360.),sigma=sigma_s)
    return vgauss

# Define the connectivity matrix
def define_connect(seednet=seednet):
    print('Define connectivity (conserved)')
    rng = np.random.RandomState(seednet)
    gau = np.ones((n,n))
    for k in range(n):
        for l in range(k+1,n):
            gau[k,l] = gauss(degdiff(0,abs(k-l)))
            gau[l,k] = gauss(degdiff(0,abs(k-l)))
    JM1 = J_m/2. + (1-g_tt/2.)*J_p*gau + std_J*rng.randn(n,n)
    JM2 = J_m/2. + (1-g_tt/2.)*J_p*gau + std_J*rng.randn(n,n)
    JMct = J_m + J_p*gau + std_J*rng.randn(n,n)
    JMx1 = J_m/2. + g_tt*J_p*gau/2. + std_J*rng.randn(n,n)
    JMx2 = J_m/2. + g_tt*J_p*gau/2. + std_J*rng.randn(n,n)
    return JM1,JM2,JMct,JMx1,JMx2

def wm_trace(CJ):
    rng = np.random.RandomState(232)
    # local variable below
    preoffer = 1500
    dur = preoffer + offeron + offeroff + targeton + go
    wmS = np.zeros(2)
    wmI = np.zeros(2)
    wmIno = np.zeros(2)
    wm = np.zeros(2)
    fwm = open(datapath+'lpfc_wm_CJ_'+str(CJ)+'_.txt','w')
    for j in range(int(dur/ts)):
        # Chosen juice neurons at OFC
        if j > int(preoffer/ts) and j < int((preoffer+offeron)/ts):
            if CJ == 0:
                Ijw = jI
            elif CJ == 1:
                Ijw = np.array([jI[1],jI[0]])
        else:
            Ijw = np.zeros(2)
        # Define input current
        wmI[0] = Ijw[0] + g_ww*wmS[0] + g_wwi*wmS[1] + wmIno[0]
        wmI[1] = Ijw[1] + g_ww*wmS[1] + g_wwi*wmS[0] + wmIno[1]
        # Firing rate
        wm = f(wmI)
        if np.mod(j,int(tp/ts)) == 0:
            fwm.write(str(wm[0])+'\t'+str(wm[1])+'\n')
        # Modified Euler
        rn = rng.randn(2)
        stemp = wmS + ts*(-wmS/tau_s+gamma*(1-wmS)*wm)
        inotemp = wmIno - (wmIno-wmIno0)*ts/tau_n + sigma_n*rn*np.sqrt(ts/tau_n)
        Itemp = inotemp + Ijw + g_ww*wmS
        Itemp[0] += g_wwi*wmS[1]
        Itemp[1] += g_wwi*wmS[0]
        wmS = wmS + ts/2.*(-wmS/tau_s + gamma*(1-wmS)*wm - stemp/tau_s + gamma*(1-stemp)*f(Itemp))
        wmIno = wmIno + ts/2.*(-(wmIno-wmIno0)/tau_n-(inotemp-wmIno0)/tau_n) + sigma_n*rn*np.sqrt(ts/tau_n)
    fwm.close()

# Trial
def trial(path,dir,CJ,npeak=1,bamp=0.,opp=None,output=0):
    # Initialization
    rng = np.random.RandomState(seedrun)
    target = np.zeros(2)
    visI = np.zeros((2,n))
    wmS = np.zeros(2)
    wmI = np.zeros(2)
    wmIno = np.zeros(2)
    wm = np.zeros(2)
    tgS = np.zeros((2,n))
    tgI = np.zeros((2,n))
    tgIno = np.zeros((2,n))
    tg = np.zeros((2,n))
    ctS = np.zeros(n)
    ctI = np.zeros(n)
    ctIno = np.zeros(n)
    ct = np.zeros(n)
    sigma_n_tg = sigma_n + std_sig*rng.randn(2,n)
    tgIno0 = ctIno0 + std_I*rng.randn(2,n)
    
    # Files for writing
    ft1 = open(path+'/lpfc_t1_net_'+str(seednet)+'_run_'+str(seedrun)+'_CJ_'+str(CJ)+'_d_'+str(dir)+'_.txt','w')
    ft2 = open(path+'/lpfc_t2_net_'+str(seednet)+'_run_'+str(seedrun)+'_CJ_'+str(CJ)+'_d_'+str(dir)+'_.txt','w')
    fct = open(path+'/lpfc_ct_net_'+str(seednet)+'_run_'+str(seedrun)+'_CJ_'+str(CJ)+'_d_'+str(dir)+'_.txt','w')

    # Target direction
    target[0] = dir
    if opp == None:
        target[1] = pylab.mod(target[0]+180,360)
    else:
        target[1] = opp
    
    # Visual neurons activated by target onset
    vgauss = visgauss(target[0],target[1],sigma_s,npeak,n,bamp)

    # Simulation
    for j in range(int(dur/ts)):
        # Chosen juice neurons at OFC
        if j > int(preoffer/ts) and j < int((preoffer+offeron)/ts):
            if CJ == 0:
                Ijw = jI
            elif CJ == 1:
                Ijw = np.array([jI[1],jI[0]])
        else:
            Ijw = np.zeros(2)
    
        # Define input current
        wmI[0] = Ijw[0] + g_ww*wmS[0] + g_wwi*wmS[1] + wmIno[0]
        wmI[1] = Ijw[1] + g_ww*wmS[1] + g_wwi*wmS[0] + wmIno[1]

        tgI = tgIno + g_wt*np.dot(wmS[np.newaxis,:].T,np.ones(n)[np.newaxis,:]) + g_st*visI + np.array([JM1.dot(tgS[0,:])/n, JM2.dot(tgS[1,:])/n])
        tgI += np.array([JMx1.dot(tgS[1,:])/n, JMx2.dot(tgS[0,:])/n])
        ctI = ctIno + JMct.dot(ctS)/n + g_tc*(tgS[0,:]+tgS[1,:])
        
        # Visual neurons activated by target onset
        if cue_on == 1:
            if j >= (preoffer+offeron+offeroff)/ts:
                visI = vgauss + std_vis*rng.randn(2,n)
            else:
                visI = np.zeros(n)
        else:
            if j >= (preoffer+offeron+offeroff)/ts and j < (preoffer+offeron+offeroff+targeton)/ts:
                visI = vgauss + std_vis*rng.randn(2,n)
            else:
                visI = np.zeros(n)

        # Firing rate
        wm = f(wmI)
        tg = f(tgI)
        ct = f(ctI)
        
        # Print files
        if np.mod(j,int(tp/ts)) == 0 and j>= tgt_start/ts:
            for k in range(n):
                ft1.write(str(tg[0,k])+'\t')
                ft2.write(str(tg[1,k])+'\t')
                fct.write(str(ct[k])+'\t')
            ft1.write('\n')
            ft2.write('\n')
            fct.write('\n')
        
        # Modified Euler
        rn = rng.randn(2)
        stemp = wmS + ts*(-wmS/tau_s+gamma*(1-wmS)*wm)
        inotemp = wmIno - (wmIno-wmIno0)*ts/tau_n + sigma_n*rn*np.sqrt(ts/tau_n)
        Itemp = inotemp + Ijw + g_ww*wmS
        Itemp[0] += g_wwi*wmS[1]
        Itemp[1] += g_wwi*wmS[0]
        wmS = wmS + ts/2.*(-wmS/tau_s + gamma*(1-wmS)*wm - stemp/tau_s + gamma*(1-stemp)*f(Itemp))
        wmIno = wmIno + ts/2.*(-(wmIno-wmIno0)/tau_n-(inotemp-wmIno0)/tau_n) + sigma_n*rn*np.sqrt(ts/tau_n)
        
        rn = rng.randn(2,n)
        stemp = tgS + ts*(-tgS/tau_s+gamma*(1-tgS)*tg)
        inotemp = tgIno - (tgIno-tgIno0)*ts/tau_n + sigma_n_tg*rn*np.sqrt(ts/tau_n)
        Itemp = inotemp + np.array([JM1.dot(tgS[0,:])/n, JM2.dot(tgS[1,:])/n]) + g_wt*np.dot(wmS[np.newaxis,:].T,np.ones(n)[np.newaxis,:]) + g_st*visI
        Itemp += np.array([JMx1.dot(tgS[1,:])/n, JMx2.dot(tgS[0,:])/n])
        tgS = tgS + ts/2.*(-tgS/tau_s + gamma*(1-tgS)*tg - stemp/tau_s + gamma*(1-stemp)*f(Itemp))
        tgIno = tgIno + ts/2.*(-(tgIno-tgIno0)/tau_n-(inotemp-tgIno0)/tau_n) + sigma_n_tg*rn*np.sqrt(ts/tau_n)
        
        rn = rng.randn(n)
        stemp = ctS + ts*(-ctS/tau_s+gamma*(1-ctS)*ct)
        inotemp = ctIno - (ctIno-ctIno0)*ts/tau_n + sigma_n_ring*rn*np.sqrt(ts/tau_n)
        Itemp = inotemp + JMct.dot(ctS)/n + g_tc*(tgS[0,:]+tgS[1,:])
        ctS = ctS + ts/2.*(-ctS/tau_s + gamma*(1-ctS)*ct - stemp/tau_s + gamma*(1-stemp)*f(Itemp))
        ctIno = ctIno + ts/2.*(-(ctIno-ctIno0)/tau_n-(inotemp-ctIno0)/tau_n) + sigma_n_ring*rn*np.sqrt(ts/tau_n)
    temp,decodedCT = resultant(ct,div=n)
    print 'CT = '+str(target[CJ])+'; Decoded CT = '+str(decodedCT)
    if anglediff(target[CJ],decodedCT) > 22.5:
        print 'Please check the result!'
    ft1.close()
    ft2.close()
    fct.close()
    if output == 1:
        return np.min([np.abs(target[CJ]-decodedCT),360-np.abs(target[CJ]-decodedCT)])

def outfigure(path,dir,CJ):
    pylab.figure(figsize=[9,6])
    pylab.subplot(411)
    item = pylab.loadtxt(datapath+'lpfc_wm_CJ_'+str(CJ)+'_.txt')
    pylab.plot(np.arange(tgt_start,dur,tp),item[-(dur-tgt_start)/tp:,0]*1000)
    pylab.plot(np.arange(tgt_start,dur,tp),item[-(dur-tgt_start)/tp:,1]*1000)
    pylab.ylabel('WM')
    pylab.yticks([0,30,60])
    pylab.xticks([])
    pylab.xlim(tgt_start,dur)
    
    pylab.subplot(412)
    item = pylab.loadtxt(path+'/lpfc_t1_net_'+str(seednet)+'_run_'+str(seedrun)+'_CJ_'+str(CJ)+'_d_'+str(dir)+'_.txt')
    pylab.pcolor(item.T*1000,rasterized=True,cmap='jet')
    pylab.clim(0,40)
    pylab.colorbar(ticks=[0,20,40],cmap='jet')
    pylab.ylabel('TG-A')
    pylab.yticks(range(0,n+1,n/2),['0$^{\circ}$','180$^{\circ}$','360$^{\circ}$'])
    pylab.xticks([])
    pylab.xlim(0,targeton+go)
    pylab.ylim(0,n)
    
    pylab.subplot(413)
    item2 = pylab.loadtxt(path+'/lpfc_t2_net_'+str(seednet)+'_run_'+str(seedrun)+'_CJ_'+str(CJ)+'_d_'+str(dir)+'_.txt')
    pylab.pcolor(item2.T*1000,rasterized=True,cmap='jet')
    pylab.clim(0,40)
    pylab.colorbar(ticks=[0,20,40])
    #pylab.colorbar(ticks=[0,20,40])
    #pylab.clim(0,20)
    #cbar = pylab.colorbar(ticks=[0,20])
    #mpl.rcParams['axes.labelsize'] = font_size+6
    #cbar.set_label('Rate (Hz)')
    pylab.ylabel('TG-B')
    pylab.yticks(range(0,n+1,n/2),['0$^{\circ}$','180$^{\circ}$','360$^{\circ}$'])
    pylab.xticks([])
    pylab.xlim(0,targeton+go)
    pylab.ylim(0,n)
    
    pylab.subplot(414)
    item = pylab.loadtxt(path+'/lpfc_ct_net_'+str(seednet)+'_run_'+str(seedrun)+'_CJ_'+str(CJ)+'_d_'+str(dir)+'_.txt')
    pylab.pcolor(item.T*1000,rasterized=True,cmap='jet')
    pylab.clim(0,40)
    pylab.colorbar(ticks=[0,20,40])
    pylab.ylabel('CT')
    pylab.yticks(range(0,n+1,n/2),['0$^{\circ}$','180$^{\circ}$','360$^{\circ}$'])
    pylab.xticks([])
    pylab.xlim(0,targeton+go)
    pylab.ylim(0,n)
    pylab.xticks([0,1000])
    pylab.subplot(411)
    pylab.colorbar(ticks=[0,20,40])
    pylab.savefig(path+'/lpfc_CJ_'+str(CJ)+'_d_'+str(dir)+'_.png')
    pylab.close('all')

def resultant(vect,m=1,div=8):
    import cmath
    if any(vect<0):
        print 'vect < 0'
    s = 0
    for j in range(div):
        s = s + vect[j]*np.exp(complex(0,m*2.*np.pi*j/div))
    s = s/sum(vect)
    if cmath.phase(s) < 0:
        return abs(s),(2*np.pi + cmath.phase(s))/m*180/np.pi
    else:
        return abs(s),cmath.phase(s)/m*180/np.pi

def prefdir(vect,div=8):
    [r1,a1] = resultant(vect,1,div)
    [r2,a2] = resultant(vect,2,div)
    if r1 >= r2:
        return a1
    else:
        if degreediff(a1,a2,1) < degreediff(a1,a2+180,1):
            return a2
        else:
            return a2+180

def decode(which_ring,t1=900,t2=1000):
    decoded_arr = []
    for CJ in range(2):
        for dir in range(0,360,45):
            fdata = pylab.loadtxt(datapath+'lpfc_'+which_ring+'_net_'+str(seednet)+'_run_'+str(seedrun)+'_CJ_'+str(CJ)+'_d_'+str(dir)+'_.txt')
            data = fdata[(preoffer+offeron+offeroff+t1)/tp:(preoffer+offeron+offeroff+t2)/tp,:]
            temp,decoded = resultant(np.mean(data,0),div=n)
            print 'Actual = '+str(dir)+'; Decoded = '+str(decoded)
            decoded_arr.append(decoded)
    return decoded_arr

def degreediff(a1,a2,mode=2):
    if mode == 0:    # [-180 180]
        return np.mod(np.mod(a2-a1,360)-180,360)-180
    elif mode == 1:    # absolute [0 360]
        return min(abs(a1-a2),360-abs(a1-a2))
    elif mode == 2:    # [-90 270]
        return np.mod(np.mod(a2-a1,360)-270,360)-90

def time_transition(path,CJ=0,dir=90):
    decoded = []
    for j in range(0,800):
        ft2 = pylab.loadtxt(path+'/lpfc_t2_net_'+str(seednet)+'_run_'+str(seedrun)+'_CJ_'+str(CJ)+'_d_'+str(dir)+'_.txt')
        t2 = ft2[j:j+200,:]
        decoded.append(prefdir(np.mean(t2,0),n))
    #pylab.figure()
    #pylab.plot(range(100,900),decoded)
    #pylab.xlabel('Time after target onset')
    #pylab.ylabel('Decoded direction')
    Ttran = np.argmax(np.abs(np.diff(decoded))) + 1
    if abs(decoded[Ttran-1]-decoded[Ttran]) < 90:
        return 0
    else:
        return 1000./(Ttran+100)

def tuning_curves(path,Ltime=400):
    tuningE0 = np.zeros((n*3,8))
    tuningE1 = np.zeros((n*3,8))
    tuningL0 = np.zeros((n*3,8))
    tuningL1 = np.zeros((n*3,8))
    for dir in range(0,360,45):
        ft1 = pylab.loadtxt(path+'/lpfc_t1_net_'+str(seednet)+'_run_'+str(seedrun)+'_CJ_0_d_'+str(dir)+'_.txt')
        ft2 = pylab.loadtxt(path+'/lpfc_t2_net_'+str(seednet)+'_run_'+str(seedrun)+'_CJ_0_d_'+str(dir)+'_.txt')
        fct = pylab.loadtxt(path+'/lpfc_ct_net_'+str(seednet)+'_run_'+str(seedrun)+'_CJ_0_d_'+str(dir)+'_.txt')
        tuningE0[:,dir/45] = np.append(np.append(np.mean(ft1[:200,:],0),np.mean(ft2[:200,:],0)),np.mean(fct[:200,:],0))
        tuningL0[:,dir/45] = np.append(np.append(np.mean(ft1[Ltime:Ltime+200,:],0),np.mean(ft2[Ltime:Ltime+200,:],0)),np.mean(fct[Ltime:Ltime+200,:],0))
        ft1 = pylab.loadtxt(path+'/lpfc_t1_net_'+str(seednet)+'_run_'+str(seedrun)+'_CJ_1_d_'+str(dir)+'_.txt')
        ft2 = pylab.loadtxt(path+'/lpfc_t2_net_'+str(seednet)+'_run_'+str(seedrun)+'_CJ_1_d_'+str(dir)+'_.txt')
        fct = pylab.loadtxt(path+'/lpfc_ct_net_'+str(seednet)+'_run_'+str(seedrun)+'_CJ_1_d_'+str(dir)+'_.txt')
        tuningE1[:,dir/45] = np.append(np.append(np.mean(ft1[:200,:],0),np.mean(ft2[:200,:],0)),np.mean(fct[:200,:],0))
        tuningL1[:,dir/45] = np.append(np.append(np.mean(ft1[Ltime:Ltime+200,:],0),np.mean(ft2[Ltime:Ltime+200,:],0)),np.mean(fct[Ltime:Ltime+200,:],0))
    return tuningE0, tuningE1, tuningL0, tuningL1

def tuning_cj(path,t_start):
    tuning0 = np.zeros((n*3,8))
    tuning1 = np.zeros((n*3,8))
    for dir in range(0,360,45):
        ft1 = pylab.loadtxt(path+'/lpfc_t1_net_'+str(seednet)+'_run_'+str(seedrun)+'_CJ_0_d_'+str(dir)+'_.txt')
        ft2 = pylab.loadtxt(path+'/lpfc_t2_net_'+str(seednet)+'_run_'+str(seedrun)+'_CJ_0_d_'+str(dir)+'_.txt')
        fct = pylab.loadtxt(path+'/lpfc_ct_net_'+str(seednet)+'_run_'+str(seedrun)+'_CJ_0_d_'+str(dir)+'_.txt')
        tuning0[:,dir/45] = np.append(np.append(np.mean(ft1[t_start:t_start+200,:],0),np.mean(ft2[t_start:t_start+200,:],0)),np.mean(fct[t_start:t_start+200,:],0))
        ft1 = pylab.loadtxt(path+'/lpfc_t1_net_'+str(seednet)+'_run_'+str(seedrun)+'_CJ_1_d_'+str(dir)+'_.txt')
        ft2 = pylab.loadtxt(path+'/lpfc_t2_net_'+str(seednet)+'_run_'+str(seedrun)+'_CJ_1_d_'+str(dir)+'_.txt')
        fct = pylab.loadtxt(path+'/lpfc_ct_net_'+str(seednet)+'_run_'+str(seedrun)+'_CJ_1_d_'+str(dir)+'_.txt')
        tuning1[:,dir/45] = np.append(np.append(np.mean(ft1[t_start:t_start+200,:],0),np.mean(ft2[t_start:t_start+200,:],0)),np.mean(fct[t_start:t_start+200,:],0))
    return tuning0, tuning1

"""
def tuning_cj2(path,t_start):
    # Not used
    tuning0 = np.zeros((n*3,8))
    tuning1 = np.zeros((n*3,8))
    for dir in range(0,360,45):
        ft1 = pylab.loadtxt(path+'/lpfc_t1_net_'+str(seednet)+'_run_0_CJ_0_d_'+str(dir)+'_.txt')
        ft2 = pylab.loadtxt(path+'/lpfc_t2_net_'+str(seednet)+'_run_0_CJ_0_d_'+str(dir)+'_.txt')
        fct = pylab.loadtxt(path+'/lpfc_ct_net_'+str(seednet)+'_run_0_CJ_0_d_'+str(dir)+'_.txt')
        tuning0[:,dir/45] = np.append(np.append(np.mean(ft1[t_start:t_start+200,:],0),np.mean(ft2[t_start:t_start+200,:],0)),np.mean(fct[t_start:t_start+200,:],0))
        ft1 = pylab.loadtxt(path+'/lpfc_t1_net_'+str(seednet)+'_run_0_CJ_1_d_'+str(dir)+'_.txt')
        ft2 = pylab.loadtxt(path+'/lpfc_t2_net_'+str(seednet)+'_run_0_CJ_1_d_'+str(dir)+'_.txt')
        fct = pylab.loadtxt(path+'/lpfc_ct_net_'+str(seednet)+'_run_0_CJ_1_d_'+str(dir)+'_.txt')
        tuning1[:,dir/45] = np.append(np.append(np.mean(ft1[t_start:t_start+200,:],0),np.mean(ft2[t_start:t_start+200,:],0)),np.mean(fct[t_start:t_start+200,:],0))
    return tuning0, tuning1

def tuning_at_a_trial(path,t_start):
    tuning0 = np.zeros((n*3,8))
    tuning1 = np.zeros((n*3,8))
    for dir in range(0,360,45):
        ft1 = pylab.loadtxt(path+'/lpfc_t1_net_'+str(seednet)+'_run_'+str(seedrun)+'_CJ_0_d_'+str(dir)+'_.txt')
        ft2 = pylab.loadtxt(path+'/lpfc_t2_net_'+str(seednet)+'_run_'+str(seedrun)+'_CJ_0_d_'+str(dir)+'_.txt')
        fct = pylab.loadtxt(path+'/lpfc_ct_net_'+str(seednet)+'_run_'+str(seedrun)+'_CJ_0_d_'+str(dir)+'_.txt')
        tuning0[:,dir/45] = np.append(np.append(np.mean(ft1[t_start:t_start+200,:],0),np.mean(ft2[t_start:t_start+200,:],0)),np.mean(fct[t_start:t_start+200,:],0))
        ft1 = pylab.loadtxt(path+'/lpfc_t1_net_'+str(seednet)+'_run_'+str(seedrun)+'_CJ_1_d_'+str(dir)+'_.txt')
        ft2 = pylab.loadtxt(path+'/lpfc_t2_net_'+str(seednet)+'_run_'+str(seedrun)+'_CJ_1_d_'+str(dir)+'_.txt')
        fct = pylab.loadtxt(path+'/lpfc_ct_net_'+str(seednet)+'_run_'+str(seedrun)+'_CJ_1_d_'+str(dir)+'_.txt')
        tuning1[:,dir/45] = np.append(np.append(np.mean(ft1[t_start:t_start+200,:],0),np.mean(ft2[t_start:t_start+200,:],0)),np.mean(fct[t_start:t_start+200,:],0))
    return tuning0, tuning1
"""

def peak_tuning(tuningE0, tuningE1, tuningL0, tuningL1):
    print 'Running peak_tunng: find the peaks of the 4 tuning curves'
    p = np.zeros((np.size(tuningE0,0),4))
    peakdata = np.zeros((np.size(tuningE0,0),4))
    #f = open(path+'peak_selected.txt','w')
    idx = []
    for j in range(np.size(tuningE0,0)):
        p[j,0] = prefdir(tuningE0[j,:])
        p[j,1] = prefdir(tuningE1[j,:])
        p[j,2] = prefdir(tuningL0[j,:])
        p[j,3] = prefdir(tuningL1[j,:])
        peakdata[j,0] = degreediff(p[j,0],p[j,1])
        peakdata[j,1] = degreediff(p[j,2],p[j,3])
        peakdata[j,2] = degreediff(p[j,0],p[j,2])
        peakdata[j,3] = degreediff(p[j,1],p[j,3])
        """
        if sum(np.isnan(peakdata[j,:])) == 0:
            f.write(str(peakdata[j,0])+'\t'+str(peakdata[j,1])+'\t'+str(peakdata[j,2])+'\t'+str(peakdata[j,3]))
            f.write('\n')
            idx.append(j)
            """
    #print 'peak_tuning: there are '+str(len(idx))+' valid sets of peaks.'
    peakdata = peakdata[np.sum(np.isnan(peakdata),1)==0]
    #return p[idx,:],peakdata[idx,:]
    return peakdata

def clustering(peakdata,eps=eps,min_samples=min_samples,epilson=epilson):
    from sklearn.cluster import DBSCAN
    db = DBSCAN(eps=eps, min_samples=min_samples).fit(peakdata)
    idx = db.labels_+1
    ncount = np.zeros(4)
    ncheck = np.zeros(4)
    uc = []
    nuc = []
    for j in range(1,max(idx)+1):
        if abs(np.mean(peakdata[idx==j,0]))<epilson and abs(np.mean(peakdata[idx==j,1]))<epilson and abs(np.mean(peakdata[idx==j,2]))<epilson and abs(np.mean(peakdata[idx==j,3]))<epilson:
            ncount[0] = peakdata[idx==j,0].size
            ncheck[0] = j
        elif abs(np.mean(peakdata[idx==j,0])-180)<epilson and abs(np.mean(peakdata[idx==j,1])-180)<epilson and abs(np.mean(peakdata[idx==j,2]))<epilson and abs(np.mean(peakdata[idx==j,3]))<epilson:
            ncount[1] = peakdata[idx==j,0].size
            ncheck[1] = j
        elif abs(np.mean(peakdata[idx==j,0]))<epilson and abs(np.mean(peakdata[idx==j,1])-180)<epilson and abs(np.mean(peakdata[idx==j,2])-180)<epilson and abs(np.mean(peakdata[idx==j,3]))<epilson:
            ncount[2] = peakdata[idx==j,0].size
            ncheck[2] = j
        elif abs(np.mean(peakdata[idx==j,0]))<epilson and abs(np.mean(peakdata[idx==j,1])-180)<epilson and abs(np.mean(peakdata[idx==j,2]))<epilson and abs(np.mean(peakdata[idx==j,3])-180)<epilson:
            ncount[3] = peakdata[idx==j,0].size
            ncheck[3] = j
        else:
            uc.append(j)
            nuc.append(peakdata[idx==j,0].size)
    return idx,ncount,ncheck,uc,nuc

def cov_ellipse(points, nstd=1, ax=None, **kwargs):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse
    pos = points.mean(axis=0)
    cov = np.cov(points, rowvar=False)
    def eigsorted(cov):
        vals, vecs = np.linalg.eigh(cov)
        order = vals.argsort()[::-1]
        return vals[order], vecs[:,order]        
    if ax is None:
        ax = plt.gca()
    vals, vecs = eigsorted(cov)
    theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))
    width, height = 2 * nstd * np.sqrt(vals)
    ellip = Ellipse(xy=pos, width=width, height=height, angle=theta, **kwargs)
    ax.add_artist(ellip)
    return pos

def compp(p1,p2,n1,n2):
    """ p1 - proportion in population 1
        p2 - proportion in population 2
        n1 - number of neurons for population 1
        n2 - number of neurons for population 2"""
    pboth=(p1*n1+p2*n2)/(n1+n2)
    z = -abs(p1-p2)/np.sqrt(pboth*(1-pboth)*(1./n1+1./n2))
    p = 2*norm.cdf(z)
    return p
