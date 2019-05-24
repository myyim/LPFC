import numpy as np
import pylab
import matplotlib as mpl
import scipy.interpolate as sp
import os.path

execfile('lpfc_para.py')
execfile('lpfc_model.py')

if not os.path.exists(figpath+'s0.pdf'):
    print 'Making supplementary figure 0...'
    seedall = [[101,51],[3,56],[42,36],[55,9],[32,332]]
    ddir = 180
    s1l = []
    s1h = []
    s2l = []
    s2h = []
    for sigma_n in [0.015,0.03]:
        sigma_n_ring = sigma_n
        for scene in [1,2]:
            if scene == 1:
                J_p,J_m,g_wt,g_tt,wmIno0,ctIno0,std_J,std_t,std_I,std_sig = para_homogeneous_scene1()
            elif scene == 2:
                J_p,J_m,g_wt,g_tt,wmIno0,ctIno0,std_J,std_t,std_I,std_sig = para_homogeneous_scene2()
            for seednet,seedrun in seedall:
                fname = 'input2pk_sig'+str(sigma)+'_ddir'+str(ddir)+'jp'+str(J_p)+'jm'+str(-J_m)+'gwt'+str(g_wt)+'gst'+str(g_st)+'gtt'+str(g_tt)+'I'+str(wmIno0)+'std'+str(std_J)+'seednet'+str(seednet)+'seedrun'+str(seedrun)
                if not os.path.exists(datapath+fname) or 1:
                    JM1,JM2,JMct,JMx1,JMx2 = define_connect()
                    if not os.path.exists(datapath+fname):
                        os.makedirs(datapath+fname)
                print 'Running simulations for '+fname
                for CJ in range(2):
                    for dir in range(0,360,45):
                        print 'dir: '+str(dir)
                        temp = trial(datapath+fname,dir,CJ,npeak=2,opp=pylab.mod(dir+ddir,360),output=1)
                        if sigma_n == 0.015:
                            if scene == 1:
                                s1l.append(temp)
                            elif scene == 2:
                                s2l.append(temp)
                        elif sigma_n == 0.03:
                            if scene == 1:
                                s1h.append(temp)
                            elif scene == 2:
                                s2h.append(temp)
                        #outfigure(datapath+fname,dir,CJ)
    pylab.figure()
    ax = pylab.subplot(111)
    pylab.plot(len(s1l)*[0],s1l,'x',c='0.7')
    pylab.plot(len(s1l)*[1],s2l,'x',c='0.7')
    pylab.plot(len(s1l)*[2],s1h,'x',c='0.7')
    pylab.plot(len(s1l)*[3],s2h,'x',c='0.7')
    pylab.plot(range(4),[np.mean(s1l),np.mean(s2l),np.mean(s1h),np.mean(s2h)],'kd')
    pylab.ylabel('Absolute error in degrees')
    pylab.xticks(range(4),('S1,low','S2,low','S1,high','S2,high'))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    pylab.savefig(figpath+'s0.pdf')
seednet = 101
seedrun = 51

def supp1():
    print 'Making supplementary figure 1...'
    fig = pylab.figure(figsize=[13,13])
    J_p,J_m,g_wt,g_tt,wmIno0,ctIno0,std_J,std_t,std_I,std_sig = para_homogeneous_scene2()
    J_p_arr,J_m_arr = para_homogeneous_scan()
    for j in range(3):
        for k in range(11):
            pylab.subplot(11,3,3*k+j+1)
            g_tt = 0.1*k
            fname = 'jp'+str(J_p_arr[j])+'jm'+str(-J_m_arr[j])+'gwt'+str(g_wt)+'gtt'+str(g_tt)+'I'+str(wmIno0)+'std'+str(std_J)+'seednet'+str(seednet)+'seedrun'+str(seedrun)
            item = pylab.loadtxt(datapath+fname+'/lpfc_t2_net_'+str(seednet)+'_run_'+str(seedrun)+'_CJ_0_d_90_.txt')
            item = item[0:int(targeton/tp)+1,:]
            pylab.pcolor(item.T*1000,rasterized=True,cmap='jet')
            pylab.clim(0,60)
            cbar = pylab.colorbar(ticks=[0,30,60])
            if j == 2 and k == 0:
                cbar.set_label('Rate (Hz)')
            else:
                cbar.remove()
            pylab.ylabel('IN-B')
            if j == 0:
                pylab.yticks(range(0,n+1,n/2),['0$^{\circ}$','180$^{\circ}$','360$^{\circ}$'])
            else:
                pylab.yticks([])
            pylab.xticks([])
            pylab.xlim(0,targeton)
            pylab.ylim(0,n)
            if k == 10:
                pylab.xticks([0,500,1000])
                pylab.xlabel('Time from target onset (ms)')
            if k == 0:
                if j == 0:
                    pylab.title('strong inh')
                elif j == 1:
                    pylab.title('reference')
                elif j == 2:
                    pylab.title('strong exc')
    #pylab.subplots_adjust(top=0.9,bottom=0.1,hspace=0.2,wspace=0.05)
    pylab.savefig(figpath+'s1.pdf')

def supp2():
    print 'Making supplementary figure 2...'
    font_size = 12
    mpl.rcParams['font.size'] = font_size
    tuningE0 = pylab.loadtxt(matlabpath+'tuningallA0.txt')
    tuningE1 = pylab.loadtxt(matlabpath+'tuningallB0.txt')
    tuningL0 = pylab.loadtxt(matlabpath+'tuningallA400.txt')
    tuningL1 = pylab.loadtxt(matlabpath+'tuningallB400.txt')
    peakdata = peak_tuning(tuningE0, tuningE1, tuningL0, tuningL1)
    dipp = pylab.loadtxt(matlabpath+'dipp_all.txt')
    fig = pylab.figure(figsize=[12,15])
    # D
    fig.text(0.05,0.47,'D',fontsize=fl_size,ha='right',va='bottom')
    pylab.subplot(437)
    pylab.scatter(peakdata[:,0],peakdata[:,1],c='k',s=30,lw=0)
    pylab.plot([-90,270],[-90,270],'k--')
    pylab.xlim(-95,275)
    pylab.ylim(-95,275)
    pylab.xticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.yticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.xlabel('$\Delta_{Early}$')
    pylab.ylabel('$\Delta_{Late}$')
    pylab.subplot(438)
    pylab.scatter(peakdata[:,2],peakdata[:,3],c='k',s=30,lw=0)
    pylab.plot([-90,270],[-90,270],'k--')
    pylab.xlim(-95,275)
    pylab.ylim(-95,275)
    pylab.xticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.yticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.xlabel('$\Delta_{A}$')
    pylab.ylabel('$\Delta_{B}$')
    # E
    fig.text(0.65,0.47,'E',fontsize=fl_size,ha='right',va='bottom')
    ax = pylab.subplot(8,6,29)
    pylab.hist(peakdata[:,0],pylab.arange(-90,271,22.5),color='k')
    #pylab.xticks(range(0,271,180),['0$^{\circ}$','180$^{\circ}$'])
    pylab.xticks([])
    pylab.yticks(range(0,160,100))
    pylab.xlim(-90,270)
    pylab.ylim(0,180)
    pylab.xlabel('$\Delta_{Early}$')
    #pylab.title('%s %.3f %s %.3f' % ('dip=',dipp[0,0],'; p=',dipp[0,1]));
    if dipp[0,1] == 0.:
        pylab.text(-75,150,'$p <$ 10$^{-5}$')
    else:
        pylab.text(-75,150,'%s %.3f' % ('$p = $',dipp[0,1]))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax = pylab.subplot(8,6,30)
    pylab.hist(peakdata[:,1],pylab.arange(-90,271,22.5),color='k')
    #pylab.xticks(range(0,271,180),['0$^{\circ}$','180$^{\circ}$'])
    pylab.xticks([])
    pylab.yticks(range(0,160,100))
    pylab.xlim(-90,270)
    pylab.ylim(0,180)
    pylab.xlabel('$\Delta_{Late}$')
    if dipp[1,1] == 0.:
        pylab.text(-75,150,'$p <$ 10$^{-5}$')
    else:
        pylab.text(-75,150,'%s %.3f' % ('$p = $',dipp[1,1]))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax = pylab.subplot(8,6,35)
    pylab.hist(peakdata[:,2],pylab.arange(-90,271,22.5),color='k')
    #pylab.xticks(range(0,271,180),['0$^{\circ}$','180$^{\circ}$'])
    pylab.xticks([])
    pylab.yticks(range(0,160,100))
    pylab.xlim(-90,270)
    pylab.ylim(0,180)
    pylab.xlabel('$\Delta_{A}$')
    if dipp[2,1] == 0.:
        pylab.text(-75,150,'$p <$ 10$^{-5}$')
    else:
        pylab.text(-75,150,'%s %.3f' % ('$p = $',dipp[2,1]))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax = pylab.subplot(8,6,36)
    pylab.hist(peakdata[:,3],pylab.arange(-90,271,22.5),color='k')
    #pylab.xticks(range(0,271,180),['0$^{\circ}$','180$^{\circ}$'])
    pylab.xticks([])
    pylab.yticks(range(0,160,100))
    pylab.xlim(-90,270)
    pylab.ylim(0,180)
    pylab.xlabel('$\Delta_{B}$')
    if dipp[3,1] == 0.:
        pylab.text(-75,150,'$p <$ 10$^{-5}$')
    else:
        pylab.text(-75,150,'%s %.3f' % ('$p = $',dipp[3,1]))
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
    # F
    n_neu = 1078
    rng = np.random.RandomState(123)
    tuningE0 = rng.rand(n_neu,8)*100
    tuningE1 = rng.rand(n_neu,8)*100
    tuningL0 = rng.rand(n_neu,8)*100
    tuningL1 = rng.rand(n_neu,8)*100
    peakdata = peak_tuning(tuningE0, tuningE1, tuningL0, tuningL1)
    if not os.path.exists(datapath+'dipp_random.txt'):
        # dipp_random.txt is created by MATLAB using dip_pythondata.m
        printdipp = 0
        fpeak = open(datapath+'peak_random.txt','w')
        for j in range(n_neu):
            for k in range(4):
                fpeak.write(str(peakdata[j,k])+'\t')
            fpeak.write('\n')
        fpeak.close()
    else:
        printdipp = 1
        dipp = pylab.loadtxt(datapath+'dipp_random.txt')
    fig.text(0.05,0.25,'F',fontsize=fl_size,ha='right',va='bottom')
    pylab.subplot(4,3,10)
    pylab.scatter(peakdata[:,0],peakdata[:,1],c='k',s=30,lw=0)
    pylab.plot([-90,270],[-90,270],'k--')
    pylab.xlim(-95,275)
    pylab.ylim(-95,275)
    pylab.xticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.yticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.xlabel('$\Delta_{Early}$')
    pylab.ylabel('$\Delta_{Late}$')
    pylab.subplot(4,3,11)
    pylab.scatter(peakdata[:,2],peakdata[:,3],c='k',s=30,lw=0)
    pylab.plot([-90,270],[-90,270],'k--')
    pylab.xlim(-95,275)
    pylab.ylim(-95,275)
    pylab.xticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.yticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.xlabel('$\Delta_{A}$')
    pylab.ylabel('$\Delta_{B}$')
    # G
    fig.text(0.65,0.25,'G',fontsize=fl_size,ha='right',va='bottom')
    ax = pylab.subplot(8,6,41)
    pylab.hist(peakdata[:,0],pylab.arange(-90,271,22.5),color='k')
    #pylab.xticks(range(0,271,180),['0$^{\circ}$','180$^{\circ}$'])
    pylab.xticks([])
    pylab.yticks(range(0,60,50))
    pylab.xlim(-90,270)
    pylab.ylim(0,120)
    pylab.xlabel('$\Delta_{Early}$')
    #pylab.title('%s %.3f %s %.3f' % ('dip=',dipp[0,0],'; p=',dipp[0,1]))
    if printdipp:
        if dipp[0,1] == 0.:
            pylab.text(-60,100,'$p <$ 10$^{-5}$')
        else:
            pylab.text(-60,100,'%s %.3f' % ('$p = $',dipp[0,1]))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax = pylab.subplot(8,6,42)
    pylab.hist(peakdata[:,1],pylab.arange(-90,271,22.5),color='k')
    #pylab.xticks(range(0,271,180),['0$^{\circ}$','180$^{\circ}$'])
    pylab.xticks([])
    pylab.yticks(range(0,60,50))
    pylab.xlim(-90,270)
    pylab.ylim(0,120)
    pylab.xlabel('$\Delta_{Late}$')
    if printdipp:
        if dipp[1,1] == 0.:
            pylab.text(-60,100,'$p <$ 10$^{-5}$')
        else:
            pylab.text(-60,100,'%s %.3f' % ('$p = $',dipp[1,1]))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax = pylab.subplot(8,6,47)
    pylab.hist(peakdata[:,2],pylab.arange(-90,271,22.5),color='k')
    #pylab.xticks(range(0,271,180),['0$^{\circ}$','180$^{\circ}$'])
    pylab.xticks([])
    pylab.yticks(range(0,60,50))
    pylab.xlim(-90,270)
    pylab.ylim(0,120)
    pylab.xlabel('$\Delta_{A}$')
    if printdipp:
        if dipp[2,1] == 0.:
            pylab.text(-60,100,'$p <$ 10$^{-5}$')
        else:
            pylab.text(-60,100,'%s %.3f' % ('$p = $',dipp[2,1]))
    pylab.text(-220,140,'Number',rotation='vertical')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax = pylab.subplot(8,6,48)
    pylab.hist(peakdata[:,3],pylab.arange(-90,271,22.5),color='k')
    #pylab.xticks(range(0,271,180),['0$^{\circ}$','180$^{\circ}$'])
    pylab.xticks([])
    pylab.yticks(range(0,60,50))
    pylab.xlim(-90,270)
    pylab.ylim(0,120)
    pylab.xlabel('$\Delta_{B}$')
    if printdipp:
        if dipp[3,1] == 0.:
            pylab.text(-60,100,'$p <$ 10$^{-5}$')
        else:
            pylab.text(-60,100,'%s %.3f' % ('$p = $',dipp[3,1]))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    eps = 55
    min_samples = 20
    epilson = 30
    tuningE0 = pylab.loadtxt(matlabpath+'tuningA0.txt')
    tuningE1 = pylab.loadtxt(matlabpath+'tuningB0.txt')
    tuningL0 = pylab.loadtxt(matlabpath+'tuningA400.txt')
    tuningL1 = pylab.loadtxt(matlabpath+'tuningB400.txt')
    peakdata = peak_tuning(tuningE0, tuningE1, tuningL0, tuningL1)
    dipp = pylab.loadtxt(matlabpath+'dipp_all.txt')
    # A
    fig.text(0.05,0.9,'A',fontsize=fl_size,ha='right',va='bottom')
    idx,ncount,ncheck,uc,nuc = clustering(peakdata,eps=eps,min_samples=min_samples,epilson=epilson)
    uc.append(0)
    fpeak = open(datapath+'peak_selected_s4.txt','w')
    for j in range(idx.size):
        for k in range(4):
            fpeak.write(str(peakdata[j,k])+'\t')
        fpeak.write('\n')
    fpeak.close()
    fidx = open(datapath+'idx_selected_s4.txt','w')
    for j in range(idx.size):
        if idx[j] == ncheck[0]:
            fidx.write('1\n')
        elif idx[j] == ncheck[1]:
            fidx.write('2\n')
        elif idx[j] == ncheck[2]:
            fidx.write('3\n')
        elif idx[j] == ncheck[3]:
            fidx.write('4\n')
        else:
            fidx.write('0\n')
    fidx.close()
    print 'Number of clusters = '+str(max(idx))
    print 'The number of unclassified clusters is '+str(len(uc))
    print 'The cluster sizes are '+str(ncount)+' , and '+str(nuc)
    print 'Please make sure that the size of all the former clusters is >= 0.05*n and that of all the latter clusters is < 0.05*n'
    pylab.subplot(431)
    pylab.scatter(peakdata[idx==0,0],peakdata[idx==0,1],c='0.5',s=30,lw=0)
    for j in range(len(uc)):
        pylab.scatter(peakdata[idx==uc[j],0],peakdata[idx==uc[j],1],c='0.5',s=30,lw=0)
    if ncheck[0] != 0:
        pylab.scatter(peakdata[idx==ncheck[0],0],peakdata[idx==ncheck[0],1],c=cTG,s=30,lw=0)
    if ncheck[1] != 0:
        pylab.scatter(peakdata[idx==ncheck[1],0],peakdata[idx==ncheck[1],1],c=cCT,s=30,lw=0)
    if ncheck[2] != 0:
        pylab.scatter(peakdata[idx==ncheck[2],0],peakdata[idx==ncheck[2],1],c=cTS1,s=30,lw=0)
    if ncheck[3] != 0:
        pylab.scatter(peakdata[idx==ncheck[3],0],peakdata[idx==ncheck[3],1],c=cTS2,s=30,lw=0)
    pylab.plot([-90,270],[-90,270],'k--')
    pylab.xlim(-95,275)
    pylab.ylim(-95,275)
    pylab.xticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.yticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.xlabel('$\Delta_{Early}$')
    pylab.ylabel('$\Delta_{Late}$')
    pylab.subplot(432)
    pylab.scatter(peakdata[idx==0,2],peakdata[idx==0,3],c='0.5',s=30,lw=0)
    for j in range(len(uc)):
        pylab.scatter(peakdata[idx==uc[j],2],peakdata[idx==uc[j],3],c='0.5',s=30,lw=0)
    if ncheck[0] != 0:
        pylab.scatter(peakdata[idx==ncheck[0],2],peakdata[idx==ncheck[0],3],c=cTG,s=30,lw=0)
    if ncheck[1] != 0:
        pylab.scatter(peakdata[idx==ncheck[1],2],peakdata[idx==ncheck[1],3],c=cCT,s=30,lw=0)
    if ncheck[2] != 0:
        pylab.scatter(peakdata[idx==ncheck[2],2],peakdata[idx==ncheck[2],3],c=cTS1,s=30,lw=0)
    if ncheck[3] != 0:
        pylab.scatter(peakdata[idx==ncheck[3],2],peakdata[idx==ncheck[3],3],c=cTS2,s=30,lw=0)
    pylab.plot([-90,270],[-90,270],'k--')
    pylab.xlim(-95,275)
    pylab.ylim(-95,275)
    pylab.xticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.yticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.xlabel('$\Delta_{A}$')
    pylab.ylabel('$\Delta_{B}$')
    # B
    fig.text(0.05,0.68,'B',fontsize=fl_size,ha='right',va='bottom')
    pylab.subplot(434)
    pylab.scatter([0,180,0],[0,180,180],marker='+',c='k',s=100,lw=2,label='Theory')
    if ncheck[0] != 0:
        pos = cov_ellipse(peakdata[idx==ncheck[0],:2], nstd=1, alpha=alpha, color=cTG)
        pylab.scatter(pos[0],pos[1],marker='x',color='0.5',s=70,lw=2,label='Data')
        pylab.legend(loc=4,scatterpoints=1,frameon=False)
    if ncheck[1] != 0:
        pos = cov_ellipse(peakdata[idx==ncheck[1],:2], nstd=1, alpha=alpha, color=cCT)
        pylab.scatter(pos[0],pos[1],marker='x',color='0.5',s=70,lw=2)
    if ncheck[2] != 0:
        pos = cov_ellipse(peakdata[idx==ncheck[2],:2], nstd=1, alpha=alpha, color=cTS1)
        pylab.scatter(pos[0],pos[1],marker='x',color='0.5',s=70,lw=2)
    if ncheck[3] != 0:
        pos = cov_ellipse(peakdata[idx==ncheck[3],:2], nstd=1, alpha=alpha, color=cTS2)
        pylab.scatter(pos[0],pos[1],marker='x',color='0.5',s=70,lw=2)
    pylab.plot([-90,-23],[-90,-23],'k--')
    pylab.plot([25,150],[25,150],'k--')
    pylab.plot([200,270],[200,270],'k--')
    pylab.xlim(-95,275)
    pylab.ylim(-95,275)
    pylab.xticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.yticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.xlabel('$\Delta_{Early}$')
    pylab.ylabel('$\Delta_{Late}$')
    pylab.subplot(435)
    pylab.scatter([0,0,180],[0,180,0],marker='+',c='k',s=100,lw=2)
    if ncheck[0] != 0:
        pos = cov_ellipse(peakdata[idx==ncheck[0],2:], nstd=1, alpha=alpha, color=cTG)
        pylab.scatter(pos[0],pos[1],marker='x',color='0.5',s=70,lw=2)
    if ncheck[1] != 0:
        pos = cov_ellipse(peakdata[idx==ncheck[1],2:], nstd=1, alpha=alpha, color=cCT)
        pylab.scatter(pos[0],pos[1],marker='x',color='0.5',s=70,lw=2)
    if ncheck[2] != 0:
        pos = cov_ellipse(peakdata[idx==ncheck[2],2:], nstd=1, alpha=alpha, color=cTS1)
        pylab.scatter(pos[0],pos[1],marker='x',color='0.5',s=70,lw=2)
    if ncheck[3] != 0:
        pos = cov_ellipse(peakdata[idx==ncheck[3],2:], nstd=1, alpha=alpha, color=cTS2)
        pylab.scatter(pos[0],pos[1],marker='x',color='0.5',s=70,lw=2)
    pylab.plot([-90,-36],[-90,-36],'k--')
    pylab.plot([40,270],[40,270],'k--')
    pylab.xlim(-95,275)
    pylab.ylim(-95,275)
    pylab.xticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.yticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.xlabel('$\Delta_{A}$')
    pylab.ylabel('$\Delta_{B}$')
    # C
    fig.text(0.65,0.68,'C',fontsize=fl_size,ha='right',va='bottom')
    ax = pylab.subplot(436)
    icount = []
    for j in range(int(max(idx))+1):
        icount.append(sum(idx==j))
    pylab.bar(range(4),ncount/1082.,color=[cTG,cCT,cTS1,cTS2])
    pylab.xticks(pylab.arange(0,4,1),('TG','CT','TS1','TS2'))
    pylab.xlim(-0.6,3.6)
    pylab.ylabel('Fraction')
    #pylab.ylim(0,150)
    #pylab.yticks(range(0,151,50))
    pylab.ylim(0,0.13)
    pylab.yticks([0,0.05,0.1],('0','0.05','0.1'))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    pylab.subplots_adjust(bottom=0.1,hspace=0.35,wspace=0.4)
    pylab.savefig(figpath+'s2.pdf')

def supp4():
    print 'Making supplementary figure 4...'
    fig = pylab.figure(figsize=[7.5,11])
    selected_peakdata = pylab.loadtxt(datapath+'peak_selected_s4.txt')
    selected_idx = pylab.loadtxt(datapath+'idx_selected_s4.txt')
    selectedID = pylab.loadtxt(matlabpath+'selectedID.txt',dtype='int')
    with open(matlabpath+'brainarea_all.txt','r') as tempfile:
        ba = tempfile.read().replace('\n', '')
    numba = []
    for j in range(len(ba)):
        if ba[j] == 'v':
            numba.append(0)
        else:
            numba.append(1)
    numba = np.array(numba)
    nA = sum(numba==0)
    nB = sum(numba==1)
    numba = np.take(numba,selectedID-1)
    print sum(numba==0),sum(numba==1)
    ncount = np.zeros((2,4))
    # A
    peakdata = selected_peakdata[numba==0,:]
    idx = selected_idx[numba==0]
    for j in range(4):
        ncount[0,j] = sum(idx==j+1)
    fig.text(0.05,0.9,'A',fontsize=fl_size,ha='right',va='bottom')
    pylab.subplot(321)
    pylab.scatter(peakdata[idx==0,0],peakdata[idx==0,1],c='0.5',s=30,lw=0)
    pylab.scatter(peakdata[idx==1,0],peakdata[idx==1,1],c=cTG,s=30,lw=0)
    pylab.scatter(peakdata[idx==2,0],peakdata[idx==2,1],c=cCT,s=30,lw=0)
    pylab.scatter(peakdata[idx==3,0],peakdata[idx==3,1],c=cTS1,s=30,lw=0)
    pylab.scatter(peakdata[idx==4,0],peakdata[idx==4,1],c=cTS2,s=30,lw=0)
    pylab.plot([-90,270],[-90,270],'k--')
    pylab.xlim(-95,275)
    pylab.ylim(-95,275)
    pylab.xticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.yticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.xlabel('$\Delta_{Early}$')
    pylab.ylabel('$\Delta_{Late}$')
    pylab.subplot(322)
    pylab.scatter(peakdata[idx==0,2],peakdata[idx==0,3],c='0.5',s=30,lw=0)
    pylab.scatter(peakdata[idx==1,2],peakdata[idx==1,3],c=cTG,s=30,lw=0)
    pylab.scatter(peakdata[idx==2,2],peakdata[idx==2,3],c=cCT,s=30,lw=0)
    pylab.scatter(peakdata[idx==3,2],peakdata[idx==3,3],c=cTS1,s=30,lw=0)
    pylab.scatter(peakdata[idx==4,2],peakdata[idx==4,3],c=cTS2,s=30,lw=0)
    pylab.plot([-90,270],[-90,270],'k--')
    pylab.xlim(-95,275)
    pylab.ylim(-95,275)
    pylab.xticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.yticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.xlabel('$\Delta_{A}$')
    pylab.ylabel('$\Delta_{B}$')
    # B
    fig.text(0.05,0.62,'B',fontsize=fl_size,ha='right',va='bottom')
    peakdata = selected_peakdata[numba==1,:]
    idx = selected_idx[numba==1]
    idx0 = np.ones(idx.size)
    for j in range(4):
        ncount[1,j] = sum(idx==j+1)
    pylab.subplot(323)
    pylab.scatter(peakdata[idx==0,0],peakdata[idx==0,1],c='0.5',s=30,lw=0)
    pylab.scatter(peakdata[idx==1,0],peakdata[idx==1,1],c=cTG,s=30,lw=0)
    pylab.scatter(peakdata[idx==2,0],peakdata[idx==2,1],c=cCT,s=30,lw=0)
    pylab.scatter(peakdata[idx==3,0],peakdata[idx==3,1],c=cTS1,s=30,lw=0)
    pylab.scatter(peakdata[idx==4,0],peakdata[idx==4,1],c=cTS2,s=30,lw=0)
    pylab.plot([-90,270],[-90,270],'k--')
    pylab.xlim(-95,275)
    pylab.ylim(-95,275)
    pylab.xticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.yticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.xlabel('$\Delta_{Early}$')
    pylab.ylabel('$\Delta_{Late}$')
    pylab.subplot(324)
    pylab.scatter(peakdata[idx==0,2],peakdata[idx==0,3],c='0.5',s=30,lw=0)
    pylab.scatter(peakdata[idx==1,2],peakdata[idx==1,3],c=cTG,s=30,lw=0)
    pylab.scatter(peakdata[idx==2,2],peakdata[idx==2,3],c=cCT,s=30,lw=0)
    pylab.scatter(peakdata[idx==3,2],peakdata[idx==3,3],c=cTS1,s=30,lw=0)
    pylab.scatter(peakdata[idx==4,2],peakdata[idx==4,3],c=cTS2,s=30,lw=0)
    pylab.plot([-90,270],[-90,270],'k--')
    pylab.xlim(-95,275)
    pylab.ylim(-95,275)
    pylab.xticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.yticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.xlabel('$\Delta_{A}$')
    pylab.ylabel('$\Delta_{B}$')
    # C
    fig.text(0.05,0.31,'C',fontsize=fl_size,ha='right',va='bottom')
    #ax = pylab.subplot(313)
    ax = fig.add_axes([0.18, 0.1, 0.68, 0.22])
    pylab.bar([0,1,2.5,3.5,5,6],[ncount[0,0]/float(nA),ncount[1,0]/float(nB),ncount[0,1]/float(nA),ncount[1,1]/float(nB),ncount[0,2]/float(nA),ncount[1,2]/float(nB)],color=[cTG,cTG,cCT,cCT,cTS1,cTS1])
    pylab.bar([5,6],[ncount[0,3]/float(nA),ncount[1,3]/float(nB)],bottom=[ncount[0,2]/float(nA),ncount[1,2]/float(nB)],color=[cTS2,cTS2])
    pval = compp(ncount[0,0]/float(nA),ncount[1,0]/float(nB),nA,nB)
    print nA, nB
    print 'TG: '+str(pval)
    if pval < 0.05:
        txt = '*'
    else:
        txt = 'NS'
    sigdiffh = 1.1*max(ncount[0,0]/float(nA),ncount[1,0]/float(nB))
    pylab.plot([0,0,1,1],[sigdiffh,sigdiffh+0.02,sigdiffh+0.02,sigdiffh],'k',lw=1.5)
    pylab.text(0.5,sigdiffh+0.02,txt,ha='center',va='bottom')
    pval = compp(ncount[0,1]/float(nA),ncount[1,1]/float(nB),nA,nB)
    print 'CT: '+str(pval)
    if pval < 0.05:
        txt = '*'
    else:
        txt = 'NS'
    sigdiffh = 1.1*max(ncount[0,1]/float(nA),ncount[1,1]/float(nB))
    pylab.plot([2.5,2.5,3.5,3.5],[sigdiffh,sigdiffh+0.02,sigdiffh+0.02,sigdiffh],'k',lw=1.5)
    pylab.text(3,sigdiffh+0.02,txt,ha='center',va='bottom')
    pval = compp(ncount[0,2]/float(nA)+ncount[0,3]/float(nA),ncount[1,2]/float(nB)+ncount[1,3]/float(nB),nA,nB)
    print 'TS: '+str(pval)
    if pval < 0.05:
        txt = '*'
    else:
        txt = 'NS'
    sigdiffh = 1.1*max(ncount[0,2]/float(nA)+ncount[0,3]/float(nA),ncount[1,2]/float(nB)+ncount[1,3]/float(nB))
    pylab.plot([5,5,6,6],[sigdiffh,sigdiffh+0.02,sigdiffh+0.02,sigdiffh],'k',lw=1.5)
    pylab.text(5.5,sigdiffh+0.02,txt,ha='center',va='bottom')
    pylab.xticks([0,1,2.5,3.5,5,6],('v','d','v','d','v','d'))
    pylab.xlabel('    TG                             CT                             TS')
    #pylab.xticks([0.5,3.5,6.5],('\nTG','\nCT','\nTS'))
    pylab.xlim(-0.9,6.6)
    pylab.ylabel('Fraction')
    pylab.yticks([0,0.2,0.4])
    pylab.ylim(0,0.25)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    pylab.subplots_adjust(bottom=0.1,hspace=0.35,wspace=0.4)
    pylab.savefig(figpath+'s4.pdf')

def supp5():
    print 'Making supplementary figure 5...'
    J_p,J_m,g_wt,g_tt,wmIno0,ctIno0,std_J,std_t,std_I,std_sig = para_heterogeneous()
    fname = 'jp'+str(J_p)+'jm'+str(-J_m)+'gwt'+str(g_wt)+'gtt'+str(g_tt)+'I'+str(wmIno0)+'std'+str(std_J)+'seednet'+str(seednet)+'seedrun'+str(seedrun)
    pylab.figure(figsize=[12,8])
    item1 = pylab.loadtxt(datapath+fname+'/lpfc_t1_net_'+str(seednet)+'_run_'+str(seedrun)+'_CJ_0_d_90_.txt')
    item2 = pylab.loadtxt(datapath+fname+'/lpfc_t2_net_'+str(seednet)+'_run_'+str(seedrun)+'_CJ_0_d_90_.txt')
    t2p = [10,100,200,300,400,500]
    for j in range(6):
        ax = pylab.subplot(2,3,j+1)
        pylab.plot(alldir[:-1],1000*item1[t2p[j]/tp,:],label='IN-A')
        pylab.plot(alldir[:-1],1000*item2[t2p[j]/tp,:],label='IN-B')
        pylab.title('$t = $' + str(t2p[j])+'ms')
        if j == 0:
            pylab.legend(loc=1,frameon=False)
        if j > 2:
            pylab.xlabel('Neuron')
        if j < 3:
            pylab.xticks([])
        else:
            pylab.xticks(range(0,361,180),['0$^{\circ}$','180$^{\circ}$','360$^{\circ}$'])
        pylab.xlim(0,361)
        pylab.ylim(0,62)
        pylab.yticks([0,30,60])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
    pylab.savefig(figpath+'s5.pdf')

def supp6():
    print 'Making supplementary figure 6...'
    fig = pylab.figure(figsize=[13,4])
    # A
    mpl.rcParams['font.size'] = font_size-3
    fig.text(0.05,0.9,'A',fontsize=fl_size,ha='right',va='bottom')
    peakdata = pylab.loadtxt(datapath+'peak_sim.txt')
    idx = pylab.loadtxt(datapath+'idx_sim.txt')
    ncheck = pylab.loadtxt(datapath+'ncheck_lenuc_sim.txt')
    ax = pylab.subplot(131)
    tempidx = (idx!=ncheck[0])*(idx!=ncheck[1])*(idx!=ncheck[2])*(idx!=ncheck[3])
    pylab.scatter(peakdata[tempidx,0],peakdata[tempidx,1],color='0.5',s=30,lw=0)
    if ncheck[0] != 0:
        pylab.scatter(peakdata[idx==ncheck[0],0],peakdata[idx==ncheck[0],1],c=cTG,s=30,lw=0)
    if ncheck[1] != 0:
        pylab.scatter(peakdata[idx==ncheck[1],0],peakdata[idx==ncheck[1],1],c=cCT,s=30,lw=0)
    if ncheck[2] != 0:
        pylab.scatter(peakdata[idx==ncheck[2],0],peakdata[idx==ncheck[2],1],c=cTS1,s=30,lw=0)
    if ncheck[3] != 0:
        pylab.scatter(peakdata[idx==ncheck[3],0],peakdata[idx==ncheck[3],1],c=cTS2,s=30,lw=0)
    pylab.plot([-90,270],[-90,270],'k--')
    pylab.xlim(-95,275)
    pylab.ylim(-95,275)
    pylab.xticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.yticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.xlabel('$\Delta_{Early}$')
    pylab.ylabel('$\Delta_{Late}$')
    ax = pylab.subplot(132)
    pylab.scatter(peakdata[tempidx,2],peakdata[tempidx,3],color='0.5',s=30,lw=0)
    if ncheck[0] != 0:
        pylab.scatter(peakdata[idx==ncheck[0],2],peakdata[idx==ncheck[0],3],c=cTG,s=30,lw=0)
    if ncheck[1] != 0:
        pylab.scatter(peakdata[idx==ncheck[1],2],peakdata[idx==ncheck[1],3],c=cCT,s=30,lw=0)
    if ncheck[2] != 0:
        pylab.scatter(peakdata[idx==ncheck[2],2],peakdata[idx==ncheck[2],3],c=cTS1,s=30,lw=0)
    if ncheck[3] != 0:
        pylab.scatter(peakdata[idx==ncheck[3],2],peakdata[idx==ncheck[3],3],c=cTS2,s=30,lw=0)
    pylab.plot([-90,270],[-90,270],'k--')
    pylab.xlim(-95,275)
    pylab.ylim(-95,275)
    pylab.xticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.yticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.xlabel('$\Delta_{A}$')
    pylab.ylabel('$\Delta_{B}$')
    # B
    fig.text(0.65,0.9,'B',fontsize=fl_size,ha='right',va='bottom')
    if not os.path.exists(datapath+'dipp_sim.txt'):
        printdipp = 0
        fpeak = open(datapath+'peak_sim.txt','w')
        for j in range(n*3):
            for k in range(4):
                fpeak.write(str(peakdata[j,k])+'\t')
            fpeak.write('\n')
        fpeak.close()
    else:
        printdipp = 1
        dipp = pylab.loadtxt(datapath+'dipp_sim.txt')
    ax = pylab.subplot(2,6,5)
    pylab.hist(peakdata[:,0],pylab.arange(-90,271,22.5),color='k')
    pylab.xticks([])
    pylab.yticks(range(0,201,100))
    pylab.xlim(-90,270)
    pylab.ylim(0,300)
    pylab.xlabel('$\Delta_{Early}$')
    if printdipp:
        if dipp[0,1] == 0.:
            pylab.text(45,200,'$p <$ 10$^{-5}$')
        else:
            pylab.text(50,200,'%s %.3f' % ('$p = $',dipp[0,1]))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax = pylab.subplot(2,6,6)
    pylab.hist(peakdata[:,1],pylab.arange(-90,271,22.5),color='k')
    #pylab.xticks(range(0,271,180),['0$^{\circ}$','180$^{\circ}$'])
    pylab.xticks([])
    pylab.yticks(range(0,201,100))
    pylab.xlim(-90,270)
    pylab.ylim(0,300)
    pylab.xlabel('$\Delta_{Late}$')
    if printdipp:
        if dipp[1,1] == 0.:
            pylab.text(-70,200,'$p <$ 10$^{-5}$')
        else:
            pylab.text(-70,200,'%s %.3f' % ('$p = $',dipp[1,1]))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax = pylab.subplot(2,6,11)
    pylab.hist(peakdata[:,2],pylab.arange(-90,271,22.5),color='k')
    pylab.xticks(range(0,271,180),['0$^{\circ}$','180$^{\circ}$'])
    pylab.yticks(range(0,201,100))
    pylab.xlim(-90,270)
    pylab.ylim(0,300)
    pylab.xlabel('$\Delta_{A}$')
    if printdipp:
        if dipp[2,1] == 0.:
            pylab.text(45,200,'$p <$ 10$^{-5}$')
        else:
            pylab.text(50,200,'%s %.3f' % ('$p = $',dipp[2,1]))
    mpl.rcParams['font.size'] = font_size
    pylab.text(-240,420,'Number',rotation='vertical')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    mpl.rcParams['font.size'] = font_size-3
    ax = pylab.subplot(2,6,12)
    pylab.hist(peakdata[:,3],pylab.arange(-90,271,22.5),color='k')
    pylab.xticks(range(0,271,180),['0$^{\circ}$','180$^{\circ}$'])
    pylab.yticks(range(0,201,100))
    pylab.xlim(-90,270)
    pylab.ylim(0,300)
    pylab.xlabel('$\Delta_{B}$')
    if printdipp:
        if dipp[3,1] == 0.:
            pylab.text(45,200,'$p <$ 10$^{-5}$')
        else:
            pylab.text(50,200,'%s %.3f' % ('$p = $',dipp[3,1]))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    pylab.subplots_adjust(hspace=0.45,wspace=0.45,bottom=0.2)
    pylab.savefig(figpath+'s6.pdf')

def supp7():
    print 'Making supplementary figure 7...'
    fig = pylab.figure(figsize=[12,3.5])
    # A
    fig.text(0.05,0.9,'A',fontsize=fl_size,ha='right',va='bottom')
    ax = pylab.subplot(131)
    J_p,J_m,g_wt,g_tt,wmIno0,ctIno0,std_J,std_t,std_I,std_sig = para_heterogeneous_neuron()
    fname = 'jp'+str(J_p)+'jm'+str(-J_m)+'gwt'+str(g_wt)+'gtt'+str(g_tt)+'I'+str(wmIno0)+'std'+str(std_J)+'stdI'+str(std_I)+'stds'+str(std_sig)+'seednet'+str(seednet)+'seedrun'+str(seedrun)
    tuningE0, tuningE1, tuningL0, tuningL1 = tuning_curves(path=datapath+fname)
    peakdata = peak_tuning(tuningE0,tuningE1,tuningL0,tuningL1)
    idx,ncount,ncheck,uc,nuc = clustering(peakdata,eps=eps_sim)
    uc.append(0)
    fpeak = open(datapath+'peak_selected.txt','w')
    for j in range(idx.size):
        for k in range(4):
            fpeak.write(str(peakdata[j,k])+'\t')
        fpeak.write('\n')
    fpeak.close()
    for j in range(len(uc)):
        pylab.scatter(peakdata[idx==uc[j],0],peakdata[idx==uc[j],1],c='0.5',s=30,lw=0)
    if ncheck[0] != 0:
        pylab.scatter(peakdata[idx==ncheck[0],0],peakdata[idx==ncheck[0],1],c=cTG,s=30,lw=0)
    if ncheck[1] != 0:
        pylab.scatter(peakdata[idx==ncheck[1],0],peakdata[idx==ncheck[1],1],c=cCT,s=30,lw=0)
    if ncheck[2] != 0:
        pylab.scatter(peakdata[idx==ncheck[2],0],peakdata[idx==ncheck[2],1],c=cTS1,s=30,lw=0)
    if ncheck[3] != 0:
        pylab.scatter(peakdata[idx==ncheck[3],0],peakdata[idx==ncheck[3],1],c=cTS2,s=30,lw=0)
    pylab.plot([-90,270],[-90,270],'k--')
    pylab.xlim(-95,275)
    pylab.ylim(-95,275)
    pylab.xticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.yticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.xlabel('$\Delta_{Early}$')
    pylab.ylabel('$\Delta_{Late}$')
    pylab.subplot(132)
    pylab.scatter(peakdata[idx==0,2],peakdata[idx==0,3],c='0.5',s=30,lw=0)
    for j in range(len(uc)):
        pylab.scatter(peakdata[idx==uc[j],2],peakdata[idx==uc[j],3],c='0.5',s=30,lw=0)
    if ncheck[0] != 0:
        pylab.scatter(peakdata[idx==ncheck[0],2],peakdata[idx==ncheck[0],3],c=cTG,s=30,lw=0)
    if ncheck[1] != 0:
        pylab.scatter(peakdata[idx==ncheck[1],2],peakdata[idx==ncheck[1],3],c=cCT,s=30,lw=0)
    if ncheck[2] != 0:
        pylab.scatter(peakdata[idx==ncheck[2],2],peakdata[idx==ncheck[2],3],c=cTS1,s=30,lw=0)
    if ncheck[3] != 0:
        pylab.scatter(peakdata[idx==ncheck[3],2],peakdata[idx==ncheck[3],3],c=cTS2,s=30,lw=0)
    pylab.plot([-90,270],[-90,270],'k--')
    pylab.xlim(-95,275)
    pylab.ylim(-95,275)
    pylab.xticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.yticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.xlabel('$\Delta_{A}$')
    pylab.ylabel('$\Delta_{B}$')
    # B
    fig.text(0.68,0.9,'B',fontsize=fl_size,ha='right',va='bottom')
    ax = pylab.subplot(133)
    pylab.bar(range(4),ncount/(3.*n),color=[cTG,cCT,cTS1,cTS2])
    pylab.xticks(pylab.arange(0,4,1),('TG','CT','TS1','TS2'))
    pylab.xlim(-0.6,3.6)
    pylab.ylabel('Fraction')
    #pylab.ylim(0,150)
    pylab.yticks(np.arange(0,0.35,0.1))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    pylab.subplots_adjust(bottom=0.15,wspace=0.4,right=0.95)
    pylab.savefig(figpath+'s7.pdf')

#supp1()
#supp2()
#supp4()
#supp5()
#supp6()

J_p = 2.5
J_m = -1.6
g_wt = 0.05
g_st = 0.243
g_tt = 0.3
bamp = 0.9
wmIno0 = 0.31
ctIno0 = 0.31
std_J = 2.0
ddir = 180
fname = 'input2pk_jp'+str(J_p)+'jm'+str(-J_m)+'gwt'+str(g_wt)+'gst'+str(g_st)+'gtt'+str(g_tt)+'b'+str(bamp)+'I'+str(wmIno0)+'std'+str(std_J)+'seednet'+str(seednet)+'seedrun'+str(seedrun)
if not os.path.exists(datapath+fname):
    JM1,JM2,JMct,JMx1,JMx2 = define_connect()
    if not os.path.exists(datapath+fname):
        os.makedirs(datapath+fname)
    print 'Running simulations for '+fname
    for CJ in range(2):
        for dir in range(0,360,45):
            print 'dir: '+str(dir)
            trial(datapath+fname,dir,CJ,npeak=2,bamp=bamp,opp=pylab.mod(dir+ddir,360))

def supp8():
    print 'Making supplementary figure 8...'
    fig = pylab.figure(figsize=[12,15])
    pylab.subplots_adjust(top=0.9,bottom=0.1,hspace=0.35,wspace=0.35)
    # A
    fig.text(0.05,0.9,'A',fontsize=fl_size,ha='right',va='bottom')
    ax = pylab.subplot(541)
    profile = visgauss(90,270,sigma_s,npeak=2,n=n,bamp=0.9)
    pylab.plot(alldir,np.append(profile[0,:],profile[0,0]),cA,linewidth=2)
    pylab.plot(alldir,np.append(profile[1,:],profile[1,0]),cB,linewidth=2)
    pylab.plot([90,270],[1,1],'k')
    pylab.plot([90],[1],'k<')
    pylab.plot([265],[1],'k>')
    pylab.text(140,1.03,'180$^{\circ}$')
    pylab.xlim(0,360)
    pylab.ylim(0,1.3)
    pylab.xticks([])
    pylab.yticks([])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    #pylab.xlabel('TG')
    pylab.title('Visual input')
    pylab.xlabel('Location')
    # B
    fig.text(0.33,0.9,'B',fontsize=fl_size,ha='right',va='bottom')
    pylab.subplot(13,3,2)
    fname = 'input2pk_jp2.5jm1.6gwt0.05gst0.243gtt0.3b0.9I0.31std2.0seednet101seedrun51'
    ymax = 60
    #fname = 'input2pk_jp2.1jm2.65gwt0.13gst0.1gtt1.0b0.9I0.31std2.0seednet101seedrun51'
    #ymax = 12
    print fname
    item = pylab.loadtxt(datapath+fname+'/lpfc_t1_net_'+str(seednet)+'_run_'+str(seedrun)+'_CJ_0_d_90_.txt')
    pylab.pcolor(item.T*1000,rasterized=True,cmap='jet')
    pylab.clim(0,ymax)
    cb = pylab.colorbar(ticks=[])
    cb.remove()
    pylab.ylabel('IN-A')
    pylab.yticks(range(0,n+1,n/2),['0$^{\circ}$','180$^{\circ}$','360$^{\circ}$'])
    pylab.xticks([])
    pylab.xlim(0,targeton/tp)
    pylab.ylim(0,n)
    pylab.title('A chosen')
    pylab.subplot(13,3,5)
    item = pylab.loadtxt(datapath+fname+'/lpfc_t2_net_'+str(seednet)+'_run_'+str(seedrun)+'_CJ_0_d_90_.txt')
    pylab.pcolor(item.T*1000,rasterized=True,cmap='jet')
    pylab.clim(0,ymax)
    cb = pylab.colorbar(ticks=[])
    cb.remove()
    pylab.ylabel('IN-B')
    pylab.yticks(range(0,n+1,n/2),['0$^{\circ}$','180$^{\circ}$','360$^{\circ}$'])
    pylab.xticks([])
    pylab.xlim(0,targeton/tp)
    pylab.ylim(0,n)
    pylab.subplot(13,3,8)
    item = pylab.loadtxt(datapath+fname+'/lpfc_ct_net_'+str(seednet)+'_run_'+str(seedrun)+'_CJ_0_d_90_.txt')
    pylab.pcolor(item.T*1000,rasterized=True,cmap='jet')
    pylab.clim(0,ymax)
    cb = pylab.colorbar(ticks=[])
    cb.remove()
    pylab.ylabel('RO')
    pylab.yticks(range(0,n+1,n/2),['0$^{\circ}$','180$^{\circ}$','360$^{\circ}$'])
    pylab.xticks(pylab.arange(0,(dur-tgt_start-go)/tp+1,500/tp),['0','500','1000'])
    pylab.xlim(0,targeton/tp)
    pylab.ylim(0,n)
    pylab.xlabel('Time from target onset (ms)')
    # C
    fig.text(0.63,0.9,'C',fontsize=fl_size,ha='right',va='bottom')
    pylab.subplot(13,3,3)
    item = pylab.loadtxt(datapath+fname+'/lpfc_t1_net_'+str(seednet)+'_run_'+str(seedrun)+'_CJ_1_d_90_.txt')
    pylab.pcolor(item.T*1000,rasterized=True,cmap='jet')
    pylab.clim(0,ymax)
    cbar = pylab.colorbar(ticks=[0,ymax/2,ymax])
    cbar.set_label('Rate (Hz)')
    pylab.ylabel('IN-A')
    pylab.yticks(range(0,n+1,n/2),['0$^{\circ}$','180$^{\circ}$','360$^{\circ}$'])
    pylab.xticks([])
    pylab.yticks([])
    pylab.xlim(0,targeton/tp)
    pylab.ylim(0,n)
    pylab.title('B chosen')
    pylab.subplot(13,3,6)
    item = pylab.loadtxt(datapath+fname+'/lpfc_t2_net_'+str(seednet)+'_run_'+str(seedrun)+'_CJ_1_d_90_.txt')
    pylab.pcolor(item.T*1000,rasterized=True,cmap='jet')
    pylab.clim(0,ymax)
    cb = pylab.colorbar(ticks=[])
    cb.remove()
    #cbar.set_label('Rate (Hz)')
    pylab.ylabel('IN-B')
    pylab.yticks(range(0,n+1,n/2),['0$^{\circ}$','180$^{\circ}$','360$^{\circ}$'])
    pylab.xticks([])
    pylab.yticks([])
    pylab.xlim(0,targeton/tp)
    pylab.ylim(0,n)
    pylab.subplot(13,3,9)
    item = pylab.loadtxt(datapath+fname+'/lpfc_ct_net_'+str(seednet)+'_run_'+str(seedrun)+'_CJ_1_d_90_.txt')
    pylab.pcolor(item.T*1000,rasterized=True,cmap='jet')
    pylab.clim(0,ymax)
    cb = pylab.colorbar(ticks=[])
    cb.remove()
    pylab.yticks([])
    #pylab.yticks(range(0,n+1,n/2),['0$^{\circ}$','180$^{\circ}$','360$^{\circ}$'])
    pylab.xticks(pylab.arange(0,(dur-tgt_start-go)/tp+1,500/tp),['0','500','1000'])
    pylab.xlim(0,targeton/tp)
    pylab.ylim(0,n)
    pylab.ylabel('RO')
    pylab.xlabel('Time from target onset (ms)')
    # D
    selectedID = [62] # 53, arbitrary, random, to be reset: 4
    tune0 = []
    tune1 = []
    peak0 = pylab.empty(11)
    peak1 = pylab.empty(11)
    fig.text(0.05,0.67,'D',fontsize=fl_size,ha='right',va='bottom')
    for j in range(3):
        tuning0, tuning1 = tuning_cj(datapath+fname,200*j)
        tune0.append(tuning0[selectedID,:]*1000)
        tune1.append(tuning1[selectedID,:]*1000)
    for k in range(1):
        for j in range(11):
            tuning0, tuning1 = tuning_cj(datapath+fname,50*j)
            peak0[j] = prefdir(tuning0[selectedID[k],:])
            peak1[j] = prefdir(tuning1[selectedID[k],:])
        ax = pylab.subplot(7,4,k*4+12)
        pylab.plot(pylab.arange(100,601,50),peak0,c=cA,lw=2,label='A chosen')
        pylab.plot(pylab.arange(100,601,50),peak1,c=cB,lw=2,label='B chosen')
        #if k == 0:
            #pylab.legend(loc=1,frameon=False)
        pylab.xlim(0,700)
        pylab.ylim(0,363)
        pylab.xticks(range(0,601,300))
        pylab.yticks(range(0,361,180),['0$^{\circ}$','180$^{\circ}$','360$^{\circ}$'])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        pylab.title('Peak')
        pylab.xlabel('Time from target onset (ms)')
        for j in range(3):
            ax = pylab.subplot(7,4,k*4+9+j)
            fc1 = sp.interp1d(range(-180,360+180,45),np.append(tune0[j][k][-4:],np.append(tune0[j][k],tune0[j][k][:4])),kind='cubic')
            fc2 = sp.interp1d(range(-180,360+180,45),np.append(tune1[j][k][-4:],np.append(tune1[j][k],tune1[j][k][:4])),kind='cubic')
            pylab.plot(range(0,360,45),tune0[j][k],'o',c=cA)
            pylab.plot(range(0,360,15),fc1(range(0,361,15))[:-1],c=cA,lw=2)
            pylab.plot(range(0,360,45),tune1[j][k],'o',c=cB)
            pylab.plot(range(0,360,15),fc2(range(0,361,15))[:-1],c=cB,lw=2)
            pylab.plot(prefdir(tune0[j][k])*pylab.ones(2),[-20,120],c=cA)
            pylab.plot(prefdir(tune1[j][k])*pylab.ones(2),[-20,120],c=cB)
            pylab.xticks(range(0,360,90),['0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
            pylab.xlim(-20,335)
            pylab.yticks(range(0,100,20))
            pylab.ylim(-ymax/12,ymax)
            #pylab.ylim(-1,10)
            pylab.xlabel('Target A')
            if j == 0:
                pylab.ylabel('Rate (Hz)')
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')
            if k == 0:
                if j == 0:
                    pylab.title('Early')
                elif j == 1:
                    pylab.title('Mid')
                elif j == 2:
                    pylab.title('Late')
    # E
    fig.text(0.05,0.48,'E',fontsize=fl_size,ha='right',va='bottom')
    tuningE0, tuningE1, tuningL0, tuningL1 = tuning_curves(path=datapath+fname)
    print fname
    peakdata = peak_tuning(tuningE0, tuningE1, tuningL0, tuningL1)
    ### Test identities
    idx,ncount,ncheck,uc,nuc = clustering(peakdata[:n,:],eps=eps_sim)
    print ncount,uc
    idx,ncount,ncheck,uc,nuc = clustering(peakdata[n:2*n,:],eps=eps_sim)
    print ncount,uc
    idx,ncount,ncheck,uc,nuc = clustering(peakdata[-n:,:],eps=eps_sim)
    print ncount,uc
    idx,ncount,ncheck,uc,nuc = clustering(peakdata,eps=eps_sim)
    uc.append(0)
    ax = pylab.subplot(437)
    pylab.plot([-90,270],[-90,270],'k--')
    for j in range(len(uc)):
        pylab.scatter(peakdata[idx==uc[j],0],peakdata[idx==uc[j],1],c='0.5',s=30,lw=0)
        if ncheck[0] != 0:
            pylab.scatter(peakdata[idx==ncheck[0],0],peakdata[idx==ncheck[0],1],c=cTG,s=30,lw=0)
    if ncheck[1] != 0:
        pylab.scatter(peakdata[idx==ncheck[1],0],peakdata[idx==ncheck[1],1],c=cCT,s=30,lw=0)
        if ncheck[2] != 0:
            pylab.scatter(peakdata[idx==ncheck[2],0],peakdata[idx==ncheck[2],1],c=cTS1,s=30,lw=0)
    if ncheck[3] != 0:
        pylab.scatter(peakdata[idx==ncheck[3],0],peakdata[idx==ncheck[3],1],c=cTS2,s=30,lw=0)
        pylab.xlim(-95,275)
        pylab.ylim(-95,275)
        pylab.xticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
        pylab.yticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
        pylab.xlabel('$\Delta_{Early}$')
        pylab.ylabel('$\Delta_{Late}$')
    pylab.subplot(438)
    pylab.plot([-90,270],[-90,270],'k--')
    pylab.scatter(peakdata[idx==0,2],peakdata[idx==0,3],c='0.5',s=30,lw=0)
    for j in range(len(uc)):
        pylab.scatter(peakdata[idx==uc[j],2],peakdata[idx==uc[j],3],c='0.5',s=30,lw=0)
    if ncheck[0] != 0:
        pylab.scatter(peakdata[idx==ncheck[0],2],peakdata[idx==ncheck[0],3],c=cTG,s=30,lw=0)
    if ncheck[1] != 0:
        pylab.scatter(peakdata[idx==ncheck[1],2],peakdata[idx==ncheck[1],3],c=cCT,s=30,lw=0)
    if ncheck[2] != 0:
        pylab.scatter(peakdata[idx==ncheck[2],2],peakdata[idx==ncheck[2],3],c=cTS1,s=30,lw=0)
    if ncheck[3] != 0:
        pylab.scatter(peakdata[idx==ncheck[3],2],peakdata[idx==ncheck[3],3],c=cTS2,s=30,lw=0)
    pylab.xlim(-95,275)
    pylab.ylim(-95,275)
    pylab.xticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.yticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.xlabel('$\Delta_{A}$')
    pylab.ylabel('$\Delta_{B}$')
    # F
    fig.text(0.65,0.48,'F',fontsize=fl_size,ha='right',va='bottom')
    ax = pylab.subplot(439)
    icount = np.zeros(4)
    for j in range(4):
        icount[j] = sum(idx==ncheck[j])
    print icount
    pylab.bar(range(4),icount/(3.*n),color=[cTG,cCT,cTS1,cTS2])
    pylab.xticks(pylab.arange(0,4,1),('TG','CT','TS1','TS2'))
    pylab.xlim(-0.6,3.6)
    pylab.ylabel('Fraction')
    pylab.ylim(0,0.45)
    pylab.yticks(pylab.arange(0,0.41,0.2))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    # G: Data
    fig.text(0.05,0.25,'G',fontsize=fl_size,ha='right',va='bottom')
    fig.text(0.32,0.25,'LPFC',fontsize=fl_size,ha='right',va='bottom')
    idx_LPFC = pylab.loadtxt(datapath+'idx_selected.txt')
    ncheck_LPFC = [1,2,3,4]
    tuneA0 = pylab.loadtxt(matlabpath+'tuningA0.txt')
    tuneB0 = pylab.loadtxt(matlabpath+'tuningB0.txt')
    tuneA400 = pylab.loadtxt(matlabpath+'tuningA400.txt')
    tuneB400 = pylab.loadtxt(matlabpath+'tuningB400.txt')
    if 1:
        print 'Baseline subtracted'
        for j in range(tuneA0.shape[0]):
            tuneA0[j] = tuneA0[j] - np.min(tuneA0[j])
            tuneB0[j] = tuneB0[j] - np.min(tuneB0[j])
            tuneA400[j] = tuneA400[j] - np.min(tuneA400[j])
            tuneB400[j] = tuneB400[j] - np.min(tuneB400[j])
    peakTG1 = [[],[]]    # Early,Late: Main peak
    peakTG2 = [[],[]]    # Opp to peak
    peakCT1 = [[],[]]
    peakCT2 = [[],[]]
    peakTS1 = [[],[]]
    peakTS2 = [[],[]]
    for j in range(len(idx_LPFC)):
        if idx_LPFC[j] == ncheck_LPFC[0]: # 1 TG 2 CT 3 4 TS
            peakTG1[0].append(np.max(tuneA0[j]))
            peakTG1[0].append(np.max(tuneB0[j]))
            peakTG2[0].append(tuneA0[j,np.mod(np.argmax(tuneA0[j])+4,8)])
            peakTG2[0].append(tuneB0[j,np.mod(np.argmax(tuneB0[j])+4,8)])
            peakTG1[1].append(np.max(tuneA400[j]))
            peakTG1[1].append(np.max(tuneB400[j]))
            peakTG2[1].append(tuneA400[j,np.mod(np.argmax(tuneA400[j])+4,8)])
            peakTG2[1].append(tuneB400[j,np.mod(np.argmax(tuneB400[j])+4,8)])
            print j+1,(tuneA400[j,np.mod(np.argmax(tuneA400[j])+4,8)]/np.max(tuneA400[j]))/(tuneA0[j,np.mod(np.argmax(tuneA0[j])+4,8)]/np.max(tuneA0[j]))
        elif idx_LPFC[j] == ncheck_LPFC[1]:
            peakCT1[0].append(np.max(tuneA0[j]))
            peakCT1[0].append(np.max(tuneB0[j]))
            peakCT2[0].append(tuneA0[j,np.mod(np.argmax(tuneA0[j])+4,8)])
            peakCT2[0].append(tuneB0[j,np.mod(np.argmax(tuneB0[j])+4,8)])
            peakCT1[1].append(np.max(tuneA400[j]))
            peakCT1[1].append(np.max(tuneB400[j]))
            peakCT2[1].append(tuneA400[j,np.mod(np.argmax(tuneA400[j])+4,8)])
            peakCT2[1].append(tuneB400[j,np.mod(np.argmax(tuneB400[j])+4,8)])
        elif idx_LPFC[j] == ncheck_LPFC[2] or idx_LPFC[j] == ncheck_LPFC[3]:
            peakTS1[0].append(np.max(tuneA0[j]))
            peakTS1[0].append(np.max(tuneB0[j]))
            peakTS2[0].append(tuneA0[j,np.mod(np.argmax(tuneA0[j])+4,8)])
            peakTS2[0].append(tuneB0[j,np.mod(np.argmax(tuneB0[j])+4,8)])
            peakTS1[1].append(np.max(tuneA400[j]))
            peakTS1[1].append(np.max(tuneB400[j]))
            peakTS2[1].append(tuneA400[j,np.mod(np.argmax(tuneA400[j])+4,8)])
            peakTS2[1].append(tuneB400[j,np.mod(np.argmax(tuneB400[j])+4,8)])
    for j in range(2):
        ax = pylab.subplot(5,4,17+j)
        pylab.plot([0,np.max([peakTS1[j],peakTS2[j]])+1],[0,np.max([peakTS1[j],peakTS2[j]])+1],'k--',lw=1)
        pylab.plot(peakTG1[j],peakTG2[j],'k.')
        pylab.plot(peakCT1[j],peakCT2[j],'k.')
        pylab.plot(peakTS1[j],peakTS2[j],'k.')
        if j == 0:
            pylab.ylabel('Opposite')
            pylab.title('Early')
        else:
            pylab.title('Late')
        pylab.xlabel('Main peak')
    # H model
    fig.text(0.5,0.25,'H',fontsize=fl_size,ha='right',va='bottom')
    fig.text(0.75,0.25,'Model',fontsize=fl_size,ha='right',va='bottom')
    peakTG1 = [[],[]]    # Early,Late: Main peak
    peakTG2 = [[],[]]    # Opp to peak
    peakCT1 = [[],[]]
    peakCT2 = [[],[]]
    peakTS1 = [[],[]]
    peakTS2 = [[],[]]
    tuneA0,tuneB0 = tuning_cj(datapath+fname,0)
    tuneA400,tuneB400 = tuning_cj(datapath+fname,400)
    tuneA0 *= 1000
    tuneB0 *= 1000
    tuneA400 *= 1000
    tuneB400 *= 1000
    if 1:
        print 'Baseline subtracted'
        for j in range(tuneA0.shape[0]):
            tuneA0[j] = tuneA0[j] - np.min(tuneA0[j])
            tuneB0[j] = tuneB0[j] - np.min(tuneB0[j])
            tuneA400[j] = tuneA400[j] - np.min(tuneA400[j])
            tuneB400[j] = tuneB400[j] - np.min(tuneB400[j])
    for j in range(len(idx)):
        if idx[j] == ncheck[0]: # 1 TG 2 CT 3 4 TS
            peakTG1[0].append(np.max(tuneA0[j]))
            peakTG1[0].append(np.max(tuneB0[j]))
            peakTG2[0].append(tuneA0[j,np.mod(np.argmax(tuneA0[j])+4,8)])
            peakTG2[0].append(tuneB0[j,np.mod(np.argmax(tuneB0[j])+4,8)])
            peakTG1[1].append(np.max(tuneA400[j]))
            peakTG1[1].append(np.max(tuneB400[j]))
            peakTG2[1].append(tuneA400[j,np.mod(np.argmax(tuneA400[j])+4,8)])
            peakTG2[1].append(tuneB400[j,np.mod(np.argmax(tuneB400[j])+4,8)])
        elif idx[j] == ncheck[1]:
            peakCT1[0].append(np.max(tuneA0[j]))
            peakCT1[0].append(np.max(tuneB0[j]))
            peakCT2[0].append(tuneA0[j,np.mod(np.argmax(tuneA0[j])+4,8)])
            peakCT2[0].append(tuneB0[j,np.mod(np.argmax(tuneB0[j])+4,8)])
            peakCT1[1].append(np.max(tuneA400[j]))
            peakCT1[1].append(np.max(tuneB400[j]))
            peakCT2[1].append(tuneA400[j,np.mod(np.argmax(tuneA400[j])+4,8)])
            peakCT2[1].append(tuneB400[j,np.mod(np.argmax(tuneB400[j])+4,8)])
        elif idx[j] == ncheck[2] or idx[j] == ncheck[3]:
            peakTS1[0].append(np.max(tuneA0[j]))
            peakTS1[0].append(np.max(tuneB0[j]))
            peakTS2[0].append(tuneA0[j,np.mod(np.argmax(tuneA0[j])+4,8)])
            peakTS2[0].append(tuneB0[j,np.mod(np.argmax(tuneB0[j])+4,8)])
            peakTS1[1].append(np.max(tuneA400[j]))
            peakTS1[1].append(np.max(tuneB400[j]))
            peakTS2[1].append(tuneA400[j,np.mod(np.argmax(tuneA400[j])+4,8)])
            peakTS2[1].append(tuneB400[j,np.mod(np.argmax(tuneB400[j])+4,8)])
    for j in range(2):
        ax = pylab.subplot(5,4,19+j)
        pylab.plot([0,np.max([peakTS1[j],peakTS2[j]])+1],[0,np.max([peakTS1[j],peakTS2[j]])+1],'k--',lw=1)
        pylab.plot(peakTG1[j],peakTG2[j],'k.')
        pylab.plot(peakCT1[j],peakCT2[j],'k.')
        pylab.plot(peakTS1[j],peakTS2[j],'k.')
        pylab.xlabel('Main peak')
        if j == 0:
            pylab.title('Early')
        else:
            pylab.title('Late')
    pylab.savefig(figpath+'s8.pdf')
    return
    """
    peakTG_LPFC = [[],[]]    # Early,Late: Opp/Main peak
    peakCT_LPFC = [[],[]]
    peakTS_LPFC = [[],[]]
    for j in range(len(idx_LPFC)):
        if idx_LPFC[j] == 1: # 1 TG 2 CT 3 4 TS
            peakTG_LPFC[0].append(tuneA0[j,np.mod(np.argmax(tuneA0[j])+4,8)]/np.max(tuneA0[j]))
            peakTG_LPFC[0].append(tuneB0[j,np.mod(np.argmax(tuneB0[j])+4,8)]/np.max(tuneB0[j]))
            peakTG_LPFC[1].append(tuneA400[j,np.mod(np.argmax(tuneA400[j])+4,8)]/np.max(tuneA400[j]))
            peakTG_LPFC[1].append(tuneB400[j,np.mod(np.argmax(tuneB400[j])+4,8)]/np.max(tuneB400[j]))
            print (tuneA400[j,np.mod(np.argmax(tuneA400[j])+4,8)]/np.max(tuneA400[j]))/(tuneA0[j,np.mod(np.argmax(tuneA0[j])+4,8)]/np.max(tuneA0[j]))
        elif idx_LPFC[j] == 2:
            peakCT_LPFC[0].append(tuneA0[j,np.mod(np.argmax(tuneA0[j])+4,8)]/np.max(tuneA0[j]))
            peakCT_LPFC[0].append(tuneB0[j,np.mod(np.argmax(tuneB0[j])+4,8)]/np.max(tuneB0[j]))
            peakCT_LPFC[1].append(tuneA400[j,np.mod(np.argmax(tuneA400[j])+4,8)]/np.max(tuneA400[j]))
            peakCT_LPFC[1].append(tuneB400[j,np.mod(np.argmax(tuneB400[j])+4,8)]/np.max(tuneB400[j]))
        elif idx_LPFC[j] >= 3:
            peakTS_LPFC[0].append(tuneA0[j,np.mod(np.argmax(tuneA0[j])+4,8)]/np.max(tuneA0[j]))
            peakTS_LPFC[0].append(tuneB0[j,np.mod(np.argmax(tuneB0[j])+4,8)]/np.max(tuneB0[j]))
            peakTS_LPFC[1].append(tuneA400[j,np.mod(np.argmax(tuneA400[j])+4,8)]/np.max(tuneA400[j]))
            peakTS_LPFC[1].append(tuneB400[j,np.mod(np.argmax(tuneB400[j])+4,8)]/np.max(tuneB400[j]))
    pylab.plot([0,1],[0,1],'k--')
    pylab.plot(peakTG_LPFC[0],peakTG_LPFC[1],'.',c=cTG)
    pylab.plot(peakCT_LPFC[0],peakCT_LPFC[1],'.',c=cCT)
    pylab.plot(peakTS_LPFC[0],peakTS_LPFC[1],'.',c=cTS1)
    pylab.xlabel('Early ratio')
    pylab.ylabel('Late ratio')
    if 1:
        print 'Baseline subtracted'
        for j in range(tuningE0.shape[0]):
            tuningE0[j] = tuningE0[j] - np.min(tuningE0[j])
            tuningE1[j] = tuningE1[j] - np.min(tuningE1[j])
            tuningL0[j] = tuningL0[j] - np.min(tuningL0[j])
            tuningL1[j] = tuningL1[j] - np.min(tuningL1[j])
    ax = pylab.subplot(4,3,11)
    peakTG = [[],[]]    # Early,Late: Opp/Main peak
    peakCT = [[],[]]
    peakTS = [[],[]]
    for j in range(len(idx)):
        if idx[j] == ncheck[0]: # 1 TG 2 CT 3 4 TS
            peakTG[0].append(tuningE0[j,np.mod(np.argmax(tuningE0[j])+4,8)]/np.max(tuningE0[j]))
            peakTG[0].append(tuningE1[j,np.mod(np.argmax(tuningE1[j])+4,8)]/np.max(tuningE1[j]))
            peakTG[1].append(tuningL0[j,np.mod(np.argmax(tuningL0[j])+4,8)]/np.max(tuningL0[j]))
            peakTG[1].append(tuningL1[j,np.mod(np.argmax(tuningL1[j])+4,8)]/np.max(tuningL1[j]))
        elif idx[j] == ncheck[1]:
            peakCT[0].append(tuningE0[j,np.mod(np.argmax(tuningE0[j])+4,8)]/np.max(tuningE0[j]))
            peakCT[0].append(tuningE1[j,np.mod(np.argmax(tuningE1[j])+4,8)]/np.max(tuningE1[j]))
            peakCT[1].append(tuningL0[j,np.mod(np.argmax(tuningL0[j])+4,8)]/np.max(tuningL0[j]))
            peakCT[1].append(tuningL1[j,np.mod(np.argmax(tuningL1[j])+4,8)]/np.max(tuningL1[j]))
        elif idx[j] == ncheck[2] or idx[j] == ncheck[3]:
            peakTS[0].append(tuningE0[j,np.mod(np.argmax(tuningE0[j])+4,8)]/np.max(tuningE0[j]))
            peakTS[0].append(tuningE1[j,np.mod(np.argmax(tuningE1[j])+4,8)]/np.max(tuningE1[j]))
            peakTS[1].append(tuningL0[j,np.mod(np.argmax(tuningL0[j])+4,8)]/np.max(tuningL0[j]))
            peakTS[1].append(tuningL1[j,np.mod(np.argmax(tuningL1[j])+4,8)]/np.max(tuningL1[j]))
    pylab.plot([0,1],[0,1],'k--')
    pylab.plot(peakTG[0],peakTG[1],'.',c=cTG)
    pylab.plot(peakCT[0],peakCT[1],'.',c=cCT)
    pylab.plot(peakTS[0],peakTS[1],'.',c=cTS1)"""

def addition1():
    peak0 = pylab.empty(11)
    peak1 = pylab.empty(11)
    tune0 = []
    tune1 = []
    idx = pylab.loadtxt(datapath+'idx_selected.txt')
    selectedID = range(idx.size) #[56,541,53]
    for j in range(3):
        tune = pylab.loadtxt(matlabpath+'tuningA'+str(j*200)+'.txt')
        tune0.append(tune[selectedID,:])
        tune = pylab.loadtxt(matlabpath+'tuningB'+str(j*200)+'.txt')
        tune1.append(tune[selectedID,:])
    for nf in range(idx.size/12+1):
        fig = pylab.figure(figsize=[18,10])
        if nf == idx.size/12:
            krange = 5
        else:
            krange = 12
        for k in range(krange):
            id = selectedID[nf*12+k]
            print id
            for j in range(11):
                tune = pylab.loadtxt(matlabpath+'tuningA'+str(j*50)+'.txt')
                peak0[j] = prefdir(tune[id,:])
                tune = pylab.loadtxt(matlabpath+'tuningB'+str(j*50)+'.txt')
                peak1[j] = prefdir(tune[id,:])
            ax = pylab.subplot(6,8,k*4+4)
            pylab.plot(pylab.arange(100,601,50),peak0,c=cA,lw=2,label='A chosen')
            pylab.plot(pylab.arange(100,601,50),peak1,c=cB,lw=2,label='B chosen')
            if k == 0:
                pylab.legend(loc=1,frameon=False)
            pylab.xlim(0,700)
            pylab.ylim(0,363)
            pylab.xticks([])
            #pylab.yticks(range(0,360,90),['0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
            pylab.yticks(range(0,361,180),['0$^{\circ}$','180$^{\circ}$','360$^{\circ}$'])
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')
            if k <= 1:
                pylab.title('Peak')
            if k >= 10:
                pylab.xticks(range(0,601,300))
                pylab.xlabel('Time (ms)')
            #pylab.xlabel('Time from target onset (ms)')
            for j in range(3):
                ax = pylab.subplot(6,8,k*4+1+j)
                fc1 = sp.interp1d(range(-180,360+180,45),np.append(tune0[j][id][-4:],np.append(tune0[j][id],tune0[j][id][:4])),kind='cubic')
                fc2 = sp.interp1d(range(-180,360+180,45),np.append(tune1[j][id][-4:],np.append(tune1[j][id],tune1[j][id][:4])),kind='cubic')
                pylab.plot(range(0,360,45),tune0[j][id],'o',c=cA)
                pylab.plot(range(0,361,15),fc1(range(0,361,15))[:],c=cA,lw=2)
                pylab.plot(range(0,360,45),tune1[j][id],'o',c=cB)
                pylab.plot(range(0,361,15),fc2(range(0,361,15))[:],c=cB,lw=2)
                pylab.plot(prefdir(tune0[j][id])*pylab.ones(2),[-20,120],c=cA)
                pylab.plot(prefdir(tune1[j][id])*pylab.ones(2),[-20,120],c=cB)
                pylab.xticks([])
                pylab.xlim(-10,370)
                if max(np.append(tune0[j][id],tune1[j][id]))-min(np.append(tune0[j][id],tune1[j][id]))>70:
                    pylab.yticks(range(0,100,50))
                    pylab.ylim(min(np.append(tune0[j][id],tune1[j][id]))-20,max(np.append(tune0[j][id],tune1[j][id]))+20)
                elif max(np.append(tune0[j][id],tune1[j][id]))-min(np.append(tune0[j][id],tune1[j][id]))>35:
                    pylab.yticks(range(0,100,20))
                    pylab.ylim(min(np.append(tune0[j][id],tune1[j][id]))-10,max(np.append(tune0[j][id],tune1[j][id]))+10)
                else:
                    pylab.yticks(range(0,100,10))
                    pylab.ylim(min(np.append(tune0[j][id],tune1[j][id]))-5,max(np.append(tune0[j][id],tune1[j][id]))+5)
                if k >= 10:
                    pylab.xlabel('Target A')
                    pylab.xticks(range(0,360,90),['0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
                    if j == 0:
                        pylab.ylabel('Rate (Hz)')
                ax.spines['right'].set_visible(False)
                ax.spines['top'].set_visible(False)
                ax.xaxis.set_ticks_position('bottom')
                ax.yaxis.set_ticks_position('left')
                if k <= 1:
                    if j == 0:
                        pylab.title('Early')
                    elif j == 1:
                        pylab.title('Mid')
                    elif j == 2:
                        pylab.title('Late')
        pylab.subplots_adjust(top=0.9,bottom=0.1,left=0.05,right=0.95,wspace=0.5)
        pylab.savefig('LPFC_tuning'+str(nf+1)+'.pdf')
        pylab.close("all")

def addition2():
    peak0 = pylab.empty(11)
    peak1 = pylab.empty(11)
    selectedID = range(n*3) #[56,541,53]
    tune0 = []
    tune1 = []
    J_p,J_m,g_wt,g_tt,wmIno0,ctIno0,std_J,std_t,std_I,std_sig = para_heterogeneous()
    #fname = 'jp'+str(J_p)+'jm'+str(-J_m)+'gwt'+str(g_wt)+'gtt'+str(g_tt)+'I'+str(wmIno0)+'std'+str(std_J)+'seednet'+str(seednet)+'seedrun'+str(seedrun)
    fname = 'input2pk_jp2.5jm1.6gwt0.05gst0.243gtt0.3b0.9I0.31std2.0seednet101seedrun51'
    for j in range(3):
        tuning0, tuning1 = tuning_cj(datapath+fname,200*j)
        tune0.append(tuning0[selectedID,:]*1000)
        tune1.append(tuning1[selectedID,:]*1000)
    for nf in range(n*3/12):
        fig = pylab.figure(figsize=[18,10])
        for k in range(12):
            id = selectedID[nf*12+k]
            print id
            for j in range(11):
                tuning0, tuning1 = tuning_cj(datapath+fname,50*j)
                peak0[j] = prefdir(tuning0[selectedID[nf*12+k],:])
                peak1[j] = prefdir(tuning1[selectedID[nf*12+k],:])
            ax = pylab.subplot(6,8,k*4+4)
            pylab.plot(pylab.arange(100,601,50),peak0,c=cA,lw=2,label='A chosen')
            pylab.plot(pylab.arange(100,601,50),peak1,c=cB,lw=2,label='B chosen')
            if k == 0:
                pylab.legend(loc=1,frameon=False)
            pylab.xlim(0,700)
            pylab.ylim(0,363)
            pylab.xticks([])
            #pylab.yticks(range(0,360,90),['0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
            pylab.yticks(range(0,361,180),['0$^{\circ}$','180$^{\circ}$','360$^{\circ}$'])
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')
            if k <= 1:
                pylab.title('Peak')
            if k >= 10:
                pylab.xticks(range(0,601,300))
                pylab.xlabel('Time (ms)')
                #pylab.xlabel('Time from target onset (ms)')
            for j in range(3):
                ax = pylab.subplot(6,8,k*4+1+j)
                fc1 = sp.interp1d(range(-180,360+180,45),np.append(tune0[j][id][-4:],np.append(tune0[j][id],tune0[j][id][:4])),kind='cubic')
                fc2 = sp.interp1d(range(-180,360+180,45),np.append(tune1[j][id][-4:],np.append(tune1[j][id],tune1[j][id][:4])),kind='cubic')
                pylab.plot(range(0,360,45),tune0[j][id],'o',c=cA)
                pylab.plot(range(0,361,15),fc1(range(0,361,15))[:],c=cA,lw=2)
                pylab.plot(range(0,360,45),tune1[j][id],'o',c=cB)
                pylab.plot(range(0,361,15),fc2(range(0,361,15))[:],c=cB,lw=2)
                pylab.plot(prefdir(tune0[j][id])*pylab.ones(2),[-20,120],c=cA)
                pylab.plot(prefdir(tune1[j][id])*pylab.ones(2),[-20,120],c=cB)
                pylab.xticks([])
                pylab.xlim(-10,370)
                if max(np.append(tune0[j][id],tune1[j][id]))-min(np.append(tune0[j][id],tune1[j][id]))>70:
                    pylab.yticks(range(0,100,50))
                    pylab.ylim(min(np.append(tune0[j][id],tune1[j][id]))-20,max(np.append(tune0[j][id],tune1[j][id]))+20)
                elif max(np.append(tune0[j][id],tune1[j][id]))-min(np.append(tune0[j][id],tune1[j][id]))>35:
                    pylab.yticks(range(0,100,20))
                    pylab.ylim(min(np.append(tune0[j][id],tune1[j][id]))-10,max(np.append(tune0[j][id],tune1[j][id]))+10)
                else:
                    pylab.yticks(range(0,100,10))
                    pylab.ylim(min(np.append(tune0[j][id],tune1[j][id]))-5,max(np.append(tune0[j][id],tune1[j][id]))+5)
                if k >= 10:
                    pylab.xlabel('Target A')
                    pylab.xticks(range(0,360,90),['0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
                    if j == 0:
                        pylab.ylabel('Rate (Hz)')
                ax.spines['right'].set_visible(False)
                ax.spines['top'].set_visible(False)
                ax.xaxis.set_ticks_position('bottom')
                ax.yaxis.set_ticks_position('left')
                if k <= 1:
                    if j == 0:
                        pylab.title('Early')
                    elif j == 1:
                        pylab.title('Mid')
                    elif j == 2:
                        pylab.title('Late')
        pylab.subplots_adjust(top=0.9,bottom=0.1,left=0.05,right=0.95,wspace=0.5)
        pylab.savefig('tuning_example'+str(nf+1)+'.pdf')
        pylab.close("all")

supp8()
#addition2()
pylab.show()
