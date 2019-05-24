"""
This code will run all the figures of the paper from scratch.
"""

import numpy as np
import pylab
import matplotlib as mpl
import scipy.interpolate as sp
import os.path

execfile('lpfc_para.py')
if not os.path.exists(datapath):
    os.makedirs(datapath)
if not os.path.exists(figpath):
    os.makedirs(figpath)
execfile('lpfc_model.py')
execfile('lpfc_data.py')

### Figure 1
def figure1():
    print 'Making figure 1...'
    fig = pylab.figure(figsize=[10,6.5])
    #vgauss = visgauss(90,270,sigma_s)
    # A
    fig.text(0.05,0.95,'A',fontsize=fl_size,ha='right',va='bottom')
    # B
    fig.text(0.05,0.65,'B',fontsize=fl_size,ha='right',va='bottom')
    # C
    fig.text(0.5,0.65,'C',fontsize=fl_size,ha='right',va='bottom')
    ax = pylab.subplot(347)
    mpl.rcParams['axes.titlesize'] = font_size
    mpl.rcParams['xtick.labelsize'] = font_size-2
    mpl.rcParams['ytick.labelsize'] = font_size-2
    mpl.rcParams['axes.labelsize'] = font_size-2
    pylab.bar(range(2),jI,color=[cA,cB])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    pylab.ylabel('Input')
    pylab.title('Chosen juice')
    #pylab.ylabel('Chosen\njuice')
    pylab.yticks([])
    pylab.xlim([-0.9,1.8])
    pylab.xticks([0.,1.],('A','B'))
    #pylab.xlabel('WM')
    # D
    fig.text(0.72,0.65,'D',fontsize=fl_size,ha='right',va='bottom')
    ax = pylab.subplot(348)
    profile = np.zeros(alldir.size)
    for j in range(alldir.size):
        profile[j] = gauss(tgt,alldir[j])
    pylab.plot(alldir,profile,cA,linewidth=2)
    for j in range(alldir.size):
        profile[j] = gauss(opp,alldir[j])
    pylab.plot(alldir,profile,cB,linewidth=2)
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
    # E
    fig.text(0.5,0.28,'E',fontsize=fl_size,ha='right',va='bottom')
    ax = pylab.subplot(428)
    item = pylab.loadtxt(datapath+'lpfc_wm_CJ_0_.txt')
    wmdur = item.size/2/tp
    pylab.plot(np.arange(0,wmdur,tp),item[:,0]*1000,cA)
    pylab.plot(np.arange(0,wmdur,tp),item[:,1]*1000,cB)
    pylab.ylabel('WM')
    pylab.xlabel('Time')
    pylab.yticks([0,30,60])
    pylab.xticks(np.cumsum([preoffer+1000,offeron,offeroff,targeton,go]),['','','','',''])
    pylab.xlim(0,dur+1000)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    pylab.savefig(figpath+'fig1.pdf')

### Figure 2
def figure2():
    print 'Making figure 2...'
    #mpl.rcParams['font.size'] = font_size
    fig = pylab.figure(figsize=[12,8])
    pylab.subplots_adjust(top=0.9,bottom=0.1,hspace=0.2,wspace=0.2)
    # A
    fig.text(0.05,0.95,'A',fontsize=fl_size,ha='right',va='bottom')
    # Interaction profile - B
    gau = np.zeros(n+1)
    for k in range(n+1):
        gau[k] = gauss(degdiff(-n/2,k))
    pylab.subplot(541)
    J_p,J_m,g_wt,g_tt,wmIno0,ctIno0,std_J,std_t,std_I,std_sig = para_homogeneous_scene1()
    pylab.plot(range(-n/2,n/2+1),J_m/2.+(1-g_tt/2.)*J_p*gau,'m',lw=2,label='within')
    pylab.plot(range(-n/2,n/2+1),J_m/2.+g_tt*J_p*gau/2.,'c',lw=2,label='between')
    pylab.plot([0,0],[-0.5,2],'k')
    pylab.plot([-n/2,n/2+1],[0,0],'k')
    pylab.axis('off')
    pylab.xticks([-n/2,0,n/2],['-180$^{\circ}$','0$^{\circ}$','180$^{\circ}$'])
    pylab.yticks([0])
    pylab.xlim(-n/2,n/2+1)
    pylab.title('Scenario I')
    pylab.xlabel(r'$\Delta \theta (^{\circ})$')
    pylab.ylabel('Strength')
    # Interaction profile - C
    ax = pylab.subplot(545)
    J_p,J_m,g_wt,g_tt,wmIno0,ctIno0,std_J,std_t,std_I,std_sig = para_homogeneous_scene2()
    pylab.plot(range(-n/2,n/2+1),J_m/2.+(1-g_tt/2.)*J_p*gau,'m',lw=2,label='within')
    pylab.plot(range(-n/2,n/2+1),J_m/2.+g_tt*J_p*gau/2.,'c',lw=2,label='between')
    pylab.legend(loc=(-0.15,0.5),frameon=False)
    pylab.plot(range(-n/2,n/2+1),J_m/2.+(1-g_tt/2.)*J_p*gau,'m--',dashes=(5,5),lw=2,label='within')
    pylab.plot([0,0],[-0.5,2],'k')
    pylab.plot([-n/2,n/2+1],[0,0],'k')
    pylab.xticks([-n/2,0,n/2],['-180$^{\circ}$','0$^{\circ}$','180$^{\circ}$'])
    pylab.yticks([0])
    pylab.xlim(-n/2,n/2+1)
    pylab.xlabel(r'$\Delta\theta$')
    pylab.title('Scenario II')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    # B
    fig.text(0.33,0.95,'B',fontsize=fl_size,ha='right',va='bottom')
    pylab.subplot(632)
    J_p,J_m,g_wt,g_tt,wmIno0,ctIno0,std_J,std_t,std_I,std_sig = para_homogeneous_scene1()
    fname = 'jp'+str(J_p)+'jm'+str(-J_m)+'gwt'+str(g_wt)+'gtt'+str(g_tt)+'I'+str(wmIno0)+'std'+str(std_J)+'seednet'+str(seednet)+'seedrun'+str(seedrun)
    item = pylab.loadtxt(datapath+fname+'/lpfc_t1_net_101_run_51_CJ_0_d_90_.txt')
    item = item[0:int(targeton/tp)+1,:]
    pylab.pcolor(item.T*1000,rasterized=True,cmap='jet')
    pylab.clim(0,60)
    cb = pylab.colorbar(ticks=[])
    cb.remove()
    pylab.barh([180],[200],20,linewidth=0,color=cEarly)
    pylab.barh([180],[200],20,[200],linewidth=0,color=cMid)
    pylab.barh([180],[200],20,[400],linewidth=0,color=cLate)
    pylab.ylabel('IN-A')
    pylab.text(-40,n+50,'Early',color=cEarly)
    pylab.text(215,n+50,'Mid',color=cMid)
    pylab.text(410,n+50,'Late',color=cLate)
    pylab.yticks(range(0,n+1,n/2),['0$^{\circ}$','180$^{\circ}$','360$^{\circ}$'])
    pylab.ylim(0,n)
    pylab.xlim(0,targeton)
    pylab.xticks([])
    pylab.subplot(635)
    item = pylab.loadtxt(datapath+fname+'/lpfc_t2_net_'+str(seednet)+'_run_'+str(seedrun)+'_CJ_0_d_90_.txt')
    item = item[0:int(targeton/tp)+1,:]
    pylab.pcolor(item.T*1000,rasterized=True,cmap='jet')
    pylab.clim(0,60)
    cb = pylab.colorbar(ticks=[])
    cb.remove()
    pylab.ylabel('IN-B')
    pylab.yticks(range(0,n+1,n/2),['0$^{\circ}$','180$^{\circ}$','360$^{\circ}$'])
    pylab.ylim(0,n)
    pylab.xlim(0,targeton)
    pylab.xticks([])
    pylab.subplot(638)
    item = pylab.loadtxt(datapath+fname+'/lpfc_ct_net_'+str(seednet)+'_run_'+str(seedrun)+'_CJ_0_d_90_.txt')
    item = item[0:int(targeton/tp)+1,:]
    pylab.pcolor(item.T*1000,rasterized=True,cmap='jet')
    pylab.clim(0,60)
    cb = pylab.colorbar(ticks=[])
    cb.remove()
    pylab.ylabel('RO')
    pylab.yticks(range(0,n+1,n/2),['0$^{\circ}$','180$^{\circ}$','360$^{\circ}$'])
    pylab.ylim(0,n)
    pylab.xlim(0,targeton)
    pylab.xticks(np.arange(0,targeton+1,500),['0','500','1000'])
    pylab.xlabel('Time from target onset (ms)')
    # C
    fig.text(0.6,0.95,'C',fontsize=fl_size,ha='right',va='bottom')
    pylab.subplot(633)
    J_p,J_m,g_wt,g_tt,wmIno0,ctIno0,std_J,std_t,std_I,std_sig = para_homogeneous_scene2()
    fname = 'jp'+str(J_p)+'jm'+str(-J_m)+'gwt'+str(g_wt)+'gtt'+str(g_tt)+'I'+str(wmIno0)+'std'+str(std_J)+'seednet'+str(seednet)+'seedrun'+str(seedrun)
    item = pylab.loadtxt(datapath+fname+'/lpfc_t1_net_'+str(seednet)+'_run_'+str(seedrun)+'_CJ_0_d_90_.txt')
    item = item[:int(targeton/tp)+1,:]
    pylab.pcolor(item.T*1000,rasterized=True,cmap='jet')
    pylab.clim(0,60)
    cbar = pylab.colorbar(ticks=[0,30,60])
    cbar.set_label('Rate (Hz)')
    pylab.barh([180],[200],20,linewidth=0,color=cEarly)
    pylab.barh([180],[200],20,[200],linewidth=0,color=cMid)
    pylab.barh([180],[200],20,[400],linewidth=0,color=cLate)
    pylab.ylabel('IN-A')
    pylab.text(-40,n+50,'Early',color=cEarly)
    pylab.text(215,n+50,'Mid',color=cMid)
    pylab.text(410,n+50,'Late',color=cLate)
    pylab.yticks([])
    pylab.xticks([])
    pylab.xlim(0,targeton)
    pylab.ylim(0,n)
    pylab.subplot(636)
    item = pylab.loadtxt(datapath+fname+'/lpfc_t2_net_'+str(seednet)+'_run_'+str(seedrun)+'_CJ_0_d_90_.txt')
    item = item[0:int(targeton/tp)+1,:]
    pylab.pcolor(item.T*1000,rasterized=True,cmap='jet')
    pylab.clim(0,60)
    cb = pylab.colorbar(ticks=[])
    cb.remove()
    pylab.ylabel('IN-B')
    pylab.yticks([])
    pylab.xticks([])
    pylab.xlim(0,targeton)
    pylab.ylim(0,n)
    pylab.subplot(639)
    item = pylab.loadtxt(datapath+fname+'/lpfc_ct_net_'+str(seednet)+'_run_'+str(seedrun)+'_CJ_0_d_90_.txt')
    item = item[0:int(targeton/tp)+1,:]
    pylab.pcolor(item.T*1000,rasterized=True,cmap='jet')
    pylab.clim(0,60)
    cb = pylab.colorbar(ticks=[])
    cb.remove()
    pylab.ylabel('RO')
    pylab.yticks([])
    pylab.xticks([])
    pylab.xlim(0,targeton)
    pylab.ylim(0,n)
    pylab.xticks(np.arange(0,targeton+1,500),['0','500','1000'])
    pylab.xlabel('Time from target onset (ms)')
    # D
    fig.text(0.15,0.35,'D',fontsize=fl_size,ha='right',va='bottom')
    ax = fig.add_axes([0.2, 0.1, 0.25, 0.25])
    item = pylab.loadtxt(datapath+fname+'/lpfc_t2_net_'+str(seednet)+'_run_'+str(seedrun)+'_CJ_0_d_90_.txt')
    pylab.plot(alldir[:-1],1000*np.mean(item[:200/tp,:],0),color=cEarly,linewidth=2,label='Early')
    pylab.plot(alldir[:-1],1000*np.mean(item[200/tp:400/tp,:],0),color=cMid,linewidth=2,label='Mid')
    pylab.plot(alldir[:-1],1000*np.mean(item[400/tp:600/tp,:],0),color=cLate,linewidth=2,label='Late')
    pylab.legend(loc=1,frameon=False)
    pylab.xlabel('Neuron')
    pylab.ylabel('Rate (Hz)')
    pylab.xticks(range(0,361,180),['0$^{\circ}$','180$^{\circ}$','360$^{\circ}$'])
    pylab.xlim(0,361)
    pylab.ylim(0,40)
    pylab.yticks([0,15,30])
    pylab.title('IN-B of Scenario II')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    # E
    fig.text(0.5,0.35,'E',fontsize=fl_size,ha='right',va='bottom')
    ax = fig.add_axes([0.55, 0.1, 0.25, 0.25])
    grange = np.arange(0.,1.1,0.1)
    if overwrite or not os.path.exists(datapath+'fig2E.txt'):
        dangle = np.zeros((3,11))
        J_p,J_m = para_homogeneous_scan()
        for j in range(3):
            if j == 0:
                kstart = 0
            elif j == 1:
                kstart = 0
            elif j == 2:
                kstart = 0
            for k in range(kstart,11):
                g_tt = 0.1*k
                fname = 'jp'+str(J_p[j])+'jm'+str(-J_m[j])+'gwt'+str(g_wt)+'gtt'+str(g_tt)+'I'+str(wmIno0)+'std'+str(std_J)+'seednet'+str(seednet)+'seedrun'+str(seedrun)
                dangle[j,k] = time_transition(datapath+fname)
                print j,k,dangle[j,k]
                #ft2 = pylab.loadtxt(datapath+fname+'/lpfc_t2_net_'+str(seednet)+'_run_'+str(seedrun)+'_CJ_0_d_90_.txt')
                #t2 = ft2[400/tp:600/tp,:]
                ##temp,decoded_t2 = resultant(np.mean(t2,0),div=n)
                #fct = pylab.loadtxt(datapath+fname+'/lpfc_ct_net_'+str(seednet)+'_run_'+str(seedrun)+'_CJ_0_d_90_.txt')
                #ct = fct[400/tp:600/tp,:]
                ##temp,decoded_ct = resultant(np.mean(ct,0),div=n)
                ##dangle[j,k] = degdiff(decoded_t2,decoded_ct)
                #dangle[j,k] = anglediff(prefdir(np.mean(t2,0),n),prefdir(np.mean(ct,0),n))
        ftran = open(datapath+'fig2E.txt','w')
        for j in range(3):
            for k in range(11):
                ftran.write(str(dangle[j,k])+'\t')
            ftran.write('\n')
        ftran.close()
    else:
        dangle = np.loadtxt(datapath+'fig2E.txt')
    pylab.plot(grange,dangle[0,:],'^',ms=10,mew=0,label='strong inh')
    pylab.plot(grange,dangle[1,:],'+',ms=10,mew=1.5,label='reference')
    pylab.plot(grange,dangle[2,:],'x',ms=8,mew=1.5,label='strong exc')
    pylab.legend(loc=6,numpoints=1,frameon=False)
    #pylab.text(0.3,30,'Transition')
    pylab.xlabel(r'$\alpha$')
    #pylab.ylabel('$\Delta p$')
    pylab.ylabel('$1/T_{tran}$ (s$^{-1}$)')
    pylab.xticks([0,0.5,1.],['0','0.5','1'])
    pylab.xlim(-0.05,1.05)
    pylab.ylim(-0.2,3.6)
    pylab.yticks(range(0,4,1))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    pylab.savefig(figpath+'fig2.pdf')

### Figure 3
def figure3():
    print 'Making figure 3...'
    fig = pylab.figure(figsize=[12,15])
    font_size = 16
    mpl.rcParams['axes.titlesize'] = font_size+2
    mpl.rcParams['xtick.labelsize'] = font_size
    mpl.rcParams['ytick.labelsize'] = font_size
    mpl.rcParams['axes.labelsize'] = font_size
    mpl.rcParams['legend.fontsize'] = font_size-2
    mpl.rcParams['font.size'] = font_size
    #eps = 65
    #min_samples = 20
    #epilson = 30
    dotsize = 60
    # A
    fig.text(0.05,0.95,'A',fontsize=fl_size,ha='right',va='bottom')
    # B
    fig.text(0.05,0.84,'B',fontsize=fl_size,ha='right',va='bottom')
    profile1 = np.zeros(alldir.size)
    profile2 = np.zeros(alldir.size)
    for j in range(alldir.size):
        profile1[j] = 0.8*gauss(tgt,alldir[j])+0.05
        profile2[j] = 0.8*gauss(opp,alldir[j])+0.05
    ax = pylab.subplot(834)
    pylab.plot(alldir,profile1,cA,linewidth=2,label='A chosen')
    pylab.plot([-1],[-1],cB,linewidth=2,label='B chosen')
    pylab.plot(alldir,profile1,cB,linestyle='--',linewidth=2,dashes=(5,5))
    pylab.legend(loc=(0.45,0.35),frameon=False)
    #pylab.plot(alldir,profile1,cB,linestyle='--',linewidth=2)
    pylab.xlim(0,360)
    pylab.ylim(0,1)
    pylab.xticks([])
    pylab.yticks([])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    pylab.ylabel('Rate')
    pylab.title('Early')
    ax = pylab.subplot(835)
    profile = np.zeros(alldir.size)
    for j in range(alldir.size):
        profile[j] = 0.8*gauss(tgt,alldir[j])+0.05
    pylab.plot(alldir,profile,cA,linewidth=2)
    pylab.plot(alldir,profile,cB,linestyle='--',linewidth=2,dashes=(5,5))
    pylab.xlim(0,360)
    pylab.ylim(0,1)
    pylab.xticks([])
    pylab.yticks([])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    pylab.title('Late')
    ax = pylab.subplot(836)
    pylab.plot([0,1],[0.2,0.2],cA,linewidth=2)
    pylab.plot([0,1],[0.2,0.2],cB,linestyle='--',linewidth=2,dashes=(5,5))
    pylab.title('Evolution')
    pylab.xlim(0,1)
    pylab.xticks([])
    pylab.ylim(0,1)
    pylab.yticks([])
    pylab.ylabel('Peak')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.text(0.09,0.785,'TG',fontsize=mpl.rcParams['axes.titlesize'],ha='right',va='bottom',rotation='vertical')
    fig.text(0.09,0.68,'CT',fontsize=mpl.rcParams['axes.titlesize'],ha='right',va='bottom',rotation='vertical')
    fig.text(0.09,0.57,'TS',fontsize=mpl.rcParams['axes.titlesize'],ha='right',va='bottom',rotation='vertical')
    # C
    fig.text(0.05,0.71,'C',fontsize=fl_size,ha='right',va='bottom')
    ax = pylab.subplot(837)
    pylab.plot(alldir,profile1,cA,linewidth=2)
    pylab.plot(alldir,profile2,cB,linewidth=2)
    pylab.xlim(0,360)
    pylab.ylim(0,1)
    pylab.xticks([])
    pylab.yticks([])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    pylab.ylabel('Rate')
    pylab.yticks([])
    ax = pylab.subplot(838)
    for j in range(alldir.size):
        profile[j] = 0.8*gauss(tgt,alldir[j])+0.05
    pylab.plot(alldir,profile,cA,linewidth=2)
    for j in range(alldir.size):
        profile[j] = 0.8*gauss(opp,alldir[j])+0.05
    pylab.plot(alldir,profile,cB,linewidth=2)
    pylab.plot([90,270],[0.85,0.85],'k')
    pylab.plot([90],[0.85],'k<')
    pylab.plot([265],[0.85],'k>')
    pylab.text(158,0.88,'180$^{\circ}$')
    pylab.xlim(0,360)
    #pylab.xticks([0,180,360])
    pylab.ylim(0,1)
    pylab.xticks([])
    pylab.yticks([])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    pylab.yticks([])
    ax = pylab.subplot(839)
    pylab.plot([0,1],[0.2,0.2],cA,linewidth=2)
    pylab.plot([0,1],[0.8,0.8],cB,linewidth=2)
    pylab.plot([0.5,0.5],[0.2,0.8],'k')
    pylab.plot([0.5],[0.21],'kv')
    pylab.plot([0.5],[0.77],'k^')
    pylab.text(0.51,0.47,'180$^{\circ}$')
    pylab.xlim(0,1)
    pylab.xticks([])
    pylab.ylim(0,1)
    pylab.yticks([])
    pylab.ylabel('Peak')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    # D
    fig.text(0.05,0.61,'D',fontsize=fl_size,ha='right',va='bottom')
    ax = pylab.subplot(8,3,10)
    profile = np.zeros(alldir.size)
    for j in range(alldir.size):
        profile[j] = 0.8*gauss(tgt,alldir[j])+0.05
    pylab.plot(alldir,profile,cA,linewidth=2)
    pylab.plot(alldir,profile,cB,linestyle='--',linewidth=2)
    pylab.xlim(0,360)
    pylab.ylim(0,1)
    pylab.xlabel(r'Target A')
    pylab.xticks([])
    pylab.yticks([])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    pylab.ylabel('Rate')
    ax = pylab.subplot(8,3,11)
    for j in range(alldir.size):
        profile[j] = 0.8*gauss(tgt,alldir[j])+0.05
    pylab.plot(alldir,profile,cA,linewidth=2)
    for j in range(alldir.size):
        profile[j] = 0.8*gauss(opp,alldir[j])+0.05
    pylab.plot(alldir,profile,cB,linewidth=2)
    pylab.xlim(0,360)
    pylab.ylim(0,1)
    pylab.xticks([])
    pylab.yticks([])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    pylab.xlabel(r'Target A')
    pylab.yticks([])
    ax = pylab.subplot(8,3,12)
    pylab.plot([0,1],[0.2,0.2],cA,linewidth=2)
    pylab.plot([0,0.3],[0.2,0.2],cB,linestyle='--',linewidth=2)
    pylab.plot([0.31,0.31,1],[0.2,0.8,0.8],cB,linewidth=2)
    pylab.xlim(0,1)
    pylab.xticks([])
    pylab.ylim(0,1)
    pylab.yticks([])
    pylab.ylabel('Peak')
    pylab.xlabel('Time from target onset')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    # E
    fig.text(0.05,0.5,'E',fontsize=fl_size,ha='right',va='bottom')
    pylab.subplot(437)
    J_p,J_m,g_wt,g_tt,wmIno0,ctIno0,std_J,std_t,std_I,std_sig = para_homogeneous_scene1()
    fname = 'jp'+str(J_p)+'jm'+str(-J_m)+'gwt'+str(g_wt)+'gtt'+str(g_tt)+'I'+str(wmIno0)+'std'+str(std_J)+'seednet'+str(seednet)+'seedrun'+str(seedrun)
    tuningE0, tuningE1, tuningL0, tuningL1 = tuning_curves(path=datapath+fname)
    peakdata = peak_tuning(tuningE0, tuningE1, tuningL0, tuningL1)
    idx,ncount,ncheck,uc,nuc = clustering(peakdata,eps=eps_sim)
    uc.append(0)
    pylab.plot([-90,270],[-90,270],'k--')
    for j in range(len(uc)):
        pylab.scatter(peakdata[idx==uc[j],0],peakdata[idx==uc[j],1],color='0.5',s=dotsize,lw=0)
    if ncheck[0] != 0:
        pylab.scatter(peakdata[idx==ncheck[0],0],peakdata[idx==ncheck[0],1],c=cTG,s=dotsize,lw=0)
    if ncheck[1] != 0:
        pylab.scatter(peakdata[idx==ncheck[1],0],peakdata[idx==ncheck[1],1],c=cCT,s=dotsize,lw=0)
    if ncheck[2] != 0:
        pylab.scatter(peakdata[idx==ncheck[2],0],peakdata[idx==ncheck[2],1],c=cTS1,s=dotsize,lw=0)
    if ncheck[3] != 0:
        pylab.scatter(peakdata[idx==ncheck[3],0],peakdata[idx==ncheck[3],1],c=cTS2,s=dotsize,lw=0)
    #pylab.scatter([0,180],[0,180],marker='+',c='k',s=100,label='Theory')
    #pylab.scatter([0,180],[0,180],marker='+',c='k',s=100,lw=2,label='Theory')
    #pylab.legend(loc=4,scatterpoints=1,frameon=False)
    pylab.xlim(-95,275)
    pylab.ylim(-95,275)
    pylab.xticks([])
    pylab.yticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.xlabel('$\Delta_{Early}$')
    pylab.ylabel('$\Delta_{Late}$')
    pylab.subplot(438)
    pylab.plot([-90,270],[-90,270],'k--')
    for j in range(len(uc)):
        pylab.scatter(peakdata[idx==uc[j],2],peakdata[idx==uc[j],3],color='0.5',s=dotsize,lw=0)
    if ncheck[0] != 0:
        pylab.scatter(peakdata[idx==ncheck[0],2],peakdata[idx==ncheck[0],3],c=cTG,s=dotsize+10,lw=0)
    if ncheck[1] != 0:
        pylab.scatter(peakdata[idx==ncheck[1],2],peakdata[idx==ncheck[1],3],c=cCT,s=dotsize-30,lw=0)
    if ncheck[2] != 0:
        pylab.scatter(peakdata[idx==ncheck[2],2],peakdata[idx==ncheck[2],3],c=cTS1,s=dotsize,lw=0)
    if ncheck[3] != 0:
        pylab.scatter(peakdata[idx==ncheck[3],2],peakdata[idx==ncheck[3],3],c=cTS2,s=dotsize,lw=0)
    pylab.xlim(-95,275)
    pylab.ylim(-95,275)
    pylab.xticks([])
    pylab.yticks([])
    pylab.xlabel('$\Delta_{A}$')
    pylab.ylabel('$\Delta_{B}$')
    # F
    fig.text(0.64,0.5,'F',fontsize=fl_size,ha='right',va='bottom')
    ax = pylab.subplot(439)
    icount = []
    for j in range(int(max(idx))+1):
        icount.append(sum(idx==j))
    print icount
    pylab.bar(range(2),[icount[1],icount[2]],color=[cTG,cCT])
    pylab.xticks(pylab.arange(0,4,1),('TG','CT','TS1','TS2'))
    pylab.xlim(-0.6,3.6)
    pylab.ylabel('Number')
    pylab.ylim(0,520)
    pylab.yticks(range(0,513,256))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    # G
    fig.text(0.05,0.27,'G',fontsize=fl_size,ha='right',va='bottom')
    pylab.subplot(4,3,10)
    J_p,J_m,g_wt,g_tt,wmIno0,ctIno0,std_J,std_t,std_I,std_sig = para_homogeneous_scene2()
    fname = 'jp'+str(J_p)+'jm'+str(-J_m)+'gwt'+str(g_wt)+'gtt'+str(g_tt)+'I'+str(wmIno0)+'std'+str(std_J)+'seednet'+str(seednet)+'seedrun'+str(seedrun)
    tuningE0, tuningE1, tuningL0, tuningL1 = tuning_curves(path=datapath+fname)
    peakdata = peak_tuning(tuningE0, tuningE1, tuningL0, tuningL1)
    idx,ncount,ncheck,uc,nuc = clustering(peakdata,eps=eps_sim)
    uc.append(0)
    pylab.plot([-90,270],[-90,270],'k--')
    for j in range(len(uc)):
        pylab.scatter(peakdata[idx==uc[j],0],peakdata[idx==uc[j],1],color='0.5',s=dotsize,lw=0)
    if ncheck[0] != 0:
        pylab.scatter(peakdata[idx==ncheck[0],0],peakdata[idx==ncheck[0],1],c=cTG,s=dotsize,lw=0)
    if ncheck[1] != 0:
        pylab.scatter(peakdata[idx==ncheck[1],0],peakdata[idx==ncheck[1],1],c=cCT,s=dotsize,lw=0)
    if ncheck[2] != 0:
        pylab.scatter(peakdata[idx==ncheck[2],0],peakdata[idx==ncheck[2],1],c=cTS1,s=dotsize+10,lw=0)
    if ncheck[3] != 0:
        pylab.scatter(peakdata[idx==ncheck[3],0],peakdata[idx==ncheck[3],1],c=cTS2,s=dotsize-30,lw=0)
    pylab.xlim(-95,275)
    pylab.ylim(-95,275)
    pylab.xticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.yticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.xlabel('$\Delta_{Early}$')
    pylab.ylabel('$\Delta_{Late}$')
    pylab.subplot(4,3,11)
    pylab.plot([-90,270],[-90,270],'k--')
    pylab.scatter(peakdata[idx==0,2],peakdata[idx==0,3],color='0.5',s=dotsize,lw=0)
    pylab.scatter(peakdata[idx==1,2],peakdata[idx==1,3],c=cTS2,s=dotsize,lw=0)
    pylab.scatter(peakdata[idx==2,2],peakdata[idx==2,3],c=cTS1,s=dotsize,lw=0)
    pylab.scatter(peakdata[idx==3,2],peakdata[idx==3,3],c=cCT,s=dotsize,lw=0)
    pylab.xlim(-95,275)
    pylab.ylim(-95,275)
    pylab.xticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.yticks([])
    pylab.xlabel('$\Delta_{A}$')
    pylab.ylabel('$\Delta_{B}$')
    # H
    fig.text(0.64,0.27,'H',fontsize=fl_size,ha='right',va='bottom')
    ax = pylab.subplot(4,3,12)
    icount = []
    for j in range(int(max(idx))+1):
        icount.append(sum(idx==j))
    print icount
    pylab.bar(range(1,4),[icount[3],icount[2],icount[1]],color=[cCT,cTS1,cTS2])
    pylab.xticks(pylab.arange(0,4,1),('TG','CT','TS1','TS2'))
    pylab.xlim(-0.6,3.6)
    pylab.ylabel('Number')
    pylab.ylim(0,280)
    pylab.yticks(range(0,257,128))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    pylab.subplots_adjust(top=0.95,bottom=0.1,hspace=0.35,wspace=0.35)
    pylab.savefig(figpath+'fig3.pdf')

### Figure 5
def figure5():
    print 'Making figure 5...'
    #eps = 65#55
    #min_samples = 20
    #epilson = 30
    fig = pylab.figure(figsize=[12,7.5])
    # A
    fig.text(0.05,0.9,'A',fontsize=fl_size,ha='right',va='bottom')
    tuningE0 = pylab.loadtxt(matlabpath+'tuningA0.txt')
    tuningE1 = pylab.loadtxt(matlabpath+'tuningB0.txt')
    tuningL0 = pylab.loadtxt(matlabpath+'tuningA400.txt')
    tuningL1 = pylab.loadtxt(matlabpath+'tuningB400.txt')
    peakdata = peak_tuning(tuningE0, tuningE1, tuningL0, tuningL1)
    dipp = pylab.loadtxt(matlabpath+'dipp_selected.txt')
    idx,ncount,ncheck,uc,nuc = clustering(peakdata)
    uc.append(0)
    fpeak = open(datapath+'peak_selected.txt','w')
    for j in range(idx.size):
        for k in range(4):
            fpeak.write(str(peakdata[j,k])+'\t')
        fpeak.write('\n')
    fpeak.close()
    fidx = open(datapath+'idx_selected.txt','w')
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
    pylab.subplot(231)
    pylab.scatter(peakdata[idx==0,0],peakdata[idx==0,1],color='0.5',s=30,lw=0)
    for j in range(len(uc)):
        pylab.scatter(peakdata[idx==uc[j],0],peakdata[idx==uc[j],1],color='0.5',s=30,lw=0)
    if ncheck[0] != 0:
        pylab.scatter(peakdata[idx==ncheck[0],0],peakdata[idx==ncheck[0],1],c=cTG,s=30,lw=0)
    if ncheck[1] != 0:
        pylab.scatter(peakdata[idx==ncheck[1],0],peakdata[idx==ncheck[1],1],c=cCT,s=30,lw=0)
    if ncheck[2] != 0:
        pylab.scatter(peakdata[idx==ncheck[2],0],peakdata[idx==ncheck[2],1],c=cTS1,s=30,lw=0)
    if ncheck[3] != 0:
        pylab.scatter(peakdata[idx==ncheck[3],0],peakdata[idx==ncheck[3],1],c=cTS2,s=30,lw=0)
    pylab.plot([-90,270],[-90,270],'k--')
    pylab.scatter(peakdata[52,0],peakdata[52,1],marker='.',c='k',s=5,lw=5)
    pylab.scatter(peakdata[52,0],peakdata[52,1],marker='o',c=cTG,s=1,lw=1)
    pylab.scatter(peakdata[476,0],peakdata[476,1],marker='.',c='k',s=5,lw=5)
    pylab.scatter(peakdata[476,0],peakdata[476,1],marker='o',c=cCT,s=1,lw=1)
    pylab.scatter(peakdata[655,0],peakdata[655,1],marker='.',c='k',s=5,lw=5)
    pylab.scatter(peakdata[655,0],peakdata[655,1],marker='o',c=cTS1,s=1,lw=1)
    pylab.xlim(-95,275)
    pylab.ylim(-95,275)
    pylab.xticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.yticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.xlabel('$\Delta_{Early}$')
    pylab.ylabel('$\Delta_{Late}$')
    pylab.subplot(232)
    pylab.scatter(peakdata[idx==0,2],peakdata[idx==0,3],color='0.5',s=30,lw=0)
    for j in range(len(uc)):
        pylab.scatter(peakdata[idx==uc[j],2],peakdata[idx==uc[j],3],color='0.5',s=30,lw=0)
    if ncheck[0] != 0:
        pylab.scatter(peakdata[idx==ncheck[0],2],peakdata[idx==ncheck[0],3],c=cTG,s=30,lw=0)
    if ncheck[1] != 0:
        pylab.scatter(peakdata[idx==ncheck[1],2],peakdata[idx==ncheck[1],3],c=cCT,s=30,lw=0)
    if ncheck[2] != 0:
        pylab.scatter(peakdata[idx==ncheck[2],2],peakdata[idx==ncheck[2],3],c=cTS1,s=30,lw=0)
    if ncheck[3] != 0:
        pylab.scatter(peakdata[idx==ncheck[3],2],peakdata[idx==ncheck[3],3],c=cTS2,s=30,lw=0)
    pylab.plot([-90,270],[-90,270],'k--')
    pylab.scatter(peakdata[52,2],peakdata[52,3],marker='.',c='k',s=5,lw=5)
    pylab.scatter(peakdata[52,2],peakdata[52,3],marker='o',c=cTG,s=1,lw=1)
    pylab.scatter(peakdata[476,2],peakdata[476,3],marker='.',c='k',s=5,lw=5)
    pylab.scatter(peakdata[476,2],peakdata[476,3],marker='o',c=cCT,s=1,lw=1)
    pylab.scatter(peakdata[655,2],peakdata[655,3],marker='.',c='k',s=5,lw=5)
    pylab.scatter(peakdata[655,2],peakdata[655,3],marker='o',c=cTS1,s=1,lw=1)
    pylab.xlim(-95,275)
    pylab.ylim(-95,275)
    pylab.xticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.yticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.xlabel('$\Delta_{A}$')
    pylab.ylabel('$\Delta_{B}$')
    # B
    fig.text(0.65,0.9,'B',fontsize=fl_size,ha='right',va='bottom')
    ax = pylab.subplot(4,6,5)
    pylab.hist(peakdata[:,0],pylab.arange(-90,271,22.5),color='k')
    #pylab.xticks(range(0,271,180),['0$^{\circ}$','180$^{\circ}$'])
    pylab.xticks([])
    pylab.yticks(range(0,160,100))
    pylab.xlim(-90,270)
    pylab.ylim(0,180)
    pylab.xlabel('$\Delta_{Early}$')
    if dipp[0,1] == 0.:
        pylab.text(-75,150,'$p <$ 10$^{-5}$')
    else:
        pylab.text(-75,150,'%s %.3f' % ('$p = $',dipp[0,1]))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax = pylab.subplot(4,6,6)
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
    ax = pylab.subplot(4,6,11)
    pylab.hist(peakdata[:,2],pylab.arange(-90,271,22.5),color='k')
    pylab.xticks(range(0,271,180),['0$^{\circ}$','180$^{\circ}$'])
    pylab.yticks(range(0,160,100))
    pylab.xlim(-90,270)
    pylab.ylim(0,180)
    pylab.xlabel('$\Delta_{A}$')
    if dipp[2,1] == 0.:
        pylab.text(-75,150,'$p <$ 10$^{-5}$')
    else:
        pylab.text(-75,150,'%s %.3f' % ('$p = $',dipp[2,1]))
    pylab.text(-240,250,'Number',rotation='vertical')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax = pylab.subplot(4,6,12)
    pylab.hist(peakdata[:,3],pylab.arange(-90,271,22.5),color='k')
    pylab.xticks(range(0,271,180),['0$^{\circ}$','180$^{\circ}$'])
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
    # C
    fig.text(0.05,0.45,'C',fontsize=fl_size,ha='right',va='bottom')
    pylab.subplot(234)
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
    pylab.plot([-90,-30],[-90,-30],'k--')
    pylab.plot([25,150],[25,150],'k--')
    pylab.plot([200,270],[200,270],'k--')
    pylab.xlim(-95,275)
    pylab.ylim(-95,275)
    pylab.xticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.yticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.xlabel('$\Delta_{Early}$')
    pylab.ylabel('$\Delta_{Late}$')
    pylab.subplot(235)
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
    # D
    fig.text(0.65,0.45,'D',fontsize=fl_size,ha='right',va='bottom')
    ax = pylab.subplot(236)
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
    pylab.subplots_adjust(bottom=0.1,hspace=0.35,wspace=0.45)
    pylab.savefig(figpath+'fig5.pdf')

### Figure 6 joint
def figure6():
    print 'Making figure 6...'
    fig = pylab.figure(figsize=[12,11])
    sID = pylab.loadtxt(matlabpath+'selectedID.txt',dtype='int')-1
    with open(matlabpath+'brainarea_all.txt','r') as tempfile:
        ba = tempfile.read().replace('\n', '')
    numba = []
    for j in range(len(ba)):
        if ba[j] == 'v':
            numba.append(0)
        else:
            numba.append(1)
    numba = np.array(numba)
    nA = sum(numba==0) # 561 #
    nB = sum(numba==1) # 521 #
    numba = np.take(numba,sID)
    ncount = np.zeros((2,4))
    # A
    selected_peakdata = pylab.loadtxt(datapath+'peak_selected.txt')
    selected_idx = pylab.loadtxt(datapath+'idx_selected.txt')
    peakdata = selected_peakdata[numba==0,:]
    idx = selected_idx[numba==0]
    for j in range(4):
        ncount[0,j] = sum(idx==j+1)
    fig.text(0.05,0.9,'A',fontsize=fl_size,ha='right',va='bottom')
    pylab.subplot(331)
    pylab.scatter(peakdata[idx==0,0],peakdata[idx==0,1],color='0.5',s=30,lw=0)
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
    pylab.subplot(332)
    pylab.scatter(peakdata[idx==0,2],peakdata[idx==0,3],color='0.5',s=30,lw=0)
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
    fig.text(0.65,0.9,'B',fontsize=fl_size,ha='right',va='bottom')
    ax = pylab.subplot(6,6,5)
    dipp = pylab.loadtxt(matlabpath+'dipp_selectedv.txt')
    pylab.hist(peakdata[:,0],pylab.arange(-90,271,22.5),color='k')
    #pylab.xticks(range(0,271,180),['0$^{\circ}$','180$^{\circ}$'])
    pylab.xticks([])
    pylab.yticks(range(0,160,50))
    pylab.xlim(-90,270)
    pylab.ylim(0,90)
    pylab.xlabel('$\Delta_{Early}$')
    if dipp[0,1] == 0.:
        pylab.text(30,70,'$p <$ 10$^{-5}$')
    else:
        pylab.text(30,70,'%s %.3f' % ('$p = $',dipp[0,1]))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax = pylab.subplot(6,6,6)
    pylab.hist(peakdata[:,1],pylab.arange(-90,271,22.5),color='k')
    #pylab.xticks(range(0,271,180),['0$^{\circ}$','180$^{\circ}$'])
    pylab.xticks([])
    pylab.yticks(range(0,160,50))
    pylab.xlim(-90,270)
    pylab.ylim(0,90)
    pylab.xlabel('$\Delta_{Late}$')
    if dipp[1,1] == 0.:
        pylab.text(-70,80,'$p <$ 10$^{-5}$')
    else:
        pylab.text(-70,70,'%s %.3f' % ('$p = $',dipp[1,1]))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax = pylab.subplot(6,6,11)
    pylab.hist(peakdata[:,2],pylab.arange(-90,271,22.5),color='k')
    #pylab.xticks(range(0,271,180),['0$^{\circ}$','180$^{\circ}$'])
    pylab.xticks([])
    pylab.yticks(range(0,160,50))
    pylab.xlim(-90,270)
    pylab.ylim(0,90)
    pylab.xlabel('$\Delta_{A}$')
    # pylab.title('%s %.3f %s %.3f' % ('dip=',dipp[1,0],'; p=',dipp[1,1]));
    if dipp[2,1] == 0.:
        pylab.text(-70,70,'$p <$ 10$^{-5}$')
    else:
        pylab.text(-70,70,'%s %.3f' % ('$p = $',dipp[2,1]))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax = pylab.subplot(6,6,12)
    pylab.hist(peakdata[:,3],pylab.arange(-90,271,22.5),color='k')
    #pylab.xticks(range(0,271,180),['0$^{\circ}$','180$^{\circ}$'])
    pylab.xticks([])
    pylab.yticks(range(0,160,50))
    pylab.xlim(-90,270)
    pylab.ylim(0,90)
    pylab.xlabel('$\Delta_{B}$')
    if dipp[3,1] == 0.:
        pylab.text(-70,75,'$p <$ 10$^{-5}$')
    else:
        pylab.text(-70,75,'%s %.3f' % ('$p = $',dipp[3,1]))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    # C
    fig.text(0.05,0.6,'C',fontsize=fl_size,ha='right',va='bottom')
    peakdata = selected_peakdata[numba==1,:]
    idx = selected_idx[numba==1]
    idx0 = np.ones(idx.size)
    for j in range(4):
        ncount[1,j] = sum(idx==j+1)
    print ncount
    pylab.subplot(334)
    pylab.scatter(peakdata[idx==0,0],peakdata[idx==0,1],color='0.5',s=30,lw=0)
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
    pylab.subplot(335)
    pylab.scatter(peakdata[idx==0,2],peakdata[idx==0,3],color='0.5',s=30,lw=0)
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
    # D
    fig.text(0.65,0.6,'D',fontsize=fl_size,ha='right',va='bottom')
    ax = pylab.subplot(6,6,17)
    dipp = pylab.loadtxt(matlabpath+'dipp_selectedd.txt')
    pylab.hist(peakdata[:,0],pylab.arange(-90,271,22.5),color='k')
    #pylab.xticks(range(0,271,180),['0$^{\circ}$','180$^{\circ}$'])
    pylab.xticks([])
    pylab.yticks(range(0,160,50))
    pylab.xlim(-90,270)
    pylab.ylim(0,70)
    pylab.xlabel('$\Delta_{Early}$')
    pylab.text(-70,55,r'$p = 7\times 10^{-5}$')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax = pylab.subplot(6,6,18)
    pylab.hist(peakdata[:,1],pylab.arange(-90,271,22.5),color='k')
    #pylab.xticks(range(0,271,180),['0$^{\circ}$','180$^{\circ}$'])
    pylab.xticks([])
    pylab.yticks(range(0,160,50))
    pylab.xlim(-90,270)
    pylab.ylim(0,70)
    pylab.xlabel('$\Delta_{Late}$')
    pylab.text(-70,60,'%s %.3f' % ('$p = $',dipp[1,1]))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax = pylab.subplot(6,6,23)
    pylab.hist(peakdata[:,2],pylab.arange(-90,271,22.5),color='k')
    pylab.xticks(range(0,271,180),['0$^{\circ}$','180$^{\circ}$'])
    pylab.yticks(range(0,160,50))
    pylab.xlim(-90,270)
    pylab.ylim(0,70)
    pylab.xlabel('$\Delta_{A}$')
    pylab.text(-70,60,'%s %.3f' % ('$p = $',dipp[2,1]))
    pylab.text(-240,120,'Number',rotation='vertical')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax = pylab.subplot(6,6,24)
    pylab.hist(peakdata[:,3],pylab.arange(-90,271,22.5),color='k')
    pylab.xticks(range(0,271,180),['0$^{\circ}$','180$^{\circ}$'])
    pylab.yticks(range(0,160,50))
    pylab.xlim(-90,270)
    pylab.ylim(0,70)
    pylab.xlabel('$\Delta_{B}$')
    # pylab.title('%s %.3f %s %.3f' % ('dip=',dipp[1,0],'; p=',dipp[1,1]));
    pylab.text(-70,60,'%s %.3f' % ('$p = $',dipp[3,1]))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    # E
    fig.text(0.2,0.3,'E',fontsize=fl_size,ha='right',va='bottom')
    #nA = 395
    #nB = 270
    print nA, nB
    ax = fig.add_axes([0.25, 0.1, 0.5, 0.2])
    pylab.bar([0,1,2.5,3.5,5,6],[ncount[0,0]/float(nA),ncount[1,0]/float(nB),ncount[0,1]/float(nA),ncount[1,1]/float(nB),ncount[0,2]/float(nA),ncount[1,2]/float(nB)],color=[cTG,cTG,cCT,cCT,cTS1,cTS1])
    pylab.bar([5,6],[ncount[0,3]/float(nA),ncount[1,3]/float(nB)],bottom=[ncount[0,2]/float(nA),ncount[1,2]/float(nB)],color=[cTS2,cTS2])
    pval = compp(ncount[0,0]/float(nA),ncount[1,0]/float(nB),nA,nB)
    print 'TG: '+str(pval)
    if pval < 0.05:
        txt = '*'
    else:
        txt = 'NS'
    sigdiffh = 1.1*max(ncount[0,0]/float(nA),ncount[1,0]/float(nB))
    pylab.plot([0.,0.,1.,1.],[sigdiffh,sigdiffh+0.02,sigdiffh+0.02,sigdiffh],'k',lw=1.5)
    pylab.text(0.5,sigdiffh+0.02,txt,ha='center',va='bottom')
    pval = compp(ncount[0,1]/float(nA),ncount[1,1]/float(nB),nA,nB)
    print 'CT: '+str(pval)
    if pval < 0.05:
        txt = '*'
    else:
        txt = 'NS'
    sigdiffh = 1.1*max(ncount[0,1]/float(nA),ncount[1,1]/float(nB))
    pylab.plot([2.5,2.5,3.5,3.5],[sigdiffh,sigdiffh+0.02,sigdiffh+0.02,sigdiffh],'k',lw=1.5)
    pylab.text(3.,sigdiffh+0.02,txt,ha='center',va='bottom')
    pval = compp(ncount[0,2]/float(nA)+ncount[0,3]/float(nA),ncount[1,2]/float(nB)+ncount[1,3]/float(nB),nA,nB)
    print 'TS: '+str(pval)
    if pval < 0.05:
        txt = '*'
    else:
        txt = 'NS'
    sigdiffh = 1.06*max(ncount[0,2]/float(nA)+ncount[0,3]/float(nA),ncount[1,2]/float(nB)+ncount[1,3]/float(nB))
    pylab.plot([5.,5.,6.,6.],[sigdiffh,sigdiffh+0.02,sigdiffh+0.02,sigdiffh],'k',lw=1.5)
    pylab.text(5.5,sigdiffh+0.02,txt,ha='center',va='bottom')
    pylab.xticks([0.,1.,2.5,3.5,5,6],('v','d','v','d','v','d'))
    pylab.xlabel('    TG                             CT                             TS')
    """
    pylab.plot([0.4,0.4,1.4,1.4],[0.1,0.12,0.12,0.1],'k',lw=1.5)
    pylab.text(0.9,0.12,'NS',ha='center',va='bottom')
    pylab.plot([2.9,2.9,3.9,3.9],[0.18,0.2,0.2,0.18],'k',lw=1.5)
    pylab.text(3.4,0.2,'NS',ha='center',va='bottom')
    pylab.plot([5.4,5.4,6.4,6.4],[0.26,0.28,0.28,0.26],'k',lw=1.5)
    pylab.text(5.9,0.28,'*',ha='center',va='bottom')
    pylab.xticks([0.4,1.4,2.9,3.9,5.4,6.4],('v','d','v','d','v','d'))
    pylab.xlabel('    TG                             CT                             TS')
    """
    #pylab.xticks([0.5,3.5,6.5],('\nTG','\nCT','\nTS'))
    pylab.xlim(-0.9,6.6)
    pylab.ylabel('Fraction')
    pylab.yticks([0,0.1,0.2])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    pylab.subplots_adjust(bottom=0.1,hspace=0.35,wspace=0.45)
    pylab.savefig(figpath+'fig6.pdf')

### Figure 4
def figure4():
    font_size = 14
    mpl.rcParams['axes.titlesize'] = font_size+2
    mpl.rcParams['xtick.labelsize'] = font_size
    mpl.rcParams['ytick.labelsize'] = font_size
    mpl.rcParams['axes.labelsize'] = font_size
    mpl.rcParams['legend.fontsize'] = font_size-2
    mpl.rcParams['font.size'] = font_size
    print 'Making figure 4...'
    print 'Note that from figure 5, idx = ncheck refer to TG, CT, TS1 and TS2, respectively.'
    idx = pylab.loadtxt(datapath+'idx_selected.txt')
    #itg = pylab.argwhere(idx==ncheck[0])
    #ict = pylab.argwhere(idx==ncheck[1])
    #its1 = pylab.argwhere(idx==ncheck[2])
    #its2 = pylab.argwhere(idx==ncheck[3])
    peak0 = pylab.empty(11)
    peak1 = pylab.empty(11)
    #selectedID_org = [13,40,66,71,62,73]
    #selectedID = [itg[12][0],itg[39][0],ict[65][0],ict[70][0],its1[61][0],its2[72][0]]
    selectedID = [52,619,476,511,655,641]
    tune0 = []
    tune1 = []
    fig = pylab.figure(figsize=[12,12])
    # A
    fig.text(0.05,0.9,'A',fontsize=fl_size,ha='right',va='bottom')
    # B
    fig.text(0.05,0.61,'B',fontsize=fl_size,ha='right',va='bottom')
    # C
    fig.text(0.05,0.34,'C',fontsize=fl_size,ha='right',va='bottom')
    for j in range(3):
        tune = pylab.loadtxt(matlabpath+'tuningA'+str(j*200)+'.txt')
        tune0.append(tune[selectedID,:])
        tune = pylab.loadtxt(matlabpath+'tuningB'+str(j*200)+'.txt')
        tune1.append(tune[selectedID,:])
    print tune0[2][2],tune1[0][5],tune1[1][1]
    for k in range(6):
        id = selectedID[k]
        for j in range(11):
            tune = pylab.loadtxt(matlabpath+'tuningA'+str(j*50)+'.txt')
            peak0[j] = prefdir(tune[id,:])
            tune = pylab.loadtxt(matlabpath+'tuningB'+str(j*50)+'.txt')
            peak1[j] = prefdir(tune[id,:])
        ax = pylab.subplot(6,4,k*4+4)
        if k == 1:
            pylab.plot(pylab.arange(100,601,50),peak0,c=cA,lw=3,label='A chosen')
        else:
            pylab.plot(pylab.arange(100,601,50),peak0,c=cA,lw=2,label='A chosen')
        pylab.plot(pylab.arange(100,601,50),peak1,c=cB,lw=2,label='B chosen')
        if k == 0:
            pylab.title('Peak')
            pylab.legend(loc=4,frameon=False)
        pylab.xlim(0,700)
        pylab.ylim(0,363)
        pylab.xticks(range(0,601,300))
        #pylab.yticks(range(0,360,90),['0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
        pylab.yticks(range(0,361,180),['0$^{\circ}$','180$^{\circ}$','360$^{\circ}$'])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        if k == 5:
            pylab.xlabel('Time from target onset (ms)')
            pylab.ylabel('Peak')
        for j in range(3):
            ax = pylab.subplot(6,4,k*4+1+j)
            fc1 = sp.interp1d(range(-180,360+180,45),np.append(tune0[j][k][-4:],np.append(tune0[j][k],tune0[j][k][:4])),kind='cubic')
            fc2 = sp.interp1d(range(-180,360+180,45),np.append(tune1[j][k][-4:],np.append(tune1[j][k],tune1[j][k][:4])),kind='cubic')
            pylab.plot(range(0,360,45),tune0[j][k],'o',c=cA)
            pylab.plot(range(0,360,15),fc1(range(0,361,15))[:-1],c=cA,lw=2)
            pylab.plot(range(0,360,45),tune1[j][k],'o',c=cB)
            pylab.plot(range(0,360,15),fc2(range(0,361,15))[:-1],c=cB,lw=2)
            if k == 1:
                pylab.plot(prefdir(tune0[j][k])*pylab.ones(2),[-20,120],lw=3,c=cA)
            else:
                pylab.plot(prefdir(tune0[j][k])*pylab.ones(2),[-20,120],c=cA)
            pylab.plot(prefdir(tune1[j][k])*pylab.ones(2),[-20,120],c=cB)
            pylab.xticks(range(0,360,90),['0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
            pylab.xlim(-20,335)
            if max(np.append(tune0[j][k],tune1[j][k]))-min(np.append(tune0[j][k],tune1[j][k]))>70:
                pylab.yticks(range(0,100,50))
                pylab.ylim(min(np.append(tune0[j][k],tune1[j][k]))-20,max(np.append(tune0[j][k],tune1[j][k]))+20)
            elif max(np.append(tune0[j][k],tune1[j][k]))-min(np.append(tune0[j][k],tune1[j][k]))>35:
                pylab.yticks(range(0,100,20))
                pylab.ylim(min(np.append(tune0[j][k],tune1[j][k]))-10,max(np.append(tune0[j][k],tune1[j][k]))+10)
            else:
                pylab.yticks(range(0,100,10))
                pylab.ylim(min(np.append(tune0[j][k],tune1[j][k]))-5,max(np.append(tune0[j][k],tune1[j][k]))+5)
            if k == 5:
                pylab.xlabel('Target A')
                if j == 0:
                    pylab.ylabel('Rate (Hz)')
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')
            if k == 0:
                if j == 0:
                    phase = 'Early'
                elif j == 1:
                    phase = 'Mid'
                elif j == 2:
                    phase = 'Late'
                pylab.title(phase+'\n'+str(200*j)+'-'+str(200*j+200)+' ms')
    pylab.subplots_adjust(hspace=0.2,wspace=0.3)
    pylab.savefig(figpath+'fig4.pdf')

#def figure7(bamp,Jp_user,Jm_user,wmI_user,std_user,gst_user,gtt_user=None,gwt_user=None):
def figure7():
    print 'Making figure 7...'
    eps = 65
    min_samples = 20
    epilson = 30
    J_p,J_m,g_wt,g_tt,wmIno0,ctIno0,std_J,std_t,std_I,std_sig = para_heterogeneous()
    """
    print "User-defined Jp and Jm"
    J_p = Jp_user
    J_m = Jm_user
    wmIno0 = wmI_user
    ctIno0 = wmI_user
    std_J = std_user
    if gtt_user != None:
        g_tt = gtt_user
    if gwt_user != None:
        g_wt = gwt_user
        """
    fig = pylab.figure(figsize=[12,15])
    pylab.subplots_adjust(top=0.9,bottom=0.1,hspace=0.35,wspace=0.35)
    rng = np.random.RandomState(7522)
    # A
    fig.text(0.05,0.9,'A',fontsize=fl_size,ha='right',va='bottom')
    ax = pylab.subplot(941)
    gau = np.zeros(n+1)
    for k in range(n+1):
        gau[k] = gauss(degdiff(-n/2,k))
    #pylab.plot(range(-n/2,n/2+1),J_m+J_tt_p*gau2+std_J*np.append(rn,rn[0]),'0.8',marker='o',ms=2,lw=0,label='within')
    pylab.fill_between(range(-n/2,n/2+1),J_m/2.+J_p*(1-g_tt/2.)*gau-std_J,J_m/2.+J_p*(1-g_tt/2.)*gau+std_J,alpha=0.3,edgecolor='w',facecolor='m')
    rn = rng.randn(n)
    pylab.plot(range(-n/2,n/2+1),J_m/2.+J_p*(1-g_tt/2.)*gau+std_J*np.append(rn,rn[0]),'m',marker='o',ms=2,lw=0,label='within')
    pylab.plot(range(-n/2,n/2+1),J_m/2.+J_p*(1-g_tt/2.)*gau,'m',lw=2,label='within')
    #print min(J_m+J_tt_p*gau2+std_J*np.append(rn,rn[0])),max(J_m+J_tt_p*gau2+std_J*np.append(rn,rn[0]))
    pylab.plot([0,0],[-8,10],'k')
    pylab.plot([-n/2,n/2+1],[0,0],'k')
    #pylab.plot([n/2+1,n/2+1],[-0.4+0.05,2.1-0.08],'k',linewidth=2)
    pylab.text(-150,6.5,'within')
    #pylab.text(n/2-4.8,0.55,'^')
    #pylab.text(n/2+5,0.6,'$\sigma$')
    #pylab.xticks([-n/2,0,n/2],['-180$^{\circ}$','0$^{\circ}$','180$^{\circ}$'])
    pylab.xticks([])
    pylab.yticks([0])
    pylab.xlim(-n/2,n/2+1)
    pylab.ylim(-8,8)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    #ax.xaxis.set_ticks_position('bottom')
    ax = pylab.subplot(945)
    #pylab.plot(range(-n/2,n/2+1),J_m+g_tt*J_tt_p*gau2+std_g*np.append(rn,rn[0]),'ko',ms=2,label='between')
    pylab.fill_between(range(-n/2,n/2+1),J_m/2.+J_p*g_tt*gau/2.-std_J,J_m/2.+J_p*g_tt*gau/2.+std_J,alpha=0.3,edgecolor='w',facecolor='c')
    rn = rng.randn(n)
    pylab.plot(range(-n/2,n/2+1),J_m/2.+J_p*g_tt*gau/2.+std_J*np.append(rn,rn[0]),'c',marker='o',ms=2,lw=0,label='between')
    pylab.plot(range(-n/2,n/2+1),J_m/2.+J_p*g_tt*gau/2.,'c-',lw=2,label='between')
    pylab.plot([0,0],[-8,10],'k')
    pylab.plot([-n/2,n/2+1],[0,0],'k')
    pylab.text(-150,6.5,'between')
    pylab.xticks([-n/2,0,n/2],['-180$^{\circ}$','0$^{\circ}$','180$^{\circ}$'])
    pylab.yticks([0])
    pylab.xlim(-n/2,n/2+1)
    pylab.ylim(-8,8)
    pylab.xlabel(r'$\Delta\theta$')
    pylab.ylabel('Strength')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    # B
    fig.text(0.33,0.9,'B',fontsize=fl_size,ha='right',va='bottom')
    pylab.subplot(13,3,2)
    #fname = 'input2pk_sig1_ddir180jp2.32jm0.8gwt0.03gtt0.9b0.0I0.3297std2.0seednet101seedrun51'
    fname = 'jp'+str(J_p)+'jm'+str(-J_m)+'gwt'+str(g_wt)+'gtt'+str(g_tt)+'I'+str(wmIno0)+'std'+str(std_J)+'seednet'+str(seednet)+'seedrun'+str(seedrun)
    #fname = 'input2pk_jp'+str(J_p)+'jm'+str(-J_m)+'gwt'+str(g_wt)+'gtt'+str(g_tt)+'b'+str(bamp)+'I'+str(wmIno0)+'std'+str(std_J)+'seednet'+str(seednet)+'seedrun'+str(seedrun)
    #fname = 'input2pk_jp'+str(J_p)+'jm'+str(-J_m)+'gwt'+str(g_wt)+'gst'+str(gst_user)+'gtt'+str(g_tt)+'b'+str(bamp)+'I'+str(wmIno0)+'std'+str(std_J)+'seednet'+str(seednet)+'seedrun'+str(seedrun)
    print 'Please check: ',fname
    item = pylab.loadtxt(datapath+fname+'/lpfc_t1_net_'+str(seednet)+'_run_'+str(seedrun)+'_CJ_0_d_90_.txt')
    pylab.pcolor(item.T*1000,rasterized=True,cmap='jet')
    pylab.clim(0,40)
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
    pylab.clim(0,40)
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
    pylab.clim(0,40)
    cb = pylab.colorbar(ticks=[0,20,40])
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
    pylab.clim(0,40)
    cbar = pylab.colorbar(ticks=[0,20,40])
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
    pylab.clim(0,40)
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
    pylab.clim(0,40)
    cb = pylab.colorbar(ticks=[])
    cb.remove()
    pylab.yticks([])
    #pylab.yticks(range(0,n+1,n/2),['0$^{\circ}$','180$^{\circ}$','360$^{\circ}$'])
    pylab.xticks(pylab.arange(0,(dur-tgt_start-go)/tp+1,500/tp),['0','500','1000'])
    pylab.xlim(0,targeton/tp)
    pylab.ylim(0,n)
    pylab.ylabel('RO')
    pylab.xlabel('Time from target onset (ms)')
    # D, E, F
    peak0 = pylab.empty(11)
    peak1 = pylab.empty(11)
    selectedID = [56,541,53] # arbitrary, random, to be reset: 4
    tune0 = []
    tune1 = []
    fig.text(0.05,0.65,'D',fontsize=fl_size,ha='right',va='bottom')
    fig.text(0.05,0.54,'E',fontsize=fl_size,ha='right',va='bottom')
    fig.text(0.05,0.43,'F',fontsize=fl_size,ha='right',va='bottom')
    for j in range(3):
        tuning0, tuning1 = tuning_cj(datapath+fname,200*j)
        tune0.append(tuning0[selectedID,:]*1000)
        tune1.append(tuning1[selectedID,:]*1000)
    for k in range(3):
        for j in range(11):
            tuning0, tuning1 = tuning_cj(datapath+fname,50*j)
            peak0[j] = prefdir(tuning0[selectedID[k],:])
            peak1[j] = prefdir(tuning1[selectedID[k],:])
        ax = pylab.subplot(7,4,k*4+12)
        pylab.plot(pylab.arange(100,601,50),peak0,c=cA,lw=2,label='A chosen')
        pylab.plot(pylab.arange(100,601,50),peak1,c=cB,lw=2,label='B chosen')
        if k == 0:
            pylab.legend(loc=1,frameon=False)
        pylab.xlim(0,700)
        pylab.ylim(0,363)
        pylab.xticks(range(0,601,300))
        #pylab.yticks(range(0,360,90),['0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
        pylab.yticks(range(0,361,180),['0$^{\circ}$','180$^{\circ}$','360$^{\circ}$'])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        if k == 0:
            pylab.title('Peak')
        if k == 2:
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
            if max(np.append(tune0[j][k],tune1[j][k]))-min(np.append(tune0[j][k],tune1[j][k]))>70:
                pylab.yticks(range(0,100,50))
                pylab.ylim(min(np.append(tune0[j][k],tune1[j][k]))-20,max(np.append(tune0[j][k],tune1[j][k]))+20)
            elif max(np.append(tune0[j][k],tune1[j][k]))-min(np.append(tune0[j][k],tune1[j][k]))>35:
                pylab.yticks(range(0,100,20))
                pylab.ylim(min(np.append(tune0[j][k],tune1[j][k]))-10,max(np.append(tune0[j][k],tune1[j][k]))+10)
            else:
                pylab.yticks(range(0,100,10))
                pylab.ylim(min(np.append(tune0[j][k],tune1[j][k]))-5,max(np.append(tune0[j][k],tune1[j][k]))+5)
            if k == 2:
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
    ###
    fpeak = open(datapath+'peak_sim.txt','w')
    for j in range(n*3):
        for k in range(4):
            fpeak.write(str(peakdata[j,k])+'\t')
        fpeak.write('\n')
    fpeak.close()
    idx,ncount,ncheck,uc,nuc = clustering(peakdata,eps=eps_sim)
    print 'ncheck' + str(ncheck)
    temp = pylab.histogram(idx[:n],range(6))
    print temp
    temp = pylab.histogram(idx[n:2*n],range(6))
    print temp
    temp = pylab.histogram(idx[-n:],range(6))
    print temp
    fidx = open(datapath+'idx_sim.txt','w')
    for j in range(n*3):
        fidx.write(str(idx[j])+'\n')
    fidx.close()
    fncheck = open(datapath+'ncheck_lenuc_sim.txt','w')
    for j in range(4):
        fncheck.write(str(int(ncheck[j]))+'\n')
    fncheck.write(str(len(uc))+'\n')
    fncheck.close()
    #ax = fig.add_axes([0.125, 0.1, 0.2, 0.3])
    # G
    fig.text(0.05,0.25,'G',fontsize=fl_size,ha='right',va='bottom')
    ax = pylab.subplot(4,3,10)
    #tuningE0, tuningE1, tuningL0, tuningL1 = tuning_curves(path=datapath+fname)
    #peakdata = peak_tuning(tuningE0,tuningE1,tuningL0,tuningL1)
    #idx,ncount,ncheck,uc,nuc = clustering(peakdata,eps=eps_sim)
    uc.append(0)
    print 'Print files disabled'
    """
    fpeak = open(datapath+'peak_selected.txt','w')
    for j in range(idx.size):
        for k in range(4):
            fpeak.write(str(peakdata[j,k])+'\t')
        fpeak.write('\n')
    fpeak.close()
    """
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
    pylab.subplot(4,3,11)
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
    """
    pylab.scatter([0,180,0],[0,180,180],marker='+',c='k',s=100,lw=2,label='Theory')
    pos = cov_ellipse(peakdata[idx==ncheck[0],:2], nstd=1, alpha=alpha, color=cTG)
    pylab.scatter(pos[0],pos[1],marker='x',color='0.5',s=100,lw=2,label='Sim')
    pylab.legend(loc=4,scatterpoints=1,frameon=False)
    pos = cov_ellipse(peakdata[idx==ncheck[1],:2], nstd=1, alpha=alpha, color=cCT)
    pylab.scatter(pos[0],pos[1],marker='x',color='0.5',s=100,lw=2)
    pos = cov_ellipse(peakdata[idx==ncheck[2],:2], nstd=1, alpha=alpha, color=cTS1)
    pylab.scatter(pos[0],pos[1],marker='x',color='0.5',s=100,lw=2)
    pos = cov_ellipse(peakdata[idx==ncheck[3],:2], nstd=1, alpha=alpha, color=cTS2)
    pylab.scatter(pos[0],pos[1],marker='x',color='0.5',s=100,lw=2)
    pylab.plot([-90,270],[-90,270],'k--')
    pylab.xlim(-95,275)
    pylab.ylim(-95,275)
    pylab.xticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.yticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.xlabel('$\Delta_{Early}$')
    pylab.ylabel('$\Delta_{Late}$')
    #ax = fig.add_axes([0.4, 0.1, 0.2, 0.3])
    ax = pylab.subplot(4,3,11)
    pylab.scatter([0,0,180],[0,180,0],marker='+',c='k',s=100,lw=2)
    pos = cov_ellipse(peakdata[idx==ncheck[0],2:], nstd=1, alpha=alpha, color=cTG)
    pylab.scatter(pos[0],pos[1],marker='x',color='0.5',s=100,lw=2)
    pos = cov_ellipse(peakdata[idx==ncheck[1],2:], nstd=1, alpha=alpha, color=cCT)
    pylab.scatter(pos[0],pos[1],marker='x',color='0.5',s=100,lw=2)
    pos = cov_ellipse(peakdata[idx==ncheck[2],2:], nstd=1, alpha=alpha, color=cTS1)
    pylab.scatter(pos[0],pos[1],marker='x',color='0.5',s=100,lw=2)
    pos = cov_ellipse(peakdata[idx==ncheck[3],2:], nstd=1, alpha=alpha, color=cTS2)
    pylab.scatter(pos[0],pos[1],marker='x',color='0.5',s=100,lw=2)
    pylab.plot([-90,270],[-90,270],'k--')
    pylab.xlim(-95,275)
    pylab.ylim(-95,275)
    pylab.xticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.yticks(range(-90,271,90),['-90$^{\circ}$','0$^{\circ}$','90$^{\circ}$','180$^{\circ}$','270$^{\circ}$'])
    pylab.xlabel('$\Delta_{A}$')
    pylab.ylabel('$\Delta_{B}$')
    """
    # H
    fig.text(0.65,0.25,'H',fontsize=fl_size,ha='right',va='bottom')
    ax = pylab.subplot(4,3,12)
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
    #pylab.subplots_adjust(top=0.9,bottom=0.1,hspace=0.3,wspace=0.3)
    pylab.savefig(figpath+'fig7.pdf')

### Figure 8
def figure8():
    print 'Making figure 8...'
    #eps = 65
    #min_samples = 20
    #epilson = 30
    import scipy.interpolate as sp
    fig = pylab.figure(figsize=[10,9])
    # A
    fig.text(0.05,0.9,'A',fontsize=fl_size,ha='right',va='bottom')
    ax = pylab.subplot(221)
    # homogeneous network
    J_p,J_m,g_wt,g_tt,wmIno0,ctIno0,std_J,std_t,std_I,std_sig = para_homogeneous_scene1()
    nTG = np.zeros(11)
    nTS = np.zeros(11)
    for k in range(11):
        g_tt = 0.1*k
        fname = 'jp'+str(J_p)+'jm'+str(-J_m)+'gwt'+str(g_wt)+'gtt'+str(g_tt)+'I'+str(wmIno0)+'std'+str(std_J)+'seednet'+str(seednet)+'seedrun'+str(seedrun)
        tuningE0, tuningE1, tuningL0, tuningL1 = tuning_curves(datapath+fname)
        peakdata = peak_tuning(tuningE0, tuningE1, tuningL0, tuningL1)
        #peakdata = peakdata[:n*2,:]
        idx,ncount,ncheck,uc,nuc = clustering(peakdata,eps=eps_sim,min_samples=min_samples,epilson=epilson)
        nTG[k] = ncount[0]
        nTS[k] = ncount[2] + ncount[3]
    pylab.plot(np.arange(0,1.1,0.1),nTG/(n*2.),'+',color='0.5',ms=10,mew=2,label='TG, hom')
    pylab.plot(np.arange(0,1.1,0.1),nTS/(n*2.),'x',color='0.5',ms=8,mew=2,label='TS, hom')
    J_p,J_m,g_wt,g_tt,wmIno0,ctIno0,std_J,std_t,std_I,std_sig = para_heterogeneous()
    if overwrite or not os.path.exists(datapath+'fig8A.txt'):
        nTG = np.zeros(11,dtype=int)
        nTS = np.zeros(11,dtype=int)
        for k in range(4,11):
            g_tt = 0.1*k
            fname = 'jp'+str(J_p)+'jm'+str(-J_m)+'gwt'+str(g_wt)+'gtt'+str(g_tt)+'I'+str(wmIno0)+'std'+str(std_J)+'seednet'+str(seednet)+'seedrun'+str(seedrun)
            tuningE0, tuningE1, tuningL0, tuningL1 = tuning_curves(datapath+fname)
            peakdata = peak_tuning(tuningE0, tuningE1, tuningL0, tuningL1)
            #peakdata = peakdata[:n*2,:]
            idx,ncount,ncheck,uc,nuc = clustering(peakdata,eps=eps_sim,min_samples=min_samples,epilson=epilson)
            nTG[k] = ncount[0]
            nTS[k] = ncount[2] + ncount[3]
        nTG = nTG[4:]
        nTS = nTS[4:]
        fprop = open(datapath+'fig8A.txt','w')
        for j in range(7):
            fprop.write(str(nTG[j])+'\t'+str(nTS[j]))
            fprop.write('\n')
        fprop.close()
    else:
        temp = np.loadtxt(datapath+'fig8A.txt')
        nTG = temp[:,0]
        nTS = temp[:,1]
    pylab.plot(np.arange(0.4,1.1,0.1),nTG/(n*2.),cTG,marker='+',ms=10,mew=2,lw=2,label='TG, het')
    pylab.plot(np.arange(0.4,1.1,0.1),nTS/(n*2.),cTS1,marker='x',ms=8,mew=2,lw=2,label='TS, het')
    pylab.legend(loc=6,numpoints=1,frameon=False)
    pylab.xlabel(r'$\alpha$')
    pylab.xticks([0,0.5,1],['0','0.5','1'])
    pylab.yticks([0,0.5,1],['0','0.5','1'])
    pylab.ylim(-0.05,1.05)
    pylab.xlim(-0.05,1.05)
    pylab.ylabel('Fraction')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    # B
    fig.text(0.5,0.9,'B',fontsize=fl_size,ha='right',va='bottom')
    ax = pylab.subplot(222)
    J_p,J_m,g_wt,g_tt,wmIno0,ctIno0,std_J,std_t,std_I,std_sig = para_heterogeneous()
    if overwrite or not os.path.exists(datapath+'fig8B.txt'):
        nTG = np.zeros(len(std_J_arr),dtype=int)
        nTS = np.zeros(len(std_J_arr),dtype=int)
        for k in range(len(std_J_arr)):
            fname = 'jp'+str(J_p)+'jm'+str(-J_m)+'gwt'+str(g_wt)+'gtt'+str(g_tt)+'I'+str(wmIno0)+'std'+str(std_J_arr[k])+'seednet'+str(seednet)+'seedrun'+str(seedrun)
            tuningE0, tuningE1, tuningL0, tuningL1 = tuning_curves(datapath+fname)
            peakdata = peak_tuning(tuningE0, tuningE1, tuningL0, tuningL1)
            #peakdata = peakdata[:n*2,:]
            idx,ncount,ncheck,uc,nuc = clustering(peakdata,eps=eps_sim,min_samples=min_samples,epilson=epilson)
            nTG[k] = ncount[0]
            nTS[k] = ncount[2] + ncount[3]
        fprop = open(datapath+'fig8B'+str(eps)+'.txt','w')
        for j in range(len(std_J_arr)):
            fprop.write(str(nTG[j])+'\t'+str(nTS[j]))
            fprop.write('\n')
        fprop.close()
    else:
        temp = np.loadtxt(datapath+'fig8B.txt')
        nTG = temp[:,0]
        nTS = temp[:,1]
    pylab.plot(std_J_arr,nTG/(n*2.),cTG,marker='+',ms=10,mew=2,lw=2,label='TG')
    pylab.plot(std_J_arr,nTS/(n*2.),cTS1,marker='x',ms=8,mew=2,lw=2,label='TS')
    pylab.legend(loc=6,numpoints=1,frameon=False)
    pylab.xticks(range(4))
    pylab.yticks([0,0.5,1],['0','0.5','1'])
    pylab.ylim(-0.05,1.05)
    pylab.xlim(-0.05,3)
    pylab.xlabel(r'$\beta$')
    pylab.ylabel('Fraction')
    #pylab.title('$g_{TT} = 0.45$')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    # C
    fig.text(0.05,0.48,'C',fontsize=fl_size,ha='right',va='bottom')
    ax = pylab.subplot(223)
    if overwrite or not os.path.exists(datapath+'fig8C.txt'):
        prop = np.zeros((11,len(std_J_arr)))
        fprop = open(datapath+'fig8C.txt','w')
        for k in range(4,11):
            g_tt = 0.1*k
            for j in range(len(std_J_arr)-1):
                fname = 'jp'+str(J_p)+'jm'+str(-J_m)+'gwt'+str(g_wt)+'gtt'+str(g_tt)+'I'+str(wmIno0)+'std'+str(std_J_arr[j])+'seednet'+str(seednet)+'seedrun'+str(seedrun)
                tuningE0, tuningE1, tuningL0, tuningL1 = tuning_curves(datapath+fname)
                peakdata = peak_tuning(tuningE0, tuningE1, tuningL0, tuningL1)
                #peakdata = peakdata[:n*2,:]
                idx,ncount,ncheck,uc,nuc = clustering(peakdata,eps=eps_sim,min_samples=min_samples,epilson=epilson)
                nTG = ncount[0]
                nTS = ncount[2] + ncount[3]
                prop[k,j] = float(nTG)/(nTG+nTS)
                fprop.write(str(prop[k,j])+'\t')
                print k,j,prop[k,j]
            fprop.write('\n')
        fprop.close()
        prop = prop[4:,:-1]
    else:
        prop = np.loadtxt(datapath+'fig8C.txt')
    pylab.pcolor(np.arange(0.35,1.1,0.1),std_J_arr-np.append(np.diff(std_J_arr),0.25)/2,prop.T,rasterized=True,cmap='jet')
    # locate 0.28
    fit = np.empty(0)
    for j in range(1,6):
        lb = sum(0.28>=prop[j,:])-1
        fit = np.append(fit,std_J_arr[lb] + np.diff(std_J_arr[lb:lb+2])*(0.28-prop[j,lb])/(prop[j,lb+1]-prop[j,lb]))
    fit = np.append(fit,std_J_arr[-1] + np.diff(std_J_arr[-2:])*(0.28-prop[-1,-1])/(0.40625-prop[-1,-1]))
    pylab.plot(np.arange(0.5,1.01,0.1),fit,'w',lw=3)
    pylab.ylim(-0.25,2.625)
    pylab.xlabel(r'$\alpha$')
    pylab.ylabel(r'$\beta$')
    pylab.xticks(np.arange(0.4,1.1,0.2))
    pylab.yticks(range(3))
    pylab.clim(0,1)
    cbar = pylab.colorbar(ticks=[0,1])
    cbar.set_label('Fraction of TG')
    # D
    fig.text(0.5,0.48,'D',fontsize=fl_size,ha='right',va='bottom')
    ax = pylab.subplot(224)
    hetero_std = np.arange(0.,0.026,0.005)
    if overwrite or not os.path.exists(datapath+'fig8D.txt'):
        fnum = open(datapath+'fig8D.txt','w')
        J_p,J_m,g_wt,g_tt,wmIno0,ctIno0,std_J,std_t,std_I,std_sig = para_heterogeneous_neuron()
        nTG = np.zeros(hetero_std.size)
        nTS = np.zeros(hetero_std.size)
        for k in range(hetero_std.size):
            if k == 0:
                fname = 'jp'+str(J_p)+'jm'+str(-J_m)+'gwt'+str(g_wt)+'gtt'+str(g_tt)+'I'+str(wmIno0)+'std'+str(std_J)+'seednet'+str(seednet)+'seedrun'+str(seedrun)
            else:
                fname = 'jp'+str(J_p)+'jm'+str(-J_m)+'gwt'+str(g_wt)+'gtt'+str(g_tt)+'I'+str(wmIno0)+'std'+str(std_J)+'stdI'+str(hetero_std[k])+'stds'+str(hetero_std[k])+'seednet'+str(seednet)+'seedrun'+str(seedrun)
            tuningE0, tuningE1, tuningL0, tuningL1 = tuning_curves(datapath+fname)
            peakdata = peak_tuning(tuningE0,tuningE1,tuningL0,tuningL1)
            #peakdata = peakdata[:n*2,:]
            idx,ncount,ncheck,uc,nuc = clustering(peakdata)
            nTG[k] = ncount[0]
            nTS[k] = ncount[2] + ncount[3]
            fnum.write(str(nTG[k])+'\t'+str(nTS[k])+'\n')
        fnum.close()
    else:
        temp = np.loadtxt(datapath+'fig8D.txt')
        nTG = temp[:,0]
        nTS = temp[:,1]
    pylab.plot(hetero_std,nTG/(n*2.),cTG,marker='+',ms=10,mew=2,lw=2,label='TG')
    pylab.plot(hetero_std,nTS/(n*2.),cTS1,marker='x',ms=8,mew=2,lw=2,label='TS')
    pylab.xticks(np.arange(0.,0.031,0.01))
    pylab.yticks([0,0.5,1],['0','0.5','1'])
    pylab.ylim(-0.05,1.05)
    pylab.xlim(-0.005,0.036)
    #pylab.xlabel('$\kappa$ (nA)')
    pylab.xlabel('$\kappa$')
    pylab.ylabel('Fraction')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    pylab.subplots_adjust(bottom=0.15,hspace=0.35,wspace=0.35)
    pylab.savefig(figpath+'fig8.pdf')

def morefigure7(bamp,Jp,Jm,wmI,std,gst,gtt,gwt):
    fname = 'input2pk_jp'+str(Jp)+'jm'+str(-Jm)+'gwt'+str(gwt)+'gst'+str(gst)+'gtt'+str(gtt)+'b'+str(bamp)+'I'+str(wmI)+'std'+str(std)+'seednet'+str(seednet)+'seedrun'+str(seedrun)
    print fname
    for CJ in [0,1]:
        for j in range(8):
            if CJ == 0:
                target = j*45
            else:
                target = pylab.mod(j*45+180,360)
            ct = pylab.loadtxt(datapath+fname+'/lpfc_ct_net_'+str(seednet)+'_run_'+str(seedrun)+'_CJ_'+str(CJ)+'_d_'+str(j*45)+'_.txt')
            temp,decodedCT = resultant(ct[-1,:],div=n)
            print 'CT = '+str(target)+'; Decoded CT = '+str(decodedCT)
            if anglediff(target,decodedCT) > 10:
                print 'Please check the result!'
    figure7(bamp,Jp,Jm,wmI,std,gst,gtt,gwt)

#figure1()
#figure2()
#figure3()
#figure5()
#figure4()
#figure6()
figure7()
#figure8()

pylab.show()
