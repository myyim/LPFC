%%% This code is written by Man Yi Yim and finalized on Aug 28, 2017,
%%% with MATLAB R2014b.
%%%
%%% This code generates the necessary data files from monkey recording
%%% data for further analysis and use in Python for the paper:
%%% Yim, Cai & Wang (2019) Transforming the choice outcome to an action
%%% plan in monkey lateral prefrontal cortex: a neural circuit model.
%%%
%%% Data files needed in Cai2014: all X_dirtuningPar.mat (2662 files),
%%% X_psyphycell.mat (2662), X_tuning.mat (2662) and X_bin100.mat 
%%% (2116, in profiles folder)
%%%
%%% Matlab codes needed in codes_matlab: arearead_DT.m, condReadout.m, readsession.m,
%%% sessionlistread.m, targetconvert.m, sessionlist_DT_PFC.m, degreediff.m,
%%% resultant.m, preferdir.m, hartigansdipsigniftest.m, hartigansdiptest.m
%%%
%%% Files generated:
%%% brainarea_all.txt
%%% (which LPFC region, v or d)
%%% dipp_all.txt, dipp_alld.txt, dipp_allv.txt, dipp_selected.txt, 
%%% dipp_selectedd.txt, dipp_selectedv.txt
%%% (dip value and p-value for the peak difference distribution)
%%% peak_significant.txt
%%% (p-value of ranksum test: A0, B0, A400, B400)
%%% pval_CJ.txt
%%% (p-value of CJ encoding in windows 5 and 6, prior to target on)
%%% selectedID.txt
%%% (selected neurons)
%%% tuningallA0.txt, tuningallA400.txt, tuningallB0.txt, tuningallB400.txt,
%%% (8-direction tuning curves when A or B chosen in different windows)
%%% tuningA0,..., tuningA500, tuningB0,..., tuningB500
%%% (same as above but for selected neurons)

%% Definition
disp('Please have all the necessary scripts and data in the MATLAB search path');
disp('Data files will be stored in the matlab folder under the directory containing this script');
clear all;
current_path = mfilename('fullpath');
currdir = fileparts(current_path);
cd(currdir);
addpath('codes_matlab');
addpath('Cai2014');
mkdir('matlab');
cd matlab;

sessions = sessionlistread('sessionlist_DT_PFC.m'); % read the session names
win = {'preoffer','postoffer','latedelay','memdelay','pretarget','posttarget','prego','reactime','prejuice','postjuice'}; % window labels in the data files
pth = 0.01; % threshold p-value
nboot = 100000; % number for bootstrapping in dip test
fba = fopen('brainarea_all.txt','w');
fpeaksig = fopen('peak_significant.txt','w');
fCJ = fopen('pval_CJ.txt','w');
fselect = fopen('selectedID.txt','w');
fallA0 = fopen('tuningallA0.txt','w');
fallA400 = fopen('tuningallA400.txt','w');
fallB0 = fopen('tuningallB0.txt','w');
fallB400 = fopen('tuningallB400.txt','w');
faction = 'w'; % defined for use in eval below
for ts = 0:50:500 % timestamp from target onset
    fname = ['tuningA',num2str(ts),'.txt'];
    eval(['fA',num2str(ts),' = fopen(fname,faction);']);
    fname = ['tuningB',num2str(ts),'.txt'];
    eval(['fB',num2str(ts),' = fopen(fname,faction);']);
end

%% Initiation
ncell = 0;  % total neuron count
nscell = 0; % spatially selective neuron count
pCJ = zeros(1,2);
pdir200allA0 = [];
pdir200allA400 = [];
pdir200allB0 = [];
pdir200allB400 = [];
for ts = 0:50:500 % timestamp from target onset
    eval(['pdir200A',num2str(ts),'=[];']);   % preferred direction for chosen juice A from ts for 200 ms
    eval(['pdir200B',num2str(ts),'=[];']);
end

%% Main code
for isession = 1:size(sessions,1)
    session = sessions{isession};
    readsession;    % read isession info
    
    for icell = 1:size(cells_td,1)     
        ncell = ncell + 1;  % total cell count
        cellname = [num2str(cells_td(icell,1)),num2str(cells_td(icell,2))];
        filename = [session, cellname, '_tuning'];eval(['load ',filename]);
        filename = [session, cellname, '_dirtuningPar'];eval(['load ',filename]);        
        [monkey, brainarea] = arearead_DT(session, cellname);    
        eval(['data = tuning.','AB','.neuract.bytrial.',win{1},';']);
        condTable = condReadout(data); %columns: [offertype, chosenjuice (%1=Achosen, -1=Bchosen), orientation, posA, chosentarget, fr]
        fprintf(fba,'%s\n',brainarea(1));
               
        % Load the rate data after target onset
        filename = [session, cellname, '_bin100'];eval(['load ',filename]);
        data100 = profile.AB.bytrial.offeron;
        r = data100(:,103);
        for its = 107:4:123 % we consider 600 ms for significance test
            r = r + data100(:,its);
        end
        r = r/6;
        
        % 4-way ANOVA (Cai & Padoa-Schioppa, 2014)
        filename = [session, cellname, '_psyphycell'];eval(['load ',filename]);
        relvalue = exp(psyphycell.sigmoidfit(1));
        vA = relvalue*condTable(:,2);
        vB = condTable(:,3);
        cv = vA.*(condTable(:,4)==1) + vB.*(condTable(:,4)==-1);
        cvi = cv > median(cv);        
        cj = condTable(:,4);	%B choice
        ori = condTable(:,5);
        tgt = condTable(:,6);
        ct = condTable(:,7);
        tgtcho = dirhem;
        hemA = ismember(tgt,tgtcho);
        [pval, tbl] = anovan(r, {cvi, ori, hemA, cj}, 'full', 2, {'chosenValue','orientation','hemA','taste'}, 'off');
        pwhere = find(pval<pth);
        neuronsig{ncell} = pwhere;
        
        % p-value of CJ 
        for iwin = 5:6 % prior to target onset, 500 ms each           
            eval(['data = tuning.','AB','.neuract.bytrial.',win{iwin},';']);
            condTable = condReadout(data); %columns: [offertype, chosenjuice (%1=Achosen, -1=Bchosen), orientation, posA, chosentarget, fr]   
            [pval, tbl] = anovan(condTable(:,8), {cvi, ori, hemA, cj}, 'full', 2, {'chosenValue','orientation','hemA','taste'}, 'off');
            pCJ(iwin-4) = pval(4);
        end %iwin
        fprintf(fCJ, '%20.8f \t %20.8f \n',pCJ(1),pCJ(2));
        
        for its = 103:16:119
            % print tuning curves into files
            r = (data100(:,its)+data100(:,its+4))/2; % we consider 200 ms for response tracking
            for k = 1:8
                posA(k) = mean(r(condTable(:,4)==1 & condTable(:,6)==k));
                posB(k) = mean(r(condTable(:,4)==-1 & condTable(:,6)==k));
            end
            if its == 103
                for j = 1:8
                    fprintf(fallA0,'%20.8f \t',posA(j));
                    fprintf(fallB0,'%20.8f \t',posB(j));
                end
                fprintf(fallA0,'\n');
                fprintf(fallB0,'\n');
            elseif its == 119
                for j = 1:8
                    fprintf(fallA400,'%20.8f \t',posA(j));
                    fprintf(fallB400,'%20.8f \t',posB(j));
                end
                fprintf(fallA400,'\n');
                fprintf(fallB400,'\n');
            end
            
            % locate the preferred direction
            eval(['pdir200allA',num2str((its-103)*25),' = [pdir200allA',num2str((its-103)*25),', preferdir(posA)];']);
            eval(['pdir200allB',num2str((its-103)*25),' = [pdir200allB',num2str((its-103)*25),', preferdir(posB)];']);
            
            % Are the two peaks statistically different?
            [temp mpos] = max(posA);
            r1 = r(condTable(:,4)==1 & condTable(:,6)==mpos);
            r2 = r(condTable(:,4)==1 & condTable(:,6)==mod(mpos+3,8)+1);
            pA = ranksum(r1,r2);
            [temp mpos] = max(posB);
            r1 = r(condTable(:,4)==-1 & condTable(:,6)==mpos);
            r2 = r(condTable(:,4)==-1 & condTable(:,6)==mod(mpos+3,8)+1);
            pB = ranksum(r1,r2);
            fprintf(fpeaksig,'%20.8f \t %20.8f \t',pA,pB);
        end % its
        disp(ncell);
        fprintf(fpeaksig,'\n');
        
        % Select spatially-tuned neurons for analysis
        if any(neuronsig{ncell}==2) || any(neuronsig{ncell}==3) || any(neuronsig{ncell}==8) || any(neuronsig{ncell}==10) || any(neuronsig{ncell}==14)
            nscell = nscell + 1;    % selected neuron count
            fprintf(fselect,'%d \n',ncell);
            for its = 103:2:123
                % print the tuning curves into files
                r = (data100(:,its)+data100(:,its+4))/2;             
                for k = 1:8
                    posA(k) = mean(r(condTable(:,4)==1 & condTable(:,6)==k));
                    posB(k) = mean(r(condTable(:,4)==-1 & condTable(:,6)==k));
                end
                if its == 103
                    for j = 1:8
                        fprintf(fA0,'%20.8f \t',posA(j));
                        fprintf(fB0,'%20.8f \t',posB(j));
                    end
                    fprintf(fA0,'\n');
                    fprintf(fB0,'\n');
                elseif its == 105
                    for j = 1:8
                        fprintf(fA50,'%20.8f \t',posA(j));
                        fprintf(fB50,'%20.8f \t',posB(j));
                    end
                    fprintf(fA50,'\n');
                    fprintf(fB50,'\n');
                elseif its == 107
                    for j = 1:8
                        fprintf(fA100,'%20.8f \t',posA(j));
                        fprintf(fB100,'%20.8f \t',posB(j));
                    end
                    fprintf(fA100,'\n');
                    fprintf(fB100,'\n');
                elseif its == 109
                    for j = 1:8
                        fprintf(fA150,'%20.8f \t',posA(j));
                        fprintf(fB150,'%20.8f \t',posB(j));
                    end
                    fprintf(fA150,'\n');
                    fprintf(fB150,'\n');
                elseif its == 111
                    for j = 1:8
                        fprintf(fA200,'%20.8f \t',posA(j));
                        fprintf(fB200,'%20.8f \t',posB(j));
                    end
                    fprintf(fA200,'\n');
                    fprintf(fB200,'\n');
                elseif its == 113
                    for j = 1:8
                        fprintf(fA250,'%20.8f \t',posA(j));
                        fprintf(fB250,'%20.8f \t',posB(j));
                    end
                    fprintf(fA250,'\n');
                    fprintf(fB250,'\n');
                elseif its == 115
                    for j = 1:8
                        fprintf(fA300,'%20.8f \t',posA(j));
                        fprintf(fB300,'%20.8f \t',posB(j));
                    end
                    fprintf(fA300,'\n');
                    fprintf(fB300,'\n');
                elseif its == 117
                    for j = 1:8
                        fprintf(fA350,'%20.8f \t',posA(j));
                        fprintf(fB350,'%20.8f \t',posB(j));
                    end
                    fprintf(fA350,'\n');
                    fprintf(fB350,'\n');
                elseif its == 119
                    for j = 1:8
                        fprintf(fA400,'%20.8f \t',posA(j));
                        fprintf(fB400,'%20.8f \t',posB(j));
                    end
                    fprintf(fA400,'\n');
                    fprintf(fB400,'\n');
                elseif its == 121
                    for j = 1:8
                        fprintf(fA450,'%20.8f \t',posA(j));
                        fprintf(fB450,'%20.8f \t',posB(j));
                    end
                    fprintf(fA450,'\n');
                    fprintf(fB450,'\n');
                elseif its == 123
                    for j = 1:8
                        fprintf(fA500,'%20.8f \t',posA(j));
                        fprintf(fB500,'%20.8f \t',posB(j));
                    end
                    fprintf(fA500,'\n');
                    fprintf(fB500,'\n');
                end
                eval(['pdir200A',num2str((its-103)*25),' = [pdir200A',num2str((its-103)*25),', preferdir(posA)];']);
                eval(['pdir200B',num2str((its-103)*25),' = [pdir200B',num2str((its-103)*25),', preferdir(posB)];']);           
            end % its
        end % if
    end % icell
end % isession

fclose(fba);
fclose(fpeaksig);
fclose(fCJ);
fclose(fselect);
fclose(fallA0);
fclose(fallA400);
fclose(fallB0);
fclose(fallB400);
for ts = 0:50:500
    eval(['fclose(fA',num2str(ts),');']);
    eval(['fclose(fB',num2str(ts),');']);
end

%% dip test for all neurons
idx = isnan(pdir200allA0+pdir200allB0+pdir200allA400+pdir200allB400);
ba = importdata('brainarea_all.txt');

f = fopen('dipp_all.txt','w');
angAe = pdir200allA0(idx==0);
angBe = pdir200allB0(idx==0);
angAl = pdir200allA400(idx==0);
angBl = pdir200allB400(idx==0);
ang = degreediff(angAe,angBe,2);
[dip, p] = hartigansdipsigniftest(ang, nboot);
fprintf(f,'%20.8f %20.8f \n',dip,p);
ang = degreediff(angAl,angBl,2);
[dip, p] = hartigansdipsigniftest(ang, nboot);
fprintf(f,'%20.8f %20.8f \n',dip,p);
ang = degreediff(angAe,angAl,2);
[dip, p] = hartigansdipsigniftest(ang, nboot);
fprintf(f,'%20.8f %20.8f \n',dip,p);
ang = degreediff(angBe,angBl,2);
[dip, p] = hartigansdipsigniftest(ang, nboot);
fprintf(f,'%20.8f %20.8f \n',dip,p);
fclose(f);

% dip test for LPFCv
f = fopen('dipp_allv.txt','w');
region_idx = logical((idx==0).*strcmp(ba,'v')');
angAe = pdir200allA0(region_idx);
angBe = pdir200allB0(region_idx);
angAl = pdir200allA400(region_idx);
angBl = pdir200allB400(region_idx);
ang = degreediff(angAe,angBe,2);
[dip, p] = hartigansdipsigniftest(ang, nboot);
fprintf(f,'%20.8f %20.8f \n',dip,p);
ang = degreediff(angAl,angBl,2);
[dip, p] = hartigansdipsigniftest(ang, nboot);
fprintf(f,'%20.8f %20.8f \n',dip,p);
ang = degreediff(angAe,angAl,2);
[dip, p] = hartigansdipsigniftest(ang, nboot);
fprintf(f,'%20.8f %20.8f \n',dip,p);
ang = degreediff(angBe,angBl,2);
[dip, p] = hartigansdipsigniftest(ang, nboot);
fprintf(f,'%20.8f %20.8f \n',dip,p);
fclose(f);

% dip test for LPFCd
f = fopen('dipp_alld.txt','w');
region_idx = logical((idx==0).*strcmp(ba,'d')');
angAe = pdir200allA0(region_idx);
angBe = pdir200allB0(region_idx);
angAl = pdir200allA400(region_idx);
angBl = pdir200allB400(region_idx);
ang = degreediff(angAe,angBe,2);
[dip, p] = hartigansdipsigniftest(ang, nboot);
fprintf(f,'%20.8f %20.8f \n',dip,p);
ang = degreediff(angAl,angBl,2);
[dip, p] = hartigansdipsigniftest(ang, nboot);
fprintf(f,'%20.8f %20.8f \n',dip,p);
ang = degreediff(angAe,angAl,2);
[dip, p] = hartigansdipsigniftest(ang, nboot);
fprintf(f,'%20.8f %20.8f \n',dip,p);
ang = degreediff(angBe,angBl,2);
[dip, p] = hartigansdipsigniftest(ang, nboot);
fprintf(f,'%20.8f %20.8f \n',dip,p);
fclose(f);

%% dip test for selected neurons
sidx = importdata('selectedID.txt');
idx = isnan(pdir200A0+pdir200B0+pdir200A400+pdir200B400);
ba_selected = ba(sidx);
bav = strcmp(ba_selected,'v');

f = fopen('dipp_selected.txt','w');
angAe = pdir200A0(idx==0);
angBe = pdir200B0(idx==0);
angAl = pdir200A400(idx==0);
angBl = pdir200B400(idx==0);
ang = degreediff(angAe,angBe,2);
[dip, p] = hartigansdipsigniftest(ang, nboot);
fprintf(f,'%20.8f %20.8f \n',dip,p);
ang = degreediff(angAl,angBl,2);
[dip, p] = hartigansdipsigniftest(ang, nboot);
fprintf(f,'%20.8f %20.8f \n',dip,p);
ang = degreediff(angAe,angAl,2);
[dip, p] = hartigansdipsigniftest(ang, nboot);
fprintf(f,'%20.8f %20.8f \n',dip,p);
ang = degreediff(angBe,angBl,2);
[dip, p] = hartigansdipsigniftest(ang, nboot);
fprintf(f,'%20.8f %20.8f \n',dip,p);
fclose(f);

% dip test for LPFCv
f = fopen('dipp_selectedv.txt','w');
region_idx = logical((idx==0).*(bav==1)');
angAe = pdir200A0(region_idx);
angBe = pdir200B0(region_idx);
angAl = pdir200A400(region_idx);
angBl = pdir200B400(region_idx);
ang = degreediff(angAe,angBe,2);
[dip, p] = hartigansdipsigniftest(ang, nboot);
fprintf(f,'%20.8f %20.8f \n',dip,p);
ang = degreediff(angAl,angBl,2);
[dip, p] = hartigansdipsigniftest(ang, nboot);
fprintf(f,'%20.8f %20.8f \n',dip,p);
ang = degreediff(angAe,angAl,2);
[dip, p] = hartigansdipsigniftest(ang, nboot);
fprintf(f,'%20.8f %20.8f \n',dip,p);
ang = degreediff(angBe,angBl,2);
[dip, p] = hartigansdipsigniftest(ang, nboot);
fprintf(f,'%20.8f %20.8f \n',dip,p);
fclose(f);

% dip test for LPFCd
f = fopen('dipp_selectedd.txt','w');
region_idx = logical((idx==0).*(bav==0)');
angAe = pdir200A0(region_idx);
angBe = pdir200B0(region_idx);
angAl = pdir200A400(region_idx);
angBl = pdir200B400(region_idx);
ang = degreediff(angAe,angBe,2);
[dip, p] = hartigansdipsigniftest(ang, nboot);
fprintf(f,'%20.8f %20.8f \n',dip,p);
ang = degreediff(angAl,angBl,2);
[dip, p] = hartigansdipsigniftest(ang, nboot);
fprintf(f,'%20.8f %20.8f \n',dip,p);
ang = degreediff(angAe,angAl,2);
[dip, p] = hartigansdipsigniftest(ang, nboot);
fprintf(f,'%20.8f %20.8f \n',dip,p);
ang = degreediff(angBe,angBl,2);
[dip, p] = hartigansdipsigniftest(ang, nboot);
fprintf(f,'%20.8f %20.8f \n',dip,p);
fclose(f);
