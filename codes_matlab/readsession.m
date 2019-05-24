%
%


clear parsession
filename = ['DT_',session];	eval([filename])

monkey = parsession.monkey;
day = parsession.day;
sessnum = parsession.sessnum;

%dirroot = ['/Experiments/Delay_Target/Data/',...
%	parsession.mnkdir, filesep, day, sessnum, filesep];

%fileroot_ML = ['DT', session(2:end)];
%load(fileroot_ML)
%fileroot_CED = session(1:end-1);
%filenums_CED = parsession.CEDfilenums;

cells_td = parsession.clusters;

%flags
flags.trialstart	= 21;
flags.fixon			= 25;
flags.fixoff		= 27;
flags.offeron		= 30;
flags.offeroff		= 31;
flags.sacctgton		= 35; 
flags.choicemade	= 39;
flags.juice1		= 43;
flags.juice2		= 44;
flags.juice3		= 45;
flags.trialend		= 60;
%
flags.poss_outcomes = [flags.juice1, flags.juice2, flags.juice3];
flags.outcome		= 50;

