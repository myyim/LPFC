function condTable = condReadout(data)
%data columns: trialnum, #A, #B, posA, tgtnumA, chosenpos(left or right)

ii=[];
[~,ii] = sort(data(:,1));
hitdata = data(ii,:);

offertype = hitdata(:,1:3); %trial number
%
chosenjuice = hitdata(:,4).*hitdata(:,6);	%1=A chosen, -1=B chosen;
%
tgtA = hitdata(:,5); %position of targetA
%
orientation = targetconvert(tgtA); %convert postion of targetA to orientation
%
chosentarget = hitdata(:,5);
ind = hitdata(:,4).*hitdata(:,6)==-1;	%B choice
chosentarget(ind) = mod(chosentarget(ind)+4,8);
chosentarget(~chosentarget) = 8;
%
fr = hitdata(:,end); %firing rate

condTable = [offertype chosenjuice orientation tgtA chosentarget fr];