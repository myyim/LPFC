function [monkey, brainarea] = arearead_DT(session, cellname)
% [monkey, brainarea] = arearead_DT(session, cellname) returns the monkey
% name and the brain area of the given cellname in the given session.
% e.g. arearead_DT('B110720a','1')

filename = ['DT_',session]; eval(filename)
monkey = parsession.monkey;
electrode = str2double(cellname(1));
locs = parsession.locations;
ind = find(locs(:,1)==electrode);
loc = locs(ind,:);

% if loc(4) == 8;
% 	brainarea = 'dlPFC external';
% end

%brainarea = 'dPS';

% if abs(loc(3)) > 14;
% 	brainarea = 'ITC';
% else
% 	brainarea = 'dPFC';
% end
% 
% x = [32 33 34 35 36 37];
% y = [14 14 14 14 14 14];
% 
% if loc(4) == 10;
% 	%m = find(loc(2) == x);
% 	if abs(loc(3)) > 14;
% 		brainarea = 'ITC';
% 	else
% 		brainarea = 'vPS';
% 	end
% 	
% end

if loc(4) == 0;
	brainarea = 'OFC';
end

if loc(4) == 9;
	brainarea = 'dlPFC';
end

if loc(4) == 9.5;
	brainarea = 'dlPFC';
end

if loc(4) == 10;
	brainarea = 'vlPFC';
end

if loc(4) == 14;
	brainarea = 'LOFC';
end

if loc(4) == 11;
	brainarea = 'vmPFC';
end

if loc(4) == 12;
	brainarea = 'BG';
end

% if isequal(brainarea,'ITC') ||  isequal(brainarea,'vPS'),
% 	brainarea = 'vPFC';
% end



