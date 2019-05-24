function [sessions] = sessionlistread(filename)

%
%
%

if 0
	filename = 'sessionlist_DT_PFC.m';
end

[fid,msg] = fopen(filename,'r');

if fid==-1 
   disp(sprintf('** error (sessionlistread) opening file %s; \n message: %s ',filename,msg))
   return
else
   
   i=0;
   while ~feof(fid)
      i = i+1;
      line = fgetl(fid);
      lineday(i,:) = line(1:7);
      linesess(i,:) = line(9);
	  lineflag(i,:) = str2num(line(11));
   end %while
   fclose(fid);
end

toan=find(lineflag(:,1)==1);

s = size(toan,1);
sessions = cell(s,1);

for i = 1:s
   sessions{i} = [lineday(toan(i),:),linesess(toan(i))];
end

