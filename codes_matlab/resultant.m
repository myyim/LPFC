function [r, ang] = resultant(vect,m)

if nargin < 2 || isempty(m)
    m = 1;
end

s = 0;
for j = 0:7
    s = s + vect(j+1)*exp(1i*pi*m*j/4);
end

s = s/sum(vect);

r = abs(s);

if angle(s) < 0
    ang = (2*pi + angle(s))/m*180/pi;
else
    ang = angle(s)/m*180/pi;
end