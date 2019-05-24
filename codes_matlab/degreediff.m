function degdiff = degreediff(a1,a2,mode)
    if nargin == 2
        mode = 0;
    end
    if mode == 0    % [-180 180]
        degdiff = mod(mod(a2-a1,360)-180,360)-180;
    elseif mode == 1    % absolute [0 360]
        degdiff = min(abs(a1-a2),360-abs(a1-a2));
    elseif mode == 2    % [-90 270]
        degdiff = mod(mod(a2-a1,360)-270,360)-90;
    end
end