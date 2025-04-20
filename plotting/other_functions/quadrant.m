function quad = quadrant(angle)
% output which quadrant the angle belong to
% Example:
% angle = [30 60 125]; % must be one-dimensional vector
% quad = quadrant(data)

for ii = 1:length(angle)
    temp = angle(ii);
switch true
    case temp >=0 & temp < 90
        quad(ii) = 1;
    case temp >=90 & temp < 180
        quad(ii) = 2;
    case temp >=180 & temp < 270
        quad(ii) = 3;
    case temp >=270 & temp < 360
        quad(ii) = 4;
end      
end