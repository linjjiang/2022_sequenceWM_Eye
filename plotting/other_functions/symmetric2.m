function sym = symmetric2(angle)
% check if the angles form a symmetric shape

ang_sort = sort(angle); % sort angles in ascending order
if (ang_sort(1) + ang_sort(2))/2 == 90 | ...
        (ang_sort(2) + ang_sort(3))/2 == 180 | ...
        (ang_sort(3) + ang_sort(4))/2 == 270 | ...
        (ang_sort(1) + ang_sort(4)-360)/2 == 0
    sym = 1;
else 
    sym = 0;
end
% x = 5.*cosd(ang_sort); y = 5.*sind(ang_sort); % calculate x & y of four corners
% h = polyshape(x,y); % create this polygon
% [cx, cy] = centroid(h); % detect the centroid of this polygon
% 
% if ~cx | ~cy % if cx == 0 or cy == 0
%     sym = 1; % symmetric
% else
%     sym = 0; % asymmetric
% end
% 
% end