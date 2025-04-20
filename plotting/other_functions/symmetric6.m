function symtype = symmetric6(angle)
% include fully asymmetric + only one axial symmetric + only one central symmetric case

ang_sort = sort(angle); % sort angles in ascending order
ax_sym = [~any((ang_sort(1) + ang_sort(2))/2 - 90);  % when the first and the second angle is axial symmetric
    ~any((ang_sort(2) + ang_sort(3))/2 - 180);
    ~any((ang_sort(3) + ang_sort(4))/2 - 270);
    ~any(ang_sort(1) + ang_sort(4)-360)];
ct_sym = [~any(ang_sort(1) + 180 - ang_sort(3)); 
    ~any(ang_sort(2) + 180 - ang_sort(4))];

if (sum(ax_sym)==1 & sum(ct_sym)==0) 
    symtype = 1; % one axial symmetric
elseif (sum(ax_sym)==0 & sum(ct_sym)==1)
    symtype = 2; % one central symmetric
elseif (sum(ax_sym)==0 & sum(ct_sym)==0) 
    symtype = 3; % fully asymmetric
else
    symtype = 0;
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