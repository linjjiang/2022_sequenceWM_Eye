function ind = near_axes(angle,s_angle)
% check if the angles form a symmetric shape
% fully symmetric

ang_sort = sort(angle); % sort angles in ascending order
if any((abs(ang_sort) <= s_angle)) | any((abs(ang_sort-90) <= s_angle)) |...
        any((abs(ang_sort-180) <= s_angle)) | any((abs(ang_sort-360) <= s_angle))
    ind = 1;
else 
    ind = 0;
end