function [merged_assignment error] = merge_segments_secure(assignment, family)


merged_assignment = [];
error = 0;
global debug_mode;

[~, ~, error] = check_allsegments(assignment, family);
if( error ~= 0 )
    disp('error in global IBD');
    return;
end

time = cputime;

[merged_assignment error] = merge_segments(assignment);

if( debug_mode == 1 )
    display(['merging costs time: ', num2str(cputime - time), ' seconds']); 
end

end

