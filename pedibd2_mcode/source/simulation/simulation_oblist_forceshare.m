

function [ref random_ibd_result random_ibs_result] = simulation_oblist_forceshare(p, s_id, input_pairlist, sm, pm)
% mandate that there is change in ibd status in the simulation

random_ibd_result(1:length(pm.sampled_markerlist(:,2))) = 0;

count = 0;

ref = find_breakpoints(random_ibd_result', pm.sampled_markerlist);

while( isempty(ref) )
    count = count + 1;
    [random_ibd_result random_ibs_result] = simulation_oblist(p, s_id, input_pairlist, sm, pm);
    ref = find_breakpoints(random_ibd_result', pm.sampled_markerlist);
end


ref = ref';
ref = [ref-1,ref];

% plot(random_ibd_result);

end

