function [pr, tpl] = generate_transition_all_loci_sc(all_paths, markerlist, de_ibd, de_tran)

tpl = [];
for i = 1:length(markerlist)
    tpl{i} = [];
    if( i == 1 )
        r = 0.01*(markerlist(i,2))/1000000;
    else
        r = 0.01*(markerlist(i,2)-markerlist(i-1,2))/1000000;
    end
    [pr, tpl{i}] = generate_transition_single_plus(all_paths, r, de_ibd, de_tran);
    
end

end