
function [input_pr_sc input_tpl_sc] = transition_sc(input_allpaths, sampled_markerlist)

%%%%%%%%transition probability%%%%%%%%%%%%%%%%%%%%%%%%%%


% test of sibling
% consider recombination rate is 1cM/Mb
% average marker interval distance 5kb

% default value for siblings, two paths of length 2 on paternal (maternal)
% side

% input_tpl{i}, is the transition from locus i-1 to locus i

%%%%%%%%%%%%background ibd%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bg_ibd = 0.0113; %estimated from the whole population
% % bg_tran(s_id) = 0.018;
bg_tran = 200;    %generations apart

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if( sum(input_allpaths) == 0 )
%     input_allpaths = [0,0,0,0,4];
else
    input_allpaths;
end
%inbreeding is fixed value, need to modify

[input_pr_sc input_tpl_sc] = generate_transition_all_loci_sc(input_allpaths, sampled_markerlist, bg_ibd, bg_tran);



end
