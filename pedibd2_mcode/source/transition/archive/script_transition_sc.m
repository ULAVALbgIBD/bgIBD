

%%%%%%%%transition probability%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

% test of sibling
% consider recombination rate is 1cM/Mb
% average marker interval distance 5kb

% default value for siblings, two paths of length 2 on paternal (maternal)
% side

% input_tpl{i}, is the transition from locus i-1 to locus i

%%%%%%%%%%%%background ibd%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bg_ibd(s_id_sc) = 0.0113; %estimated from the whole population
% % bg_tran(s_id) = 0.018;
bg_tran(s_id_sc) = 200;    %generations apart

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

input_allpaths = [0,0,0,0,4];
%inbreeding is fixed value, need to modify

[input_pr_sc, input_tpl_sc] = generate_transition_all_loci_sc(input_allpaths, sampled_markerlist, bg_ibd(s_id_sc), bg_tran(s_id_sc));

clear pt1 pt2 de_r de_pr de_tp temp;

%%
