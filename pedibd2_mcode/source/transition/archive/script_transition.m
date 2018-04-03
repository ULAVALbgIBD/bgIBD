

%%%%%%%%transition probability%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

% test of sibling
% consider recombination rate is 1cM/Mb
% average marker interval distance 5kb

% default value for siblings, two paths of length 2 on paternal (maternal)
% side
pt1 = [0,2];
pt2 = [0,2];

pt1 = sum(input_allpaths(2:3,:));
pt2 = sum(input_allpaths(4:5,:));

% input_tpl{i}, is the transition from locus i-1 to locus i

%%%%%%%%%%%%background ibd%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bg_ibd(s_id) = 0.0113; %estimated from the whole population
% % bg_tran(s_id) = 0.018;
bg_tran(s_id) = 200;    %generations apart

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% transition probability is not correct for parent-child relationship
[input_pr, input_tpl] = generate_transition_all_loci(input_allpaths, sampled_markerlist, bg_ibd(s_id), bg_tran(s_id));

clear pt1 pt2 de_r de_pr de_tp temp;

%%
