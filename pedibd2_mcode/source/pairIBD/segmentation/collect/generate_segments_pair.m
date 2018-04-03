% expand to matrix relationship from pairwise relationship

function [ pair_all error ] = generate_segments_pair( vit, likelihood, likelihoodIBD1, map )

global debug_mode;

pair_all = [];
error = 0;

% 1: ibd0
% 2: ibd1
% 3: ibd2
% make sure this assignment is correct

time = cputime;

markers =  generate_consistent_interval_viterbi( vit );
if( isempty(markers) )
    error = 1;
    disp('error in segmentation');
    return;
end

[num_m c] = size(markers);
if( num_m <= 0 || c ~= 3 )
    error = 1;
    disp('error in segmentation');
    return;
end

nloc = markers(num_m, 2);
n_ind = length(map);

stat_relation_prob(1:n_ind,1:n_ind,1:3) = 0;
stat_relation(1:n_ind,1:n_ind) = 0;
prob_IBD1(1:n_ind,1:n_ind,1:2,1:2) = 0;

vit_relation(1:n_ind, 1:n_ind) = 0;

pair_vit = zeros(num_m, n_ind, n_ind);
pair_pos = zeros(num_m, n_ind, n_ind, 3);
pair_posIBD1 = zeros(num_m, n_ind, n_ind, 2, 2);
pair_max = zeros(num_m, n_ind, n_ind);

for l = 1:num_m

    m_l = markers(l,3);
    
    for i = 1:n_ind
        vit_relation(i,i) = 2;       
        for j = i+1:n_ind
     
            vit_relation(i,j) = vit(map(i,j), m_l,1)-1;
            vit_relation(j,i) = vit(map(j,i), m_l,1)-1;           

        end
    end
    
    pair_vit(l, 1:n_ind, 1:n_ind) = vit_relation;
    
end



for l = 1:num_m

    s_l = markers(l,1);
    t_l = markers(l,2);
    
    interval = (t_l - s_l + 1);
    if( interval < 1 )
        error = 1;
        disp('error in segmentation');
        return;
    end
    
    for i = 1:n_ind
        stat_relation_prob(i,i,1:2) = 0;
        stat_relation_prob(i,i,3) = 1;
        stat_relation(i,i) = 2;        
        prob_IBD1(i,i,1,1) = 1;
        prob_IBD1(i,i,1,2) = 0;
        prob_IBD1(i,i,2,1) = 0;
        prob_IBD1(i,i,2,2) = 1;
        
        for j = i+1:n_ind

            stat_relation_prob(i,j,1:3) = sum(likelihood(map(i,j), s_l:t_l,1:3), 2)./interval;
            stat_relation_prob(j,i,1:3) = stat_relation_prob(i,j,1:3);
            
            [temp ix] = max( stat_relation_prob(i,j,:) );
            stat_relation(i,j) = ix - 1;
            stat_relation(j,i) = stat_relation(i,j);
            prob_IBD1(i,j,1:2,1:2) = sum(likelihoodIBD1(map(i,j), s_l:t_l,1:2,1:2), 2)./interval;
            prob_IBD1(j,i,1,1) = prob_IBD1(i,j,1,1);
            prob_IBD1(j,i,1,2) = prob_IBD1(i,j,2,1);
            prob_IBD1(j,i,2,1) = prob_IBD1(i,j,1,2);
            prob_IBD1(j,i,2,2) = prob_IBD1(i,j,2,2);
        end
    end
    
    % sum is much faster than mean ?
    
    pair_pos(l, 1:n_ind, 1:n_ind, 1:3) = stat_relation_prob;
    pair_posIBD1(l, 1:n_ind, 1:n_ind, 1:2, 1:2) = prob_IBD1;
    pair_max(l, 1:n_ind, 1:n_ind) = stat_relation;
    
end

pair_all.pair_vit = pair_vit;
pair_all.pair_pos = pair_pos;
pair_all.pair_posIBD1 = pair_posIBD1;
pair_all.pair_max = pair_max;

pair_all.intervals = markers;

if( debug_mode == 1 )
    display(['          partition finalizing costs time: ', num2str(cputime - time), ' seconds']);  
end

disp(['total ', num2str(length(pair_all.intervals(:,1))), ' chromosomal regions']);

end













