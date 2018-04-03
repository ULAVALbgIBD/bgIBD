function [pr, tp] = generate_transition_plus(all_paths, r, ref_prio, ref_tran)

ref_single = ref_prio;


% this relationship is correct even for parent-child transmission, on the
% allele level

[pr1p2p, tr1p2p] = generate_transition_single_plus(all_paths(2,:), r, ref_single, ref_tran);
[pr1p2m, tr1p2m] = generate_transition_single_plus(all_paths(3,:), r, ref_single, ref_tran);
[pr1m2p, tr1m2p] = generate_transition_single_plus(all_paths(4,:), r, ref_single, ref_tran);
[pr1m2m, tr1m2m] = generate_transition_single_plus(all_paths(5,:), r, ref_single, ref_tran);
% 
% 
% 
% 
% [1p2p, 1p2m, 1m2p, 1m2m] included states
states = [
    1, 1, 1, 1; %1 ibd0
    2, 1, 1, 1; %2 ibd1
    3, 1, 1, 1; %3
    1, 2, 1, 1; %4
    1, 3, 1, 1; %5
    1, 1, 2, 1; %6
    1, 1, 3, 1; %7
    1, 1, 1, 2; %8
    1, 1, 1, 3; %9
    2, 1, 1, 2; %10 ibd2
    3, 1, 1, 2; %11
    2, 1, 1, 3; %12
    3, 1, 1, 3; %13
    1, 2, 2, 1; %14
    1, 3, 2, 1; %15
    1, 2, 3, 1; %16
    1, 3, 3, 1; %17
    ];


[num, fields] = size(states);
if( fields ~= 4 )
    disp('error in ibd states');
    return;
end
tr = (zeros(num,num));
tp = (zeros(num,num));

prio = (zeros(1,num));
pr = (zeros(1,num));

for i = 1:length(states)
    prio(i) = pr1p2p(states(i,1)) * pr1p2m(states(i,2)) * pr1m2p(states(i,3)) * pr1m2m(states(i,4));
    if( pr1p2p(2) + pr1p2m(2) == 1 )
        if( states(i,1) ~= 2 && states(i,2) ~= 2 )
            prio(i) = 0;
        end
    end
    if( pr1p2p(2) + pr1m2p(2) == 1 )
        if( states(i,1) ~= 2 && states(i,3) ~= 2 )
            prio(i) = 0;
        end
    end
    if( pr1p2m(2) + pr1m2m(2) == 1 )
        if( states(i,2) ~= 2 && states(i,4) ~= 2 )
            prio(i) = 0;
        end
    end
    if( pr1m2p(2) + pr1m2m(2) == 1 )
        if( states(i,3) ~= 2 && states(i,4) ~= 2 )
            prio(i) = 0;
        end
    end
end
if( sum(prio) > 0 )
    pr = prio./sum(prio);
else
    pr(1:num) = 0;
end

for i = 1:num
    for j = 1:num
        tr(i,j) = 1 * tr1p2p(states(i,1), states(j,1)) * tr1p2m(states(i,2), states(j,2)) * tr1m2p(states(i,3), states(j,3)) * tr1m2m(states(i,4), states(j,4)); 
        if( pr(i) == 0 || pr(j) == 0 )
            tr(i,j) = 0;
        end
    end
    if( sum(tr(i,:)) > 0 )
        tp(i,:) = tr(i,:)./sum(tr(i,:));
    else
        tp(i,:) = 0;
    end
end










end













