function [p, tp] = generate_transition_single(pt, r)

%modify when pt = 0

ibd0 = 1;
ibd1 = 2;

tp(ibd0:ibd1, ibd0:ibd1) = 0;

p(ibd0:ibd1) = 0;

[sum_prio, sum_tran] = path_tran(pt, r); 

p(ibd1) = sum_prio;
p(ibd0) = 1 - p(ibd1);
tp(ibd1,ibd1) = sum_tran/sum_prio;


tp(ibd1,ibd0) = 1 - tp(ibd1,ibd1);
tp(ibd0,ibd1) = (p(ibd1) - tp(ibd1,ibd1)*p(ibd1))/p(ibd0);
tp(ibd0,ibd0) = 1 - tp(ibd0,ibd1);


end