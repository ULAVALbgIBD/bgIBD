function [p, tp] = generate_transition_single_plus(pt, r, ref_prio, ref_tran)

%modify when pt = 0

ibd0 = 1;
ibd1 = 2;
ibd1_b = 3;

tp(ibd0:ibd1_b, ibd0:ibd1_b) = 0;

p(ibd0:ibd1_b) = 0;



[sum_prio, sum_tran] = path_tran(pt, r); 

p(ibd1) = sum_prio;
if( p(ibd1_b) >= 1 )
    %impossible in most cases
    disp('error in transmission probability, wrong inheritance');
end
if( p(ibd1) == 0 )
    tp(ibd1,ibd1) = 0;
else
    tp(ibd1,ibd1) = sum_tran/sum_prio;
end
tp(ibd1,ibd0) = 1 - tp(ibd1,ibd1);
tp(ibd1,ibd1_b) = 0;



p(ibd1_b) = ref_prio * (1-p(ibd1));
if( p(ibd1_b) == 0 )
    tp(ibd1_b,ibd1_b) = 0;
else
    tp(ibd1_b,ibd1_b) = ((1+exp(-2*r))/2)^ref_tran;
end
tp(ibd1_b,ibd0) = 1 - tp(ibd1_b,ibd1_b);
tp(ibd1_b,ibd1) = 0;


p(ibd0) = 1 - p(ibd1) - p(ibd1_b);

tp(ibd0,ibd1) = ( tp(ibd1,ibd0)  * p(ibd1) ) / ( p(ibd0) );
tp(ibd0,ibd1_b) = ( p(ibd1_b) *  tp(ibd1_b,ibd0) ) / p(ibd0);
tp(ibd0,ibd0) = 1 - tp(ibd0,ibd1) - tp(ibd0,ibd1_b);

if( tp(ibd0,ibd0) < 0 )
    disp('error in transmission probability, wrong inheritance');
end


end



