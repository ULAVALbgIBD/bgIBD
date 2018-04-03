function [pr, tp] = generate_transition(pt1, pt2, r)

ibd0 = 1;
ibd1 = 2;
ibd2 = 3;

[p1, t1] = generate_transition_single(pt1, r);
[p2, t2] = generate_transition_single(pt2, r);

tp(ibd2,ibd2) = t1(ibd1,ibd1)*t2(ibd1,ibd1);
tp(ibd2,ibd1) = t1(ibd1,ibd0)*t2(ibd1,ibd1) + t1(ibd1,ibd1)*t2(ibd1,ibd0);
tp(ibd2,ibd0) = t1(ibd1,ibd0)*t2(ibd1,ibd0);


temp1 = (p1(ibd1)*p2(ibd0))*(t1(ibd1,ibd1)*t2(ibd0,ibd1));
temp2 = (p1(ibd0)*p2(ibd1))*(t1(ibd0,ibd1)*t2(ibd1,ibd1));
tp(ibd1,ibd2) = ( temp1 + temp2 ) / (p1(ibd1)*p2(ibd0) + p1(ibd0)*p2(ibd1));

temp1 = (p1(ibd1)*p2(ibd0))*(t1(ibd1,ibd1)*t2(ibd0,ibd0));
temp2 = (p1(ibd0)*p2(ibd1))*(t1(ibd0,ibd0)*t2(ibd1,ibd1));
temp3 = (p1(ibd1)*p2(ibd0))*(t1(ibd1,ibd0)*t2(ibd0,ibd1)); %could be ignored
temp4 = (p1(ibd0)*p2(ibd1))*(t1(ibd0,ibd1)*t2(ibd1,ibd0)); %could be ignored
tp(ibd1,ibd1) = ( temp1 + temp2 + temp3 + temp4 ) / (p1(ibd1)*p2(ibd0) + p1(ibd0)*p2(ibd1));
% decompose a conditional probability, should it be mutually exclusive?

temp1 = (p1(ibd1)*p2(ibd0))*(t1(ibd1,ibd0)*t2(ibd0,ibd0));
temp2 = (p1(ibd0)*p2(ibd1))*(t1(ibd0,ibd0)*t2(ibd1,ibd0));
tp(ibd1,ibd0) = ( temp1 + temp2 ) / (p1(ibd1)*p2(ibd0) + p1(ibd0)*p2(ibd1));


tp(ibd0,ibd0) = t1(ibd0,ibd0)*t2(ibd0*ibd0);
tp(ibd0,ibd1) = t1(ibd0,ibd1)*t2(ibd0*ibd0) + t1(ibd0,ibd0)*t2(ibd0,ibd1);
tp(ibd0,ibd2) = t1(ibd0,ibd1)*t2(ibd0,ibd1);


pr(ibd0) = p1(ibd0)*p2(ibd0);
pr(ibd1) = p1(ibd0)*p2(ibd1) + p1(ibd1)*p2(ibd0);
pr(ibd2) = p1(ibd1)*p2(ibd1);


end