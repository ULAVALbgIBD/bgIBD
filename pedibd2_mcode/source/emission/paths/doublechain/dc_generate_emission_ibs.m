% this function could be replaced by calling
% dc_generate_emission_genotype
% dc_generate_emission_ibs_genotype

function ep = dc_generate_emission_ibs(freq, ms, te)

ibd0 = 1;
ibd1 = 2;
ibd2 = 3;

ibs0 = 1;
ibs1 = 2;
ibs2 = 3;
missing = 4;

nonmissing = (1-ms)^2;

e(1:3,1:3) = 0;
ep = e;

p1 = freq(1);
p2 = freq(2);

e(ibd0,ibs0) = 2 * p1^2 * p2^2;
e(ibd0,ibs1) = 4 * p1^3 * p2 + 4 * p1 * p2^3;
e(ibd0,ibs2) = p1^4 + (2*p1*p2)^2 + p2^4;

e(ibd1,ibs0) = 0;
e(ibd1,ibs1) = 2 * p1 * p2;
e(ibd1,ibs2) = p1 * p1 + p2 * p2;

e(ibd2,ibs0) = 0;
e(ibd2,ibs1) = 0;
e(ibd2,ibs2) = 1;

% typing error emission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p1 = 0.5;
p2 = 0.5;

e0 = 2 * p1^2 * p2^2;
e1 = 4 * p1^3 * p2 + 4 * p1 * p2^3;
e2 = p1^4 + (2*p1*p2)^2 + p2^4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ep(ibd0,ibs0) = (1-te)*e(ibd0,ibs0) + te*e0;
ep(ibd0,ibs1) = (1-te)*e(ibd0,ibs1) + te*e1;
ep(ibd0,ibs2) = (1-te)*e(ibd0,ibs2) + te*e2;

ep(ibd1,ibs0) = (1-te)*e(ibd1,ibs0) + te*e0;
ep(ibd1,ibs1) = (1-te)*e(ibd1,ibs1) + te*e1;
ep(ibd1,ibs2) = (1-te)*e(ibd1,ibs2) + te*e2;

ep(ibd2,ibs0) = (1-te)*e(ibd2,ibs0) + te*e0;
ep(ibd2,ibs1) = (1-te)*e(ibd2,ibs1) + te*e1;
ep(ibd2,ibs2) = (1-te)*e(ibd2,ibs2) + te*e2;

% add effect of missing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ep(ibd0,ibs0) = nonmissing * ep(ibd0,ibs0);
ep(ibd0,ibs1) = nonmissing * ep(ibd0,ibs1);
ep(ibd0,ibs2) = nonmissing * ep(ibd0,ibs2);
ep(ibd0,missing) = 1 - nonmissing;

ep(ibd1,ibs0) = nonmissing * ep(ibd1,ibs0);
ep(ibd1,ibs1) = nonmissing * ep(ibd1,ibs1);
ep(ibd1,ibs2) = nonmissing * ep(ibd1,ibs2);
ep(ibd1,missing) = 1 - nonmissing;

ep(ibd2,ibs0) = nonmissing * ep(ibd2,ibs0);
ep(ibd2,ibs1) = nonmissing * ep(ibd2,ibs1);
ep(ibd2,ibs2) = nonmissing * ep(ibd2,ibs2);
ep(ibd2,missing) = 1 - nonmissing;

end


