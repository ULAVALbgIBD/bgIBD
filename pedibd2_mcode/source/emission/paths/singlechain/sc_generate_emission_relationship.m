function ep = sc_generate_emission_relationship(a1, a2, direct_source, genotype, ms, te)

% direct source, left side is paternal, minus means paternal
% i: -i, i

% divide into three areas, exclusively belongs to a1, a2, and shared by a1
% and a2

%single chain

ibd0 = 1;
ibd1 = 2;

g11 = 1;
g12 = 2;
g21 = 3;
g22 = 4;
g00 = 5;

nonmissing = (1-ms);

e(1:2,1:4) = 0;
ep = e;

p1 = freq(1);
p2 = freq(2);

e(ibd0,g11) = p1^2;
e(ibd0,g12) = 2 * p1 * p2;
e(ibd0,g22) = p2^2;

e(ibd1,g11) = p1;
e(ibd1,g12) = 0;
e(ibd1,g22) = p2;

% typing error emission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p1 = 0.5;
p2 = 0.5;

e11 = p1^2;
e12 = 2*p1*p2;
e22 = p2^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ep(ibd0,g11) = (1-te)*e(ibd0,g11) + te*e11;
ep(ibd0,g12) = (1-te)*e(ibd0,g12) + te*e12;
ep(ibd0,g22) = (1-te)*e(ibd0,g22) + te*e22;

ep(ibd1,g11) = (1-te)*e(ibd1,g11) + te*e11;
ep(ibd1,g12) = (1-te)*e(ibd1,g12) + te*e12;
ep(ibd1,g22) = (1-te)*e(ibd1,g22) + te*e22;

% add effect of missing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ep(ibd0,g11) = nonmissing * ep(ibd0,g11);
ep(ibd0,g12) = nonmissing * ep(ibd0,g12);
ep(ibd0,g22) = nonmissing * ep(ibd0,g22);
ep(ibd0,g00) = 1 - nonmissing;

ep(ibd1,g11) = nonmissing * ep(ibd1,g11);
ep(ibd1,g12) = nonmissing * ep(ibd1,g12);
ep(ibd1,g22) = nonmissing * ep(ibd1,g22);
ep(ibd1,g00) = 1 - nonmissing;

end


