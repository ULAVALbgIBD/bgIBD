% double chain

function ep = dc_generate_emission_genotype(freq, ms, te, emission_option)

allele_map = emission_option.allele2genotype;
forward_map = emission_option.genotype2pair;
missing_code = emission_option.pair_missing_code;

ibd0 = 1;
ibd1 = 2;
ibd2 = 3;

g11 = allele_map(1,1);
g12 = allele_map(1,2);
g22 = allele_map(2,2);
g0 = missing_code;

nonmissing = (1-ms)^2;

e(1:3,1:10) = 0;
ep = e;

p1 = freq(1);
p2 = freq(2);

e(ibd0, forward_map(g11,g11)) = p1^2 * p1^2;
e(ibd0, forward_map(g11,g12)) = p1^2 * (2 * p1 * p2);
e(ibd0, forward_map(g11,g22)) = p1^2 * p2^2;
e(ibd0, forward_map(g12,g11)) = (2 * p1 * p2) * p1^2;
e(ibd0, forward_map(g12,g12)) = (2 * p1 * p2) * (2 * p1 * p2);
e(ibd0, forward_map(g12,g22)) = (2 * p1 * p2) * p2^2;
e(ibd0, forward_map(g22,g11)) = p2^2 * p1^2;
e(ibd0, forward_map(g22,g12)) = p2^2 * (2 * p1 * p2);
e(ibd0, forward_map(g22,g22)) = p2^2 * p2^2;

e(ibd1, forward_map(g11,g11)) = p1 * p1^2;
e(ibd1, forward_map(g11,g12)) = p1 * (p1 * p2);
e(ibd1, forward_map(g11,g22)) = 0;
e(ibd1, forward_map(g12,g11)) = p1 * (p2 * p1);
e(ibd1, forward_map(g12,g12)) = p1 * p2^2 + p1^2 * p2;
e(ibd1, forward_map(g12,g22)) = (p1 * p2) * p2;
e(ibd1, forward_map(g22,g11)) = 0;
e(ibd1, forward_map(g22,g12)) = (p2 * p1) * p2;
e(ibd1, forward_map(g22,g22)) = p2 * p2^2;

e(ibd2, forward_map(g11,g11)) = p1 * p1;
e(ibd2, forward_map(g11,g12)) = 0;
e(ibd2, forward_map(g11,g22)) = 0;
e(ibd2, forward_map(g12,g11)) = 0;
e(ibd2, forward_map(g12,g12)) = (2 * p1 * p2);
e(ibd2, forward_map(g12,g22)) = 0;
e(ibd2, forward_map(g22,g11)) = 0;
e(ibd2, forward_map(g22,g12)) = 0;
e(ibd2, forward_map(g22,g22)) = p2 * p2;


% typing error emission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p1 = 0.5;
p2 = 0.5;

er(1:9) = 0;

er(forward_map(g11,g11)) = p1^2 * p1^2;
er(forward_map(g11,g12)) = p1^2 * (2 * p1 * p2);
er(forward_map(g11,g22)) = p1^2 * p2^2;
er(forward_map(g12,g11)) = (2 * p1 * p2) * p1^2;
er(forward_map(g12,g12)) = (2 * p1 * p2) * (2 * p1 * p2);
er(forward_map(g12,g22)) = (2 * p1 * p2) * p2^2;
er(forward_map(g22,g11)) = p2^2 * p1^2;
er(forward_map(g22,g12)) = p2^2 * (2 * p1 * p2);
er(forward_map(g22,g22)) = p2^2 * p2^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = [g11,g12,g22]
    for j = [g11,g12,g22]
        ep(ibd0, forward_map(i,j)) = (1-te) * e(ibd0, forward_map(i,j)) + te * er(forward_map(i,j));
        ep(ibd1, forward_map(i,j)) = (1-te) * e(ibd1, forward_map(i,j)) + te * er(forward_map(i,j));
        ep(ibd2, forward_map(i,j)) = (1-te) * e(ibd2, forward_map(i,j)) + te * er(forward_map(i,j));
    end 
end

% add effect of missing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = [g11,g12,g22]
    for j = [g11,g12,g22]
        ep(ibd0, forward_map(i,j)) = nonmissing * ep(ibd0, forward_map(i,j));
        ep(ibd1, forward_map(i,j)) = nonmissing * ep(ibd1, forward_map(i,j));
        ep(ibd2, forward_map(i,j)) = nonmissing * ep(ibd2, forward_map(i,j));
    end 
end

ep(ibd0,g0) = 1 - nonmissing;
ep(ibd1,g0) = 1 - nonmissing;
ep(ibd2,g0) = 1 - nonmissing;

end


