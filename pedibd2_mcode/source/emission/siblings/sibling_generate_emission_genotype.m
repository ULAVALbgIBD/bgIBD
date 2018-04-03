% double chain

function ep = sibling_generate_emission_genotype(f_geno, m_geno, freq, ms, te, emission_option)

% here only consider true ibd given both parents are genotyped
% bypass all other kinship, bypass all other situations

allele_map = emission_option.allele2genotype;
forward_map = emission_option.genotype2pair;
missing_code = emission_option.pair_missing_code;

g11 = allele_map(1,1);
g12 = allele_map(1,2);
g22 = allele_map(2,2);
g0 = missing_code;

ibd0 = 1;
ibd_l = 2;
ibd_r = 3;
ibd_lr = 4;


nonmissing = (1-ms)^2;


e = sibling_core_genotype(f_geno, m_geno, freq, emission_option);



if( f_geno(1) ~= 0 && f_geno(2) ~= 0 )
    father_genotyped = 1;
else
    father_genotyped = 0;
end

if( m_geno(1) ~= 0 && m_geno(2) ~= 0 )
    mother_genotyped = 1;
else
    mother_genotyped = 0;
end



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

te = 0;
for i = [g11,g12,g22]
    for j = [g11,g12,g22]
        k = forward_map(i,j);
        ep(ibd0, k) = (1 - te) * e(ibd0, k) + ( te ) * er(k);
        ep(ibd_l, k) = (1 - te) * e(ibd_l, k) + ( te ) * er(k);
        ep(ibd_r, k) = (1 - te) * e(ibd_r, k) + ( te ) * er(k); 
        ep(ibd_lr, k) = (1 - te) * e(ibd_lr, k) + ( te ) * er(k);
    end
end

% add effect of missing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = [g11,g12,g22]
    for j = [g11,g12,g22]
        k = forward_map(i,j);
        ep(ibd0, k) = nonmissing * ep(ibd0, k);
        ep(ibd_l, k) = nonmissing * ep(ibd_l, k);
        ep(ibd_r, k) = nonmissing * ep(ibd_r, k);
        ep(ibd_lr, k) = nonmissing * ep(ibd_lr, k);
    end
end

ep(ibd0,g0) = 1 - nonmissing;
ep(ibd_l,g0) = 1 - nonmissing;
ep(ibd_r,g0) = 1 - nonmissing;
ep(ibd_lr,g0) = 1 - nonmissing;

end





