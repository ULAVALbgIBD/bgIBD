% double chain

function ep = sibling_generate_emission_ibs(f_geno, m_geno, freq, ms, te, emission_option, allele2ibd, nf)

% here only consider true ibd given both parents are genotyped
% bypass all other kinship, bypass all other situations

ibd0 = 1;
ibd_l = 2;
ibd_r = 3;
ibd_lr = 4;


ibs0 = 1;
ibs1 = 2;
ibs2 = 3;
ibs3 = 4; %missing

nonmissing = (1-ms)^2;

perror(1) = te;
perror(2) = 1-(1-te)^2;
perror(4) = 1-(1-te)^4;
perror(3) = 1-(1-te)^3;



%e2 = sibling_core(f_geno, m_geno, freq);
e = sibling_core_genotype(f_geno, m_geno, freq, emission_option, allele2ibd, nf);
e = emission_genotype2ob(e, emission_option);


if( f_geno(1) == 0 || f_geno(2) == 0 )
    father_missing = 1;
else
    father_missing = 0;
end

if( m_geno(1) == 0 || m_geno(2) == 0 )
    mother_missing = 1;
else
    mother_missing = 0;
end



% typing error emission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p1 = 0.5;
p2 = 0.5;

er(1:3) = 0;

er(ibs0) = (2 * p1^2 * p2^2);
er(ibs1) = 2 * (p1^2 + p2^2) * (2 * p1 * p2);
er(ibs2) = p1^2 * p1^2 + p2^2 * p2^2 + (2 * p1 * p2) * (2 * p1 * p2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = ibs0:ibs2
    ep(ibd0, i) = (1-perror(4-father_missing-mother_missing)) * e(ibd0, i) + perror(4-father_missing-mother_missing) * er(i);
    ep(ibd_l, i) = (1-perror(3-mother_missing)) * e(ibd_l, i) + perror(3-mother_missing) * er(i);
    ep(ibd_r, i) = (1-perror(3-father_missing)) * e(ibd_r, i) + perror(3-father_missing) * er(i); 
    ep(ibd_lr, i) = (1-perror(2)) * e(ibd_lr, i) + perror(2) * er(i); 
end

% add effect of missing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = ibs0:ibs2
        ep(ibd0, i) = nonmissing * ep(ibd0, i);
        ep(ibd_l, i) = nonmissing * ep(ibd_l, i);
        ep(ibd_r, i) = nonmissing * ep(ibd_r, i);
        ep(ibd_lr, i) = nonmissing * ep(ibd_lr, i);
end

ep(ibd0,ibs3) = 1 - nonmissing;
ep(ibd_l,ibs3) = 1 - nonmissing;
ep(ibd_r,ibs3) = 1 - nonmissing;
ep(ibd_lr, ibs3) = 1 - nonmissing;

end





