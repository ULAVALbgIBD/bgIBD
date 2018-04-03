% double chain

function ep = sibling_generate_emission_ibs_side_unpaired(f_geno, m_geno, freq, ms, te)

% non-paired distant sharing p-m
% assuming
% sharing no parental alleles
% sharing is through sharing between parents


ibd1 = 1;
ibd2 = 2;


ibs0 = 1;
ibs1 = 2;
ibs2 = 3;
ibs3 = 4; %missing

nonmissing = (1-ms)^2;

te2 = 1-(1-te)^2;
te4 = 1-(1-te)^4;
te3 = 1-(1-te)^3;


e = sibling_side_unpaired(f_geno, m_geno, freq);



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
    ep(ibd1, i) = (1-te4) * e(ibd1, i) + te4 * er(i);
    ep(ibd2, i) = (1-te4) * e(ibd2, i) + te4 * er(i);
end

% add effect of missing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = ibs0:ibs2
        ep(ibd1, i) = nonmissing * ep(ibd1, i);
        ep(ibd2, i) = nonmissing * ep(ibd2, i);
end

ep(ibd1,ibs3) = 1 - nonmissing;
ep(ibd2,ibs3) = 1 - nonmissing;

end





