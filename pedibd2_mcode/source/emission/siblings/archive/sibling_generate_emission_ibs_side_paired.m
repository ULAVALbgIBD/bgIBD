% double chain

function ep = sibling_generate_emission_ibs_side_paired(f_geno, m_geno, freq, ms, te)

% paired distant sharing p-p
% assuming
% sharing no parental alleles
% sharing is through sharing within parents


ibd_l = 1;
ibd_r = 2;
ibd_lr = 3;


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
    ep(ibd_l, i) = (1-te4) * e(ibd_l, i) + te4 * er(i);
    ep(ibd_r, i) = (1-te4) * e(ibd_r, i) + te4 * er(i);
    ep(ibd_lr, i) = (1-te4) * e(ibd_lr, i) + te4 * er(i);
end

% add effect of missing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = ibs0:ibs2
        ep(ibd_l, i) = nonmissing * ep(ibd_l, i);
        ep(ibd_r, i) = nonmissing * ep(ibd_r, i);
end

ep(ibd_l,ibs3) = 1 - nonmissing;
ep(ibd_r,ibs3) = 1 - nonmissing;
ep(ibd_lr,ibs3) = 1 - nonmissing;

end