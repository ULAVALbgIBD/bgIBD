function [ep error] = sibling_generate_emission(f_geno, m_geno, freq, ms, te, emission_option)



ep = [];
error = 0;


father_genotyped = (f_geno(1)~=0) && (f_geno(2)~=0);
mother_genotyped = (m_geno(1)~=0) && (m_geno(2)~=0);


ep_merge_genotype = sibling_generate_emission_genotype(f_geno, m_geno, freq, 0, 0, emission_option);
ep_merge_oblist = emission_genotype2ob(ep_merge_genotype, emission_option);


[nrow, ncol] = size(ep_merge_oblist);

ibd0 = 1;
ibdL = 2;
ibdR = 3;
ibdLR = 4;
ibs0 = 1;
ibs1 = 2;
ibs2 = 3;
none = 4;

te = te * 4;

if( nrow == 4 && ncol == 4 )
    minor_af = min(freq(1:2));
    te2 = minor_af / 10;
    ep_merge_oblist(ibdLR, ibs0) = te2 / 10;
    ep_merge_oblist(ibdLR, ibs1) = te2;
    ep_merge_oblist(ibdLR, ibs2) = ep_merge_oblist(ibdLR, ibs2) * (1 - ms - te2 - te2/10);
    ep_merge_oblist(ibdLR, none) = ms;
    
    te1 = minor_af * minor_af / 40;
    te3 = 0.05;
    if( mother_genotyped )
        ep_merge_oblist(ibdL, ibs0) = te1;
        ep_merge_oblist(ibdL, ibs1) = ep_merge_oblist(ibdL, ibs1) * (1 - ms - te1 - te3) + 0.5 * te3;
        ep_merge_oblist(ibdL, ibs2) = ep_merge_oblist(ibdL, ibs2) * (1 - ms - te1 - te3) + 0.5 * te3;
        ep_merge_oblist(ibdL, none) = ms;
    else
        ep_merge_oblist(ibdL, ibs0) = te1;
        ep_merge_oblist(ibdL, ibs1) = ep_merge_oblist(ibdL, ibs1) * (1 - ms - te1);
        ep_merge_oblist(ibdL, ibs2) = ep_merge_oblist(ibdL, ibs2) * (1 - ms - te1);
        ep_merge_oblist(ibdL, none) = ms;
    end
    
    if( father_genotyped )
        ep_merge_oblist(ibdR, ibs0) = te1;
        ep_merge_oblist(ibdR, ibs1) = ep_merge_oblist(ibdR, ibs1) * (1 - ms - te1 - te3) +  0.5 * te3;
        ep_merge_oblist(ibdR, ibs2) = ep_merge_oblist(ibdR, ibs2) * (1 - ms - te1 - te3) +  0.5 * te3;
        ep_merge_oblist(ibdR, none) = ms;
    else
        ep_merge_oblist(ibdR, ibs0) = te1;
        ep_merge_oblist(ibdR, ibs1) = ep_merge_oblist(ibdR, ibs1) * (1 - ms - te1);
        ep_merge_oblist(ibdR, ibs2) = ep_merge_oblist(ibdR, ibs2) * (1 - ms - te1);
        ep_merge_oblist(ibdR, none) = ms;
    end
    
    if( father_genotyped || mother_genotyped )
        ep_merge_oblist(ibd0, ibs0) = ep_merge_oblist(ibd0, ibs0) * (1 - ms - te3) + 0.125 * te3;
        ep_merge_oblist(ibd0, ibs1) = ep_merge_oblist(ibd0, ibs1) * (1 - ms - te3) + 0.500 * te3;
        ep_merge_oblist(ibd0, ibs2) = ep_merge_oblist(ibd0, ibs2) * (1 - ms - te3) + 0.375 * te3;
        ep_merge_oblist(ibd0, none) = ms;        
    else
        ep_merge_oblist(ibd0, ibs0) = ep_merge_oblist(ibd0, ibs0) * (1 - ms - te) + 0.125 * te;
        ep_merge_oblist(ibd0, ibs1) = ep_merge_oblist(ibd0, ibs1) * (1 - ms - te) + 0.500 * te;
        ep_merge_oblist(ibd0, ibs2) = ep_merge_oblist(ibd0, ibs2) * (1 - ms - te) + 0.375 * te;
        ep_merge_oblist(ibd0, none) = ms;
    end
    
else
    error = 1;
    disp('emission states not supported');
    return;
end

% [1p2p, 1p2m, 1m2p, 1m2m] included states
states = [
    1, 1, 1, 1; %ibd0:1
    2, 1, 1, 1; %ibd_l:2
    3, 1, 1, 1; %ibd0:3
    1, 2, 1, 1; %ina:4 
    1, 3, 1, 1; %ibd0:5
    1, 1, 2, 1; %ina:6
    1, 1, 3, 1; %ibd0:7
    1, 1, 1, 2; %ibd_r:8
    1, 1, 1, 3; %ibd0:9
    2, 1, 1, 2; %ibd_lr:10
    3, 1, 1, 2; %ibd_r:11
    2, 1, 1, 3; %ibd_l:12
    3, 1, 1, 3; %ibd0:13
    1, 2, 2, 1; %ina:14
    1, 3, 2, 1; %ina:15
    1, 2, 3, 1; %ina:16
    1, 3, 3, 1; %ibd0:17
    ];

ibd0 = 1;
ibd_l = 2;
ibd_r = 3;
ibd_lr = 4;

% inaccessible states in prio probability of siblings, non-inbreeding
% situation
% emission is "don't care"
mapping([4,5,6,7,14,15,16,17]) = ibd0;

if( father_genotyped )
    if( mother_genotyped )
        % in case both parents genotyped
        % background ibd is explicit, absorbed in ibd0
        mapping([1,3,9,13]) = ibd0;
        mapping([2,12]) = ibd_l;
        mapping([8,11]) = ibd_r;
        mapping([10]) = ibd_lr;
        % inbreeding suppressed
        mapping([4,5,6,7]) = ibd0; % parents sharing 1 ibd
        mapping([14,15,16,17]) = ibd0; % parents sharing 2 ibd
    else
        % if only father is genotyped
        % paternal side background ibd is explict
        % only consider maternal side background ibd
        mapping([1]) = ibd0;
        mapping([2]) = ibd_l;
        mapping([3]) = ibd0;
        mapping([8]) = ibd_r;
        mapping([9]) = ibd_r; %bg_r
        mapping([10]) = ibd_lr; 
        mapping([11]) = ibd_r; %bg_l
        mapping([12]) = ibd_lr; 
        mapping([13]) = ibd_r; %bg_lr
        % inbreeding suppressed
        mapping([4,5,6,7]) = ibd0; % parents sharing 1 ibd
        mapping([14,15,16,17]) = ibd0; % parents sharing 2 ibd
    end
else
    if( mother_genotyped )
        %  absord maternal side background ibd, as ibd0
        mapping([1]) = ibd0;
        mapping([2]) = ibd_l;
        mapping([3]) = ibd_l;
        mapping([8]) = ibd_r;
        mapping([9]) = ibd0;
        mapping([10]) = ibd_lr;
        mapping([11]) = ibd_lr;
        mapping([12]) = ibd_l;
        mapping([13]) = ibd_l;
        % inbreeding suppressed
        mapping([4,5,6,7]) = ibd0; % parents sharing 1 ibd
        mapping([14,15,16,17]) = ibd0; % parents sharing 2 ibd       
    else
        mapping([1]) = ibd0;
        mapping([2:9]) = ibd_l; 
        % assign to ibd_r is also OK, same value if use uniform allele
        % frequency for both parents
        mapping([10:17]) = ibd_lr;
    end    
end



ep_oblist(1:length(mapping),1:ncol) = 0;

for i = 1:length(mapping)
    for j = 1:ncol
        ep_oblist(i,j) = ep_merge_oblist(mapping(i), j);
    end
end

ep = ep_oblist;


end


