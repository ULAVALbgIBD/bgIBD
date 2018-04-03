
function [ep error] = generate_emission(freq, ms, te, emission_option)


ep = [];
error = 0;

ep_merge_genotype = dc_generate_emission_genotype(freq, 0, 0, emission_option);
ep_merge_oblist = emission_genotype2ob(ep_merge_genotype, emission_option);

[nrow, ncol] = size(ep_merge_oblist);

ibd0 = 1;
ibd1 = 2;
ibd2 = 3;
ibs0 = 1;
ibs1 = 2;
ibs2 = 3;
none = 4;

if( nrow == 3 && ncol == 4 )
    
    temp1 = sum(ep_merge_oblist(ibd0, 2:3));
    temp2 = sum(ep_merge_oblist(ibd1, 2:3));
    temp3 = sum(ep_merge_oblist(1:2, 2));
    ep_merge_oblist(ibd0, 2) = temp3/2;
    ep_merge_oblist(ibd0, 3) = temp1 - temp3/2;
    ep_merge_oblist(ibd1, 2) = temp3/2;
    ep_merge_oblist(ibd1, 3) = temp2 - temp3/2;
    
    minor_af = min(freq(1:2));
    te2 = minor_af / 10; 
    ep_merge_oblist(ibd2, ibs0) = te2 / 10;
    ep_merge_oblist(ibd2, ibs1) = te2;
    ep_merge_oblist(ibd2, ibs2) = ep_merge_oblist(ibd2, ibs2) * (1 - ms - te2 - te2/10);
    ep_merge_oblist(ibd2, none) = ms;
    
    te1 = minor_af * minor_af / 40;
    ep_merge_oblist(ibd1, ibs0) = te1;
    ep_merge_oblist(ibd1, ibs1) = ep_merge_oblist(ibd1, ibs1) * (1 - ms - te1);
    ep_merge_oblist(ibd1, ibs2) = ep_merge_oblist(ibd1, ibs2) * (1 - ms - te1);
    ep_merge_oblist(ibd1, none) = ms;
    
    ep_merge_oblist(ibd0, ibs0) = ep_merge_oblist(ibd0, ibs0) * (1 - ms - te) + 0.125 * te;
    ep_merge_oblist(ibd0, ibs1) = ep_merge_oblist(ibd0, ibs1) * (1 - ms - te) + 0.500 * te;
    ep_merge_oblist(ibd0, ibs2) = ep_merge_oblist(ibd0, ibs2) * (1 - ms - te) + 0.375 * te;
    ep_merge_oblist(ibd0, none) = ms;
    
else
    error = 1;
    disp('emission states not supported');
    return;
end



% [1p2p, 1p2m, 1m2p, 1m2m] included states
states = [
    1, 1, 1, 1; %ibd0
    2, 1, 1, 1; %ibd1
    3, 1, 1, 1;
    1, 2, 1, 1;
    1, 3, 1, 1;
    1, 1, 2, 1;
    1, 1, 3, 1;
    1, 1, 1, 2;
    1, 1, 1, 3;
    2, 1, 1, 2; %ibd2
    3, 1, 1, 2;
    2, 1, 1, 3;
    3, 1, 1, 3;
    1, 2, 2, 1;
    1, 3, 2, 1;
    1, 2, 3, 1;
    1, 3, 3, 1;
    ];

mapping(1:17) = 1; %ibd0
mapping(2:9) = 2;   %ibd1;
mapping(10:17) = 3; %ibd2;
% count both foreground and background ibd as ibd

ep_oblist(1:length(mapping),1:ncol) = 0;

for i = 1:length(mapping)
    for j = 1:ncol
        ep_oblist(i,j) = ep_merge_oblist(mapping(i), j);
    end
end

ep = ep_oblist;


end


