
% exact allele count for allele 1 and 2, first field for allele 1

function [ pm_freq error ] = allele_frq( genotype, member_weight, num_markers )

error = 0;
pm_freq = [];

[nind c] = size(member_weight);
if( nind <= 0 || c ~= 1 )
    error = 1;
    disp('error in family structures');
    return;
end

[nrows, ncols, c] = size(genotype);
if( nrows ~= num_markers || ncols ~= nind || c ~= 2 )
    error = 1;
    disp('error in generating allele frequency');
    return;
end


allele1 = sum(genotype == 1, 3);
allele2 = sum(genotype == 2, 3);
allele0 = sum(genotype == 0, 3);

all_freq = zeros(num_markers, 3);
all_freq(1:num_markers, 1) = allele1 * member_weight;
all_freq(1:num_markers, 2) = allele2 * member_weight;
all_freq(1:num_markers, 3) = allele0 * member_weight;

freq12 = zeros(num_markers, 1);
freq123 = zeros(num_markers, 1);
freq12(1:num_markers) = sum(all_freq(1:num_markers, 1:2), 2);
freq123(1:num_markers) = sum(all_freq(1:num_markers, 1:3), 2);

all_af = zeros(num_markers, 3);
all_af(freq12 ~= 0, 1) = all_freq(freq12 ~= 0, 1)./freq12(freq12 ~= 0);
all_af(freq12 ~= 0, 2) = all_freq(freq12 ~= 0, 2)./freq12(freq12 ~= 0);
all_af(freq12 ~= 0, 3) = all_freq(freq12 ~= 0, 3)./freq123(freq12 ~= 0);

all_af(freq12 == 0, 1:2) = 0.5;
all_af(freq12 == 0, 3) = 0;

pm_freq.all_freq = all_freq;
pm_freq.all_af = all_af;


end

