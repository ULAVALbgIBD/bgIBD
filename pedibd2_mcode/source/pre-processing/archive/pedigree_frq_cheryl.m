function [ all_pedigree_frq ] = pedigree_frq_cheryl(pedigree, all_data, sampled_frq, all_families)


frq = sampled_frq(:,1);
num = length(all_data(:,1));
loci = length(all_data(1,7:end))/2;

gen(1:num,1:loci,1:2) = 0;
gen(1:num,1:loci,1) = all_data(1:num,7:2:end);
gen(1:num,1:loci,2) = all_data(1:num,8:2:end);

% all id mapped to the correct location in the genotype file
% untyped individuals are paddled with all 0's

all_pedigree_frq = generate_pedigree_frq(pedigree, gen, frq, all_families);

end

