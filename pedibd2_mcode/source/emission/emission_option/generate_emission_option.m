function emission_option = generate_emission_option(option)


% option 1, full genotype
% option 2, ibs
% opition3, split ibds 2, heterozygous and homozygous


[emission_option.allele2genotype] = allele_map();
[emission_option.genotype2pair, emission_option.pair_missing_code] = genotype_code(emission_option.allele2genotype);
emission_option.pair2genotype = genotype_reverse_code(emission_option.genotype2pair);
[emission_option.pair2ob, emission_option.ob_missing_code] = pair2obMap(emission_option, option);

emission_option.sibling_ibd_states = ibd_interface_siblingpair();
emission_option.nuclear_emission = generate_nuclear_family(emission_option);
emission_option.emission_option = option;


end