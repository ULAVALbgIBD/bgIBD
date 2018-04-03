function emission =  generate_nf_emission(genotype1, genotype2, emission_option, a2i)

    [trans genotype] = enumerate_transmission(genotype1, genotype2);
    emission = count_configuration(trans, genotype, emission_option, a2i);

end