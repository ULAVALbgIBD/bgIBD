function output_nf = generate_nuclear_family(emission_option)

output_nf = [];

allele2genotype = emission_option.allele2genotype;
genotype2pair = emission_option.genotype2pair;
a2i = emission_option.sibling_ibd_states;


for i1 = 1:2
    for i2 = 1:2
        for j1 = 1:2
            for j2 = 1:2
                genotype1 = [i1,i2];
                genotype2 = [j1,j2];
                emission = generate_nf_emission(genotype1,genotype2,emission_option,a2i);
                pair = code_pair([i1,i2,j1,j2], emission_option);
                output_nf{pair} = emission;
            end
        end
    end
end




end