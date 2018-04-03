function pair = code_pair(genotype, emission_option)

allele2genotype = emission_option.allele2genotype;
genotype2pair = emission_option.genotype2pair;
pair_missing = emission_option.pair_missing_code;

pair = [];

if( genotype(1) == 0 || genotype(2) == 0 || genotype(3) == 0 || genotype(4) == 0 )
    pair = pair_missing;
else

    [nrow, ncol] = size(allele2genotype);
    
    if( ~all(ismember(genotype,[1:nrow])) )
        disp('error in genotype coding');
    end

    genotype1 = allele2genotype(genotype(1),genotype(2));
    genotype2 = allele2genotype(genotype(3),genotype(4));
    pair = genotype2pair(genotype1, genotype2);

end

end