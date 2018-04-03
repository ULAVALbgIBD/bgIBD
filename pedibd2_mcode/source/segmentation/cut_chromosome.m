function [pairs triples_view error] ...
    = cut_chromosome(family_range, kinship, kinship2ex, ...
    family_genotype, viterbi, posterior, posteriorIBD1, oblist, parameters)

    error = 0;
    pairs = [];
    triples_view = [];
    
    time = cputime;
    
    [pairs, corrected_viterbi, error] = generate_consistent_pairs(oblist, viterbi, posterior, posteriorIBD1, kinship2ex, family_range, parameters);
    if( error ~= 0 )
        disp('error in cut the chromosome');
        return;
    end
    [~, triples_view, error] = valid_triples(family_range, kinship, kinship2ex, family_genotype, corrected_viterbi, pairs.intervals);
    if( error ~= 0 )
        disp('error in cut the chromosome');
        return;
    end
  
end




