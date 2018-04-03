
function [family_genotype parameters error] = generate_allele_frequency(pedigree, input, map, expanded_genotype)
    
    global debug_mode;
    error = 0;
    
    disp('pre-processing input data...');
    if( debug_mode == 1 )
        disp('estimating allele frequency ...');
    end    
    time = cputime;
    
    [parameters error] = generate_parameters(pedigree.founder_source, map, expanded_genotype);
    if( error ~= 0 )
        return;
    end
    
    
    [family_genotype error] = divide_families(pedigree.structure, expanded_genotype, input.family_range);
    if( error ~= 0 )
        return;
    end
    
    display(['allele frequency costs time: ', num2str(cputime - time), ' seconds']); 
    
end


