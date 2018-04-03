

function [input output pedigree parameters error] = ...
    dbFamily(pedigree_info, ...
    all_data, map, option)

input = [];
output = [];
pedigree = [];
parameters = [];
error = 0;

[input pedigree expanded_genotype error] = ...
    process_input(pedigree_info, ...
    all_data, map.physical_map);
    
if( error == 1 )
    disp('error in input file: skipped');
    return;
end

% mendelian erroneous loci removed
[input.family_genotype parameters error]= ...
    generate_allele_frequency(pedigree, ...
    input, map, expanded_genotype);

if( error == 1 )
    disp('error in generating parameters');
    return;
end


[input output error] = ...
    process_allfamilies(input, parameters, option);


if( error == 0 && option(3) == 1 )
    [output.combined_score, error] = ...
        score_allfamilies(output.merged_assignment, ...
        input.family_range);
    
end

end











