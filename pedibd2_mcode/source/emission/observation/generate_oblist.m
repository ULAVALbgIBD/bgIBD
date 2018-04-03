

function [ob emission_option error] = generate_oblist(family_genotype, sp1, sp2, sampled_markerlist, ob_option)

error = 0;
ob = [];
global debug_mode;

emission_option = generate_emission_option(ob_option);

[rows, cols, d3] = size(family_genotype);
if( cols <= 0 || sp1 > cols || sp2 > cols || d3 ~= 2 )
    error = 1;
    disp('error in genotype data');
    return;
end

[nloc, c] = size( sampled_markerlist );

if( rows ~= nloc || nloc <= 0  || c ~= 2 )
    error = 1;
    disp('error in genotype data');
    return;
end

geno1(1:nloc,1:2) = family_genotype(1:nloc, sp1, 1:2);
geno2(1:nloc,1:2) = family_genotype(1:nloc, sp2, 1:2);


input_genotype_pair = generate_genotype(geno1, geno2, emission_option);

if( isempty(input_genotype_pair) )
    error = 1;
    disp('error in genotype data, data not coded with 1 2 0');
    ob = [];
    return;
end

if( length(input_genotype_pair) ~= nloc )
    error = 1;
    disp('error in genotype data');
    return;
end

ob = pair2ob(input_genotype_pair, emission_option.pair2ob);



end


























