
function [epl epl_log error] = emission_pair(pm, family_genotype, sd1, sd2, family_structure, emission_option)

error = 0;
epl = [];
epl_log = [];

sampled_af = pm.sampled_af;
pm_missing_rate = pm.missing_rate;
pm_typing_error = pm.typing_error;
sampled_markerlist = pm.sampled_markerlist;



[d1, d2, d3] = size(family_genotype);
[nind, cols] = size(family_structure);
[nloc, c] = size( sampled_markerlist );

if( d2 <= 0 || d2 ~= nind )
    error = 1;
    disp('error in family structures');
    return;
end

if( d1 <= 0 || d1 ~= nloc )
    error = 1;
    disp('error in genotype data');
    return;
end

if( c ~= 2 || d3 ~= 2 )
    error = 1;
    disp('error in generating emission');
    return;
end

father_id1 = family_structure(sd1, 3);
mother_id1 = family_structure(sd1, 4);
father_id2 = family_structure(sd2, 3);
mother_id2 = family_structure(sd2, 4);

genotyped = family_structure(1:nind, 7) == 1;

father_geno = zeros(nloc, 2);
mother_geno = zeros(nloc, 2);

sibpair = 0;

if( father_id1 ~= 0 && father_id1 == father_id2 && genotyped(father_id1) )  
    father_geno(1:nloc,1:2) = family_genotype(1:nloc, father_id1, 1:2);
    sibpair = 1;
end

if( mother_id1 ~= 0 && mother_id1 == mother_id2 && genotyped(mother_id1) )  
    mother_geno(1:nloc,1:2) = family_genotype(1:nloc, mother_id1, 1:2);
    sibpair = 1;
end


[epl epl_log error] = generate_emission_all_loci(sampled_af, pm_missing_rate, pm_typing_error, sampled_markerlist, father_geno, mother_geno, sibpair, emission_option);
if( error ~= 0 )
    disp('error in generating emission states');
    return;
end


end


