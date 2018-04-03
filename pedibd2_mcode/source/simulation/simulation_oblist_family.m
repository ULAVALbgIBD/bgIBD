
function [family, random_ibd, random_genotype, ref_genotype] ...
    =  simulation_oblist_family(f, pedigree_all_missing, sm, parameters)


% generate genotypes for a family


pm = parameters;

temp =  pedigree_all_missing(:,1) == f ;
family = pedigree_all_missing(temp,:);
markers = length(pm.sampled_markerlist(:,1));

random_ibd = generate_random_ibd(family, markers, sm.recfrac);
random_ibd = condense_random_ibd(random_ibd, markers);

ref_genotype = generate_random_genotype(random_ibd, sm.haplotypes, sm.founder_source);
random_genotype = generate_random_missing_error(ref_genotype, pm.missing_rate, pm.typing_error);

 
end


