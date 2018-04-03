

function oblist = oblist_sc(all_data, sp, pm)


input_genotype_data_sc = genotype_sc(all_data(sp, 7:end));
oblist = input_genotype_data_sc(pm.sampled_markerlist(:,1));


end




















