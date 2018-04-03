function epl = id_generate_emission(is, all_af, mr, te)


    ms = all_af(3);
    ms = mr;

    e = id_condensed_generate_emission_genotype(all_af(1:2), is, ms, te);
    
    epl = id_condensed_generate_emission_ibs_genotype(e);


end


