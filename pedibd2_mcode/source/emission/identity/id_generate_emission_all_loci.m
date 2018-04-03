function epl = id_generate_emission_all_loci(num, parameters)

ct = 0;
t = cputime;

epl = [];

[is transition] = hasse_matrix(num);

for i = 1:length(parameters.sampled_af(:,1));
    emi = id_generate_emission(is, parameters.sampled_af(i,1:3), parameters.missing_rate, parameters.typing_error);    
    epl{i} = id_expand_emission(is, emi);
end


ct = ct + cputime - t

end


