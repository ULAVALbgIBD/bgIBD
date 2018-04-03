function emission = count_configuration(trans, genotype, emission_option, a2i)

genotype2pair = emission_option.genotype2pair;

n_rows = max(max(max(max(a2i))));
n_col = max(max(genotype2pair));


emission(1:n_rows, 1:n_col) = 0;

for i = 1:length(trans(:,1))
    temp1 = a2i(trans(i,1),trans(i,2),trans(i,3),trans(i,4));
    temp2 = code_pair(genotype(i,:), emission_option);
    emission(temp1,temp2) = emission(temp1,temp2)+1;
end

for i = 1:n_rows
    emission(i,:) = emission(i,:)./sum(emission(i,:));
end

end