function ep_merge_oblist = emission_genotype2ob(ep_merge_genotype, emission_option)

ep_merge_oblist = [];

map = emission_option.pair2ob;
[nrow, ncol] = size(ep_merge_genotype);
if( nrow <= 0 || ncol <= 0 )
    disp('error in genotype emission probability');
end


ep_merge_oblist(1:nrow, 1:max(map)) = 0;


for i = 1:nrow
    for j = 1:ncol
        temp = map(j);
        ep_merge_oblist(i,temp) = ep_merge_oblist(i,temp) + ep_merge_genotype(i,j);
    end
end


end


