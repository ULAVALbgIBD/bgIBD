function [ pairs reverse_pair ] = generate_pairlist( family_range, child_range)

count = 0;

input_pairlist = [];
input_pairlist_reverse = [];

for i = 1:length(family_range)
    for j = i+1:length(family_range)
        count = count + 1;
        input_pairlist(count, 1:2) = [family_range(i), family_range(j)];
        input_pairlist_reverse(i,j) = count;
        input_pairlist_reverse(j,i) = count;
        % order is ignored, may not be correct for allele specific IBD
    end
end


pairs = input_pairlist;  % projects to the external ids of the whole pedigree
reverse_pair = input_pairlist_reverse;



end

