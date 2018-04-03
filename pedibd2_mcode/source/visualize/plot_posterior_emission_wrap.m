function plot_posterior_emission_wrap(input, parameters, f, id1, id2, limit)

global debug_mode;
debug_mode = 1;

temp = input.family_range{f}.structure(input.family_range{f}.family_range,8);
in1 = find(temp == id1);
in2 = find(temp == id2);
if( length(in1) ~= 1 || length(in2) ~= 1 || id1 == id2 )
    return;
end

s_id = input.family_range{f}.reverse_pairs(in1, in2);
input_pair = input.family_range{f}.pairs(s_id,1:2);

coor = parameters.sampled_markerlist(:,1);
if( length(limit) ~= 2 )
    limit = [1,length(parameters.sampled_markerlist(:,1))];
else
    if( limit(1) > limit(2) || limit(1) < 1 || limit(2) > length(parameters.sampled_markerlist(:,1)) )
        return;
    end
end


all_inheritance = input.all_inheritance{f};
family_range = input.family_range{f};
family_genotype = input.family_genotype{f};

analyse_result(family_range.structure, all_inheritance, parameters, family_genotype, input_pair, limit)

clear global debug_mode;

end