function [ result error overflow ] = generate_consistent_assignment_likelihood( pairs, triples, range )

error = 0;
overflow = 0;
result = [];

if( isempty(pairs) || isempty(range) || isempty(range.family_range) || isempty(triples) )
    error = 1;
    result = [];
    return;
end



% skip for first child
[map error] = generate_map(range, 0);
if( error ~= 0 )
    disp('error in family structures');
    return;
end
map.relevance = range.allele_source.relevance;


stat_interval = pairs.intervals;
[stat_assignment_all, stat_assignment_result error overflow] = build_global_inheritance_likelihood(pairs, triples, map);
if( error ~= 0 )
    disp('error in generating global IBD');
    return;
end
if( overflow ~= 0 )
    return;
end

stat_interval(:,4:4+length(stat_assignment_result(1,:))-1) = stat_assignment_result;

result.intervals = stat_interval;
result.alleles_all = stat_assignment_all.all;


end

