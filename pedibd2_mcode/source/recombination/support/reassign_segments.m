function [new_alleles recombination error] = reassign_segments(assignment, family, parameters)

error = 0;
new_alleles = [];
recombination = [];

[alleles_all, intervals, error] = check_allsegments(assignment, family);
if( error ~= 0 )
    disp('error in processing global IBD');
    return;
end

disp(['family ', num2str(family.family_id), ': inferring recombination positions ...']);

[sibship error] = sibship_isomorphism(assignment, family);
if( error ~= 0 )
    disp('error in processing global IBD');
    return;
end


[recombination error] = extract_recombination(sibship, family.structure, intervals, parameters);
if( error ~= 0 )
    disp('error in extracting recombination points');
    return;
end

% some segments are not reordered because they are not valid
% recombination too close to each other
[new_alleles error] = reduce_recombination(family.structure, recombination.valid, alleles_all);
if( error ~= 0 )
    disp('error in global IBD');
    return;
end
assignment.alleles_all = new_alleles;

[new_alleles, ~, error] = check_allsegments(assignment, family);
if( error ~= 0 )
    disp('error in processing global IBD');
    return;
end



end



