function [output_alleles, intervals, error] = check_allsegments(assignment, family)

error = 0;
output_alleles = [];
intervals = [];

if( isempty(family) )
    error = 1;
    disp('error in family structure');
    return;
end

if( isempty(assignment) )
    error = 1;
    disp('error in global IBD');
    return;
end

full_range = family.pedigree_range_full;
genotyped_range = family.family_range;

nIND = length(full_range);

error = check_1processedfamily(family);
if( error ~= 0 )
    error = 1;
    disp('family ', num2str(family.family_id), ': error in family structure');
    return;
end

% check whether global IBD segments are valid
% check whether genotype number match family members

intervals = assignment.intervals(:,1:2);
[nSEG fields] = size(intervals);

if( fields < 2 )
    error = 1;
    disp('error in global IBD');
    return;
end

if( nSEG <= 0 )
    error = 1;
    disp('error in global IBD')
    return;
end

% check global IBD
alleles_all = assignment.alleles_all;

if( ndims(alleles_all) ~= 3 || any(size(alleles_all) ~= [nSEG, nIND, 2]) )
    error = 1;
    disp('error in global IBD, not all nSEG processed');
    return;
end


for i = 1:nSEG
    seg_assignment = reshape(alleles_all(i, 1:nIND, 1:2), [nIND, 2]);
    [error] = check_1segment(seg_assignment, family.structure);
    if( error ~= 0 )
        error = 1;
        disp('error in global IBD');
        return;
    end
end

output_alleles = zeros(nSEG, nIND, 2);
for i = 1:nSEG
    seg_assignment = reshape(alleles_all(i, 1:nIND, 1:2), [nIND, 2]);
    error = check_assignment(seg_assignment, nIND, genotyped_range);
    if( error ~= 0 )
        error = 1;
        disp(['error in global IBD, segment ', num2str(i)]);
        return;
    end
    output_alleles(i, 1:nIND, 1:2) = seg_assignment;
end

output_alleles = alleles_all;

end

















