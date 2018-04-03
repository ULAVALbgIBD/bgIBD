function error = check_assignment(alleles, num_ind, genotyped_range)

% check if all genotyped individuals are assigned

error = 0;
[rows, cols] = size(alleles);

if( rows ~= num_ind || cols ~= 2 )
    error = 1;
    return;
end


ids = unique(alleles);

ids = ids(ids>0);

if( ~isempty(ids) )
    if( max(ids) > num_ind )
        error = 1;
        return;
    end
end

% mapped to holes
if( any(any(alleles(ids,1:2)==0)) )
    error = 1;
    return;
end

% genotyped individual no IBD generated
if( any( any(alleles(genotyped_range,1:2)==0, 2), 1 ) )
    error = 1;
    return;
end


end


