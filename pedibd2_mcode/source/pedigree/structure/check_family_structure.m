function [error] = check_family_structure(pedigree, families)

% check whether the family indexing is correct, whether parents are all
% proceeding children

error = 0;

if( isempty(families) )
    disp('error in processing family structures');
    error = 1;
    return;
end

[rows, cols] = size(pedigree);

if( rows <= 0 || cols < 12 )
    error = 1;
    disp('error in pedigree structure');
    return;
end

nonrep(1:rows) = 0;

for i = 1:length(families)
    if( isempty(families{i}) )
        disp('error in processing family structures');
        error = 1;
        return;
    end
    [ error ] = check_1processedfamily(families{i});
    if( error ~= 0 )
        disp('error in family structures');
        return;
    end
    
    full_range = families{i}.pedigree_range_full;
    genotyped_range = families{i}.family_range;
    family_structure = families{i}.structure;
    nind = length(full_range);
    if( nind <= 0 )
        error = 1;
        disp('error in processing family structures');
        return;
    end
    [r1, c1] = size( family_structure );
    if( r1 <= 0 || c1 < 12 )
        error = 1;
        disp('error in processing family structures');
        return;
    end
    if( r1 ~= nind )
        error = 1;
        disp('error in processing family structures');
        return;
    end
    if( any(full_range <= 0) || any(full_range > rows) )
        error = 1;
        disp('error in pedigree indexing');
        return; 
    end
    for k = 1:length(full_range)-1
        if( full_range(k+1) ~= full_range(k) + 1 )
            error = 1;
            disp('error in pedigree indexing');
            return;           
        end
    end
    if( any(nonrep(full_range)>0) )
        error = 1;
        disp('error in partition families');
        return;
    else
        % partition
        nonrep(full_range) = 1;
    end
    ref_family = pedigree(full_range,:);
    [r2, c2] = size( ref_family );
    if( nind ~= r1 || nind ~= r2 || c1 ~= c2 || c1 < 12 )
        error = 1;
        disp('error in processing family structures');
        return;
    end
    % check the imported fields are correct
    if( nnz(family_structure(1:nind,1:12)-ref_family(1:nind,1:12)) > 0 )
        error = 1;
        disp('error in processing family structures');
        return;
    end    
    if( ~isempty(genotyped_range) )
        if( ~all(family_structure(genotyped_range,7) == 1) )
            error = 1;
            disp('error in processing family structures');
            return;         
        end
        if( length(unique(genotyped_range)) ~= length(genotyped_range) )
            error = 1;
            disp('error in processing family structures');
            return;          
        end
        if( max(genotyped_range) > nind || min(genotyped_range) < 1 )
            error = 1;
            disp('error in processing family structures');
            return;              
        end
    end
    % check family structure
    for j = 1:nind
        if( family_structure(j,2) ~= j )
            error = 1;
            disp('error in processing family structures');
            return; 
        end
        if( family_structure(j,3) > j || family_structure(j,4) > j )
            % parents not proceeding children;
            error = 1;
            disp('error in processing family structures');
            return; 
        end
    end
end

if( any( nonrep <= 0 ) )
    error = 1;
    disp('error in partitioning family structures');
    return;
end

end







