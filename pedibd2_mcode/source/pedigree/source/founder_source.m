function [ source error ] = founder_source( all_families )

% parents must proceed children in order


count = 0;
error = 0;


% direct gives the nearest genotyped alleles, if no genotyped, tracing
% until founders

if( isempty(all_families) )
    error = 1;
    disp('no family structure found');
    return;
end

maxid = 0;
for i = 1:length(all_families)
    global_scope = all_families{i}.pedigree_range_full;
    if( isempty(global_scope) )
        error = 1;  
        disp('no member in family ', num2str(all_families{i}.family_id));
        return;
    end
    if( any( global_scope <= 0 ) )
        error = 1;
        disp('error in family partitioning');
        return;
    end
    maxid = max([global_scope,maxid]);
end

source = zeros(maxid,1);
for i = 1:length(all_families)
    global_scope = all_families{i}.pedigree_range_full;
   
    weight = all_families{i}.founder_source;
    if( any(source(global_scope)>0) )
        error = 1;
        disp('error in partitioning family');
        return;
    end
    source(global_scope) = weight;
    % record in the pedigree-wise list
end



end

