function [ source error ] = allele_source( all_families )

% parents must proceed children in order

count = 0;
error = 0;


num_ind = 0;
maxid = 0;
for i = 1:length(all_families)
    global_scope = all_families{i}.pedigree_range_full;
    if( any( global_scope <= 0 ) )
        error = 1;
        disp('error in family partitioning');
        return;
    end
    num_ind = num_ind + length(all_families{i}.pedigree_range_full);
    maxid = max([global_scope,maxid]);
end

if( maxid ~= num_ind )
    error = 1;
    disp('error in family partitioning');
    return;
end

relevance = zeros(num_ind,1);

% direct gives the nearest genotyped alleles, if no genotyped, tracing
% until founders

for i = 1:length(all_families)
    global_scope = all_families{i}.pedigree_range_full;
    
    local_relevance = all_families{i}.allele_source.relevance;
    
    if( error ~= 0 )
        disp('error in processing families');
        error = 1;
        return;
    end
    if( length(local_relevance) ~= length(global_scope) )
        disp('error in processing families');
        error = 1;
        return;
    end
    if( any(relevance(global_scope)>0) )
        disp('error in partitioning families');
        error = 1;
        return;
    end
    
    relevance(global_scope) = local_relevance;
    
end

source.relevance = relevance;

end

