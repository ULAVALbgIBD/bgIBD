function [config error] = generate_family_display(alleles_all, family)

config = [];
error = 0;

[nrow, ncol] = size(alleles_all);
if( nrow <= 0 )
    disp('error in the allele assignment');
    error = 1;
    return;
end
if( ncol ~= 2 )
    disp('error in the allele assignment');
    error = 1;
    return;
end


typed_range = family.family_range;
full_range = family.pedigree_range_full;

structure = family.structure;
sex = family.coding.sex;
disease = family.coding.disease_status;
id_col = family.coding.original_id;

map2id = structure(:,id_col);

if( length(alleles_all(:,1)) ~= length(full_range) )
    disp('error in the allele assignment');
    error = 1;
    return;
end
if( isempty(typed_range) )
    % no genotyped individuals, cannot display
    return;
end


alleles = alleles_all(typed_range,1:2);

if( any( abs(alleles) > length(map2id) ) )
    disp('error in the allele assignment');
    error = 1;
    return;
end
if( any( alleles == 0 ) )
    disp('not all individuals assigned alleles');
    error = 1;
    return;
end


% assign external id to the allele map of a family
% unique_alleles order will be used as coordinates later
unique_alleles = unique(alleles);
[temp, order] = sort(abs(unique_alleles+0.1), 'ascend'); % minus, paternal allele will proceed
unique_alleles = unique_alleles(order);

for i = 1:length(typed_range)
    id_string{i} = num2str(map2id(typed_range(i)));
    for j = 1:2
        a = alleles(i,j);
        if( a > 0 )
            allele_string{i,j} = [num2str(map2id(a)), 'm'];
        else
            allele_string{i,j} = [num2str(map2id(-a)), 'p'];
        end
    end
end

for i = 1:length(unique_alleles)
    a = unique_alleles(i);
    if( a > 0 )
        allele_name{i} = [num2str(map2id(a)), 'm'];
    else
        allele_name{i} = [num2str(map2id(-a)), 'p'];
    end
end


config.affected = (structure(typed_range,disease) == 2);
config.non_affected = (structure(typed_range,disease) == 1);
config.unknown = (structure(typed_range,disease) == 0); 
config.ismale = (structure(typed_range,sex) == 1);
config.isfemale = (structure(typed_range,sex) == 2);
config.allele_string = allele_string;
config.alleles = alleles;
config.unique_alleles = unique_alleles;
config.allele_name = allele_name;
config.id_string = id_string;
config.num_alleles = length(unique_alleles);
config.num_ind = length(typed_range);

config.group1 = config.affected & config.ismale;
config.group2 = config.affected & config.isfemale;
config.group3 = config.non_affected & config.ismale;
config.group4 = config.non_affected & config.isfemale;
config.group5 = config.unknown & config.ismale;
config.group6 = config.unknown & config.isfemale;


end








