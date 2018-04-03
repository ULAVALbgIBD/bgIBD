

function [output_pedigree, error] = pedigree_structure(input_structure)

error = 0;

output_pedigree.coding.sex = 5;
output_pedigree.coding.disease_status = 6;
output_pedigree.coding.genotyped = 7;
output_pedigree.coding.original_id = 8;
output_pedigree.coding.generation = 9;
output_pedigree.coding.original_family = 12;

%not every individual in a family is connected

%input_pedigree
%column 7 indicates whether the individual is missing
%column 8~ other nFIELDS to be linked using column 11


%output_pedigree
%column 1 original family id
%column 2~4 newly generated sequential order
%column 7 indicates whether the individual is missing
%column 8 gives its original id
%column 9 generation
%column 10 inner family, to be generated later in family_size()
%column 11 gives sequential order in original input_structure
%column 12 gives the original family id, currently same as 1

%all_data and pedigree_out do not have the same order.

%make founder list
%duplicate ids in a family

 
    % temp1, order in the original input_structure
    % temp, family structure of temp1
    % temp3, newly generated family structure
    % pedigree_out, stack temp3 of each family

    global_count = 0;
    pedigree_all = [];
    generation = 9;
    
    f_list = unique(input_structure(:,1));
    
    for i = f_list'
        
        family_index = find( input_structure(:,1) == i );

        if( isempty(family_index) )
            continue;
        end

        family = input_structure(family_index, 1:7);
        
        [processed_family error] = process1familystructure(family);
        if( error ~= 0 )
            disp('error in processing families');
            return;
        end
               
        [nIND, nFIELDS] = size(processed_family);
        if( nIND < length(family_index) || nFIELDS ~= 12 )
            error = 1;
            disp('error in processing families');
            return;
        end
        for j = 1:nIND
            if( processed_family(j,11) ~= 0 )
                if( processed_family(j,11) > length(family_index) )
                    error = 1;
                    disp('error in processing families');
                    return;
                else
                    processed_family(j,11) = family_index(processed_family(j,11));
                end
            end
        end
        
        
        % stack pedigrees together
        pedigree_all(global_count+1:global_count+nIND, 1:12) = processed_family;
        global_count = global_count + nIND;
        
    end
    
    output_pedigree.structure = pedigree_all;


end














