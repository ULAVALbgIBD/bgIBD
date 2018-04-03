



function [output_all_inheritance error] = descent_graph_2alleles_all(family_range)

    global all_inheritance;
    global echo;
    
    error = 0;
    output_all_inheritance = [];

    pedigree_range_full = family_range.pedigree_range_full;
    family_structure = family_range.structure;
    
    [rows, cols] = size( family_structure );
    if( rows <= 0 || rows ~= length(pedigree_range_full) || cols < 12 )
        error = 1;
        disp('error in family structure');
        return;
    end
    nIND = rows;
    
    
    % family should be strictly indexed by the id
    family = family_structure(1:nIND,3:4);
    print_id = family_structure(1:nIND, 8);
    print_f = family_structure(1, 12);
    
    fast_index = zeros(2*nIND, 2*nIND);
       
    parental = zeros(2*nIND, 2);
    ancestor = false(2*nIND, 2*nIND);
    related = false(2*nIND, 2*nIND);
    
    % expand relationship to 2*nIND bitmap
    for i = 1:length(family(:,1))
        if( family(i,1) ~= 0 )
            %paternal allele source, two alleles of father
            parental(i*2-1, 1:2) = [family(i,1)*2-1, family(i,1)*2];
            ancestor(family(i,1)*2-1, i*2-1) = true;
            ancestor(family(i,1)*2, i*2-1) = true;
        else
            parental(i*2-1, 1:2) = 0;
        end
        if( family(i,2) ~= 0 )
            %maternal allele source, two alleles of mother
            parental(i*2, 1:2) = [family(i,2)*2-1, family(i,2)*2];
            ancestor(family(i,2)*2-1, i*2) = true;
            ancestor(family(i,2)*2, i*2) = true;
        else
            parental(i*2, 1:2) = 0;
        end
    end
    
    % nIND * 2 iteration, to make sure all relationship updated
    for k = 1 : (nIND * 2)
        for i = 1 : (nIND * 2)
            for j = 1 : (nIND * 2)
                if( ancestor(i, j) )
                    ancestor(i, :) = ancestor(i, :) | ancestor(j, :);
                end
            end
        end
    end
    
    for i = 1 : (nIND * 2)
        for j = 1 : (nIND * 2)
            if( any(ancestor(:,i) & ancestor(:,j), 1) || ancestor(i,j) )
                related(i, j) = true;
            end
        end
    end
    
    
    all_inheritance.list = [];
    all_inheritance.content = [];
    
    % build a fast index, applies to two alleles only, more alleles will
    % cost too much memory
    
    all_inheritance.fast_index = fast_index;    
    
    for id1 = 1:nIND
        for id2 = 1:nIND
            % consider all individuals
            if( echo == 1 )
                disp(['          family ', num2str(print_f), ': calculating kinship: ', num2str(print_id(id1)), ', ', num2str(print_id(id2))]); 
            end
            vec = [id1*2-1,id2*2-1];
            [~] = descent_graph_dp(parental, ancestor, related, vec );
            vec = [id1*2,id2*2-1];
            [~] = descent_graph_dp(parental, ancestor, related, vec );
            vec = [id1*2-1,id2*2];
            [~] = descent_graph_dp(parental, ancestor, related, vec );
            vec = [id1*2,id2*2];
            [~] = descent_graph_dp(parental, ancestor, related, vec );
        end
    end
          
    output_all_inheritance = all_inheritance;
    
    clear global all_inheritance;
    
end









