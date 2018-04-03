
%make sure parents are relative indexed by id
% ancestor and related are all allele based

% genF
% column 1~v: (identity) whether exist or not
% column v+1: number
% column v+2 ~ up: terms and their powers



function output_genF = descent_graph( parental_al, ancestor, related, vec )
    % find an slot to do recurrence
    slots = length(vec);

    is_ances(1:slots) = 0; %not ancestor of any individual
    have_rel(1:slots) = 0; %do not have any relatives
    for i = 1:slots
        for j = 1:slots
            if( ancestor(vec(i),vec(j)) == 1 ) % i is an ancestor of j
                is_ances(i) = 1;
            end
        end
    end
    for i = 1:slots
        for j = 1:slots
            if( vec(i) ~= vec(j) && related(vec(i),vec(j)) == 1 )
                have_rel(i) = 1;
                have_rel(j) = 1;
            end
        end
    end
    slot_of_intes = 0;
    for i = 1:slots
        % pick an allele to do recurrence
        if( is_ances(i) == 0 && have_rel(i) == 1 )
            slot_of_intes = i;
            break;
        end
    end
    if( slot_of_intes == 0 )
        % arrived at a terminal state
        output_genF(1,1) = 1;
        output_genF(1,2:2+slots-1) = identity_coding(vec);
        max_power_loc = slots + 1 + term_coding(1:length(vec));
        output_genF(1,slots+1+1:max_power_loc) = 0;
    else
        al_of_intes = vec(slot_of_intes);
        id = find( vec == al_of_intes );
        % slots of alleles to be replaced in the recurrence
        vecp = vec;
        vecp(id) = parental_al(al_of_intes, 1);
        genFp = descent_graph( parental_al, ancestor, related, vecp);
        vecm = vec;
        vecm(id) = parental_al(al_of_intes, 2);
        genFm = descent_graph( parental_al, ancestor, related, vecm);
        
        output_genF = descent_graph_combine(genFp, genFm);
        power_loc = (slots+1) + term_coding(id);
        output_genF(:,power_loc) = output_genF(:,power_loc) + 1;
    end
end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

