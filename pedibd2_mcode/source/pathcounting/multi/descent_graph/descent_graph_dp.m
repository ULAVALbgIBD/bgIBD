
%make sure parents are relative indexed by id
% ancestor and related are all allele based

% genF

% column 1: number of each term
% coefficient

% column 2~v+1: (identity) whether exist or not
% generic, only indicates which one equals which one

% column v+2 ~ up: terms and their powers
% terms encodes a sub-root of the inheritance tree
% from current alleles to reach the identity state 
% specified in 2~v+1
% each term encode the length between a sub-root to
% its parent, sub-root


function [output_genF] = descent_graph_dp(parental_al, ancestor, related, vec )

    global all_inheritance;
    

    
    % check whether the vec is already calculated
    [index, output_genF] = lookup_inheritance(vec);
    if( index ~= 0 )
        return;
    end
    
    % find a slot to do recurrence
    slots = length(vec);

    is_ances = false(slots,1); %not ancestor of any individual
    have_rel = false(slots,1); %do not have any relatives
    for i = 1:slots
        for j = 1:slots
            if( ancestor(vec(i),vec(j)) ) % i is an ancestor of j
                is_ances(i) = true;
            end
        end
    end
    for i = 1:slots
        for j = 1:slots
            if( vec(i) ~= vec(j) && related(vec(i),vec(j)) )
                have_rel(i) = true;
                have_rel(j) = true;
            end
        end
    end
    slot_of_intes = 0;
    for i = 1:slots
        % pick an allele to do recurrence
        if( ~is_ances(i) && have_rel(i) )
            slot_of_intes = i;
            break;
        end
        % terminal states are all identical alleles
        % or unrelated alleles
    end
    if( slot_of_intes == 0 )
        % arrived at a terminal state
        
        output_genF(1,1) = 1;
        output_genF(1,2:2+slots-1) = identity_coding(vec);  
        %transform to generic states
        % final identity state of the alleles, generic terminal states only
        
        max_power_loc = slots + 1 + term_coding(1:length(vec));
        % set the maximum number of bits required for storing powers
        % all alleles appear, is the largest possible number
        
        output_genF(1,slots+1+1:max_power_loc) = 0;
    else
        al_of_intes = vec(slot_of_intes);
        id = find( vec == al_of_intes );
        % slots of alleles to be replaced in the recurrence
        % if multiple do recurrence at the same time
        vecp = vec;
        vecp(id) = parental_al(al_of_intes, 1);
        [genFp] = descent_graph_dp(parental_al, ancestor, related, vecp);
        vecm = vec;
        vecm(id) = parental_al(al_of_intes, 2);
        [genFm] = descent_graph_dp(parental_al, ancestor, related, vecm);
        
        output_genF = descent_graph_combine(genFp, genFm);
        power_loc = (slots+1) + term_coding(id);
        % increment the correponding power by 1
        output_genF(:,power_loc) = output_genF(:,power_loc) + 1;
    end
    
    % update all_inheritance, by adding the new generating function
    [~] = add_inheritance(vec, output_genF);
    
end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

