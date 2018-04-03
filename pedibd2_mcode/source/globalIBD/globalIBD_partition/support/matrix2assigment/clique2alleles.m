function [assignment error] = clique2alleles(clique_config, family_range, kinship2ex)
    % assign ancestral allele origin
    % ancestral allele must be real local id in whole range of the family
    % 1:nIND
    error = 0;
    assignment = [];
    
    clique_assignment = clique_config.assignment;
    rep_list = clique_config.rep_list;   
    
    genotyped = family_range.family_range;
    nIND = length(family_range.pedigree_range_full);
    nGENO = length(genotyped);
    family = family_range.structure;
    [r, fields] = size(family);
    if( r <= 0 || r ~= nIND || fields < 12 )
        error = 1;
        return;
    end
    
    nchr = length(rep_list);
    if( nchr <= 0 )
        error = 1;
        disp('error in global IBD');
        return;
    end
    
    [rows cols] = size(clique_assignment);
    if( isempty(clique_assignment) || rows <= 0 || rows ~= nchr || cols ~= nGENO )
        error = 1;
        disp('global IBD inference error by diploidy removal');
        return;
    end
    
    if( nGENO > nIND || any(genotyped > nIND) )
        error = 1;
        disp('error in family structure');
        return;
    end
    
    if( length(unique(rep_list)) ~= nchr )
        error = 1;
        disp('index reuse, pairwise ibd inconsistency');
        return;
    end
    
    if( any(size(kinship2ex) ~= [nIND, nIND, 2, 2]) )
        error = 1;
        disp('error in family structure');
        return;
    end    
        
    compact_assignment(1:nGENO,1:2) = 0;
  
    for i = 1:nGENO
        temp = find( clique_assignment(1:nchr,i) > 0 );
        if( length(temp) ~= 2 )
            error = 1;
            disp('error in graph partition, not separating correctly');
            return;
        end
        compact_assignment(i,1:2) = temp;
    end
        
    % switch paternal and maternal alleles as forced by pedigree
    assignment = zeros(nIND,2);
    assignment(genotyped,1:2) = compact_assignment(1:nGENO, 1:2);
    [assignment, ~, error] = phase_alleles(assignment, family, kinship2ex);
    if( error ~= 0 )
        disp('error in determining parental sources');
        return;
    end
    
end










