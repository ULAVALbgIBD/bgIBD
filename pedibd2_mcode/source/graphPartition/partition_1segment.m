function [alleles_all error] = ...
    partition_1segment(vit, ml, ...
    kinship2ex, ...
    singleLINEAL, ...
    posIBD1, pair_pos, ...
    triples, ...
    range, verify, debug)

%%
    alleles_all = [];
    error = 0;
    % pick one choice of viterbi decoding or maximum likelihood
    if( isempty(range) )
        error = 1;
        return;
    end
    nIND = length(range.family_range);
    if( nIND <= 0 )
        error = 1;
        return;
    end
    [r, c] = size(vit);
    if( r ~= nIND || c ~= nIND )
        error = 1;
        return;
    end
    [r, c] = size(ml);
    if( r ~= nIND || c ~= nIND )
        error = 1;
        return;
    end
    [d1, d2, d3, d4] = size(posIBD1);
    if( d1 ~= nIND || d2 ~= nIND || d3 ~= 2 || d4 ~= 2 )
        error = 1;
        return;
    end
    [d1, d2, d3] = size(pair_pos);
    if( d1 ~= nIND || d2 ~= nIND || d3 ~= 3 )
        error = 1;
        return;
    end
%%
    ambig_triples = triples.ambig_triples;
    num_triples = triples.num_triples;
    triple_viewpair = triples.triple_viewpair;
    [d1, d2, d3] = size(ambig_triples);
    if( d1 ~= nIND || d2 ~= nIND || d3 ~= nIND )
        error = 1;
        disp('error in triplets');
        return;
    end
    [d1, d2, d3] = size(triple_viewpair);
    if( d1 ~= nIND || d2 ~= nIND || d3 ~= nIND )
        error = 1;
        disp('error in triplets');
        return;
    end
    if( max(max(max(ambig_triples))) > num_triples )
        error = 1;
        disp('error in triplets');
        return;
    end
    
    printid(1:nIND) = range.structure(range.family_range,8);
    [clique_config, error] ...
        = GraphPartition(vit, ...
        kinship2ex, ...
        singleLINEAL, ...
        posIBD1, pair_pos, ...
        triples, ...
        range, verify, debug, printid);    
    if( error ~= 0 )
        disp('error in segregating homologous chromosomes');
        return;
    end
    
    [alleles_all error] = clique2alleles(clique_config, range, kinship2ex);
    if( error ~= 0 )
        disp('error in segregating homologous chromosomes');
        return;
    end
end


















