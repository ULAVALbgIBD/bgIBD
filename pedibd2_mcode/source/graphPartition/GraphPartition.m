
function [clique_config error] ...
    = GraphPartition(pairs, ...
    kinship2ex, ...
    singleLINEAL, ...
    posteriorIBD1, posterior, ...
    triples, ...
    range, verify, debug, printid)

%%

error = 0;
clique_config = [];

[nGENO, c] = size(pairs);
if( nGENO <= 0 || nGENO ~= c )
    error = 1;
    disp('error in input pairs');
    return;
end
nIND = length(range.pedigree_range_full);
if( nIND <= 0 )
    error = 1;
    disp('error in family structures');
    return;
end
genotyped = range.family_range;
if( length(genotyped) ~= nGENO || any(genotyped > nIND) )
    error = 1;
    disp('error in family structures');
    return;
end
family = range.structure;
[r, c] = size(family);
if( r <= 0 || r ~= nIND || c < 12 )
    error = 1;
    disp('error in family structures');
    return;
end
kinship2 = kinship2ex(genotyped, genotyped, 1:2, 1:2);

%%
tolerance = 0.5;
spectrum(1:nGENO,1:nGENO) = 0;
for i = 1:nGENO

    index = i;
    [final_p_clique, final_m_clique, final_pm_clique, error] ...
        = segregate_homologous(index, pairs, triples, ...
        kinship2, ...
        singleLINEAL(genotyped, genotyped, 1:2), ...
        posterior, posteriorIBD1, tolerance);
    
    if( error ~= 0 )
        disp('error in chromosome segregating');
        return;
    end
    
    spectrum(index * 2 - 1, [final_p_clique, final_pm_clique]) = index * 2 - 1;
    spectrum(index * 2, [final_m_clique, final_pm_clique]) = index * 2;

end


% due to consolidation, not always founder paternal allele is assigned

[T error] = consolidate_homology(spectrum, range, kinship2, debug, printid);
if( error ~= 0 )
    disp('error in forming chromosome spectrum');
    return;
end

% get the haplotype assign only at itself partition i*2-1, i*2, lines
% if this line is not merged with other lines, will be left with 
% its own haplotype coding
chr_dist(1:2*nGENO,1:nGENO) = 0;
for i = 1:nGENO
    p = T(2*i-1);
    m = T(2*i);
    if( p == m )
        disp(['          homologous chromosomes fused, at individual ', num2str(printid(i))]);
        if( p == i*2 - 1 )
            m = i*2;
        else
            p = i*2 - 1;
        end
    end
    chr_dist(p,i) = 1;
    chr_dist(m,i) = 1;
end


clique_config.rep_list = [1:2*nGENO];
clique_config.num_cliques = 2*nGENO;
clique_config.assignment = chr_dist;


if( debug )
    figure;
    if( isempty(verify) )
        display_cliques(chr_dist, pairs, printid);
    else
        subplot(1,2,1);        
        title('graph to viterbi');
        display_cliques(chr_dist, pairs, printid);
        subplot(1,2,2);
        title('graph to likelihood');
        display_cliques(chr_dist, allele2pair(verify), printid);
    end
end

end

















