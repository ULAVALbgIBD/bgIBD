
function [corrected_viterbi error] = correct_viterbi(viterbi, kinship2ex, range, map)

error = 0;
corrected_viterbi = [];
global debug_mode;

if( isempty(range) )
    error = 1;
    disp('error in family structures');
    return;
end
genotyped = range.family_range;
nGENO = length(genotyped);
nIND = length(range.pedigree_range_full);
if( nIND <= 0 || nGENO > nIND || any(genotyped > nIND) )
    error = 1;
    disp('error in pedigree structures');
    return;
end
family = range.structure;
[r, c] = size(family);
if( r <= 0 || r ~= nIND || c < 12 )
    error = 1;
    disp('error in family structures');
    return;
end
reverse_pairs = range.reverse_pairs;
[d1, d2] = size(reverse_pairs);
if( d1 <= 0 || d1 ~= d2 || d1 ~= nGENO )
    error = 1;
    disp('error in family structure');
    return;
end
[d1, d2, d3, d4] = size(kinship2ex);
if( d1 <= 0 || d1 ~= nIND || d1 ~= d2 || d3 ~= 2 || d4 ~= 2 )
    error = 1;
    disp('error in kinship');
    return;
end
pairs = range.pairs;
[nPAIRS, c] = size(pairs);
if( nPAIRS <= 0 || c ~= 2 )
    error = 1;
    disp('error in family structures');
    return;
end

nLOC = length(map);
if( nLOC <= 0 )
    error = 1;
    disp('error in marker list');
    return;
end

[len r c] = size(viterbi);
if( len <= 0 || len ~= nPAIRS )
    error = 1;
    disp('error in viterbi decoding');
    return;
end
if( r <= 0 || r ~= nLOC || c ~= 2 )
    error = 1;
    disp('error in viterbi decoding');
    return;
end


time = cputime;

anc(1:nIND,1:nIND) = false;
for i = 1:nIND
    father = family(i,3);
    mother = family(i,4);
    if( father ~= 0 )
        anc(i,father) = true;
        anc(i,1:nIND) = anc(i,1:nIND) | anc(father,1:nIND);
    end
    if( mother ~= 0 )
        anc(i,mother) = true;
        anc(i,1:nIND) = anc(i,1:nIND) | anc(mother,1:nIND);
    end
end

rel(1:nIND,1:nIND) = false;
for i = 1:nIND
    for j = 1:nIND
        if( anc(i,j) || anc(j,i) || any(anc(i,1:nIND) & anc(j,1:nIND)) )
            rel(i,j) = true;
        end
    end
end

% nearest genotyped shield
% trace to nearest genotyped ancestral source, if none, trace to founders
% the shield is flawless
% assuming ancestor always preceding offspring
[ngs, error] = genotyped_ancestor(family);

epairs(1:nIND, 1:nIND) = 0;
for i = 1:nGENO
    id1 = genotyped(i);
    for j = 1:nGENO
        id2 = genotyped(j);
        epairs(id1, id2) = reverse_pairs(i,j);
    end
end


ekin = kinship2ex;
isGENO(1:nIND) = (family(1:nIND, 7) == 1);
for i = 1:nIND
    if( ~isGENO(i) )
        continue;
    end

    for j = 1:nIND
        if( ~isGENO(j) )
            continue;
        end
        if( i == j )
            continue;
        end
        % skip situations where i is an ancestor of j
        if( anc(j, i) )
            continue;
        end
        % skip situation where j is on the shield of i
        if( ngs(i, j) )
            continue;
        end
        shield = ngs(i, 1:nIND);
        if( ~any(shield) )
            continue;
        end
        filter1 = shield & rel(j, 1:nIND);

        % filter1 is all possible relatives lying between i and j
        if( ~any(filter1) )
            continue;
        end 
        
        if( ~all(isGENO(filter1)) )
            continue;
            % filter1 has holes
        end
        
        p0 = epairs(i,j);
        allpairs = epairs(filter1,j);
        if( p0 <= 0 || p0 > nPAIRS || any(allpairs <= 0) || any(allpairs > nPAIRS) )
            error = 1;
            disp('error in family structures');
            return;
        end
        maskbit = viterbi(allpairs, 1:nLOC, 1);
        curbit(1:nLOC) = viterbi(p0, 1:nLOC, 1) > 1;
        suppress(1:nLOC) = curbit & ~any(maskbit, 1);
        viterbi(p0,suppress, 1) = 1;
        
    end
end

corrected_viterbi = viterbi;

if( debug_mode == 1 )
    display(['          ibd correction costs time: ', num2str(cputime - time), ' seconds']); 
end

end










