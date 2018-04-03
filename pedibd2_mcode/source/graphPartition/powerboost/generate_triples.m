function [triples, error] = generate_triples(family_range, kinship, kinship2)


triples = [];
error = 0;

global debug_mode;

if( isempty(family_range) )
    error = 1;
    return;
end

nIND = length(family_range.pedigree_range_full);
if( nIND <= 0 )
    error = 1;
    disp('error in family structures');
    return;
end

family = family_range.structure;
[r, c] = size(family);
if( r <= 0 || r ~= nIND || c < 12 )
    error = 1;
    disp('error in family structures');
    return;
end

genotyped = family_range.family_range;
nGENO = length(genotyped);
if( nGENO < 1 )
    triples.ambig_triples = [];
    triples.deter_triples = [];
    triples.total_triples = 0;
    return;
end
if( nGENO < 3 && nGENO > 0 )
    temp(1:nGENO,1:nGENO,1:nGENO) = 0;
    triples.ambig_triples = temp;
    triples.deter_triples = temp;
    triples.total_triples = 0;
    return;
end
if( any(genotyped > nIND) || any(genotyped < 0 ) )
    error = 1;
    disp('error in family structures');
    return;
end

[r, c] = size(kinship);
if( r ~= nGENO || c ~= nGENO )
    error = 1;
    disp('error in kinship');
    return;
end
[d1, d2, d3, d4] = size(kinship2);
if( d1 ~= nGENO || d2 ~= nGENO || d3 ~= 2 || d4 ~= 2 )
    error = 1;
    disp('error in kinship');
    return;
end

% from ip and im to j and k
% there must be independent paths
% there must be a path from both j, k to either ip or im
% there must be paths linking j and k on the other alleles
% i, j, k

time = cputime;

isvalid = false(nGENO, nGENO, nGENO);
determined = zeros(nGENO, nGENO, nGENO);
for i = 1:nGENO
    for j = 1:nGENO
        for k = 1:nGENO
            if( i == j || j == k || i == k )
                continue;
            end
            if( kinship(i,j) == 0 || kinship(j,k) == 0 || kinship(i,k) == 0 )
                continue;
            end
            sit1 = false;
            sit2 = false;
            for p = 1:2
                for q = 1:2
                    for r = 1:2
                        if( kinship2(i,j,p,q) ~= 0 && kinship2(j,k,q,r) ~= 0 && kinship2(k,i,r,p) ~= 0 )
                            sit1 = true;
                        end
                        p2 = mod(p,2) + 1;
                        q2 = mod(q,2) + 1;
                        r2 = mod(r,2) + 1;
                        if( kinship2(i,j,p,q2) ~= 0 && kinship2(j,k,q,r2) ~= 0 && kinship2(k,i,r,p2) ~= 0 )
                            sit2 = true;
                        end
                    end
                end
            end
            if( sit1 && sit2 )
                isvalid(i,j,k) = true;
            end
            if( sit1 && ~sit2 )
                determined(i,j,k) = 1;
            end
            if( ~sit1 && sit2 )
                determined(i,j,k) = -1;
            end
        end
    end
end

for i = 1:nGENO
    for j = 1:nGENO
        for k = 1:nGENO
            if( isvalid(i,j,k) ~= isvalid(i,k,j) )
                error = 1;
                return;
            end
            if( isvalid(j,i,k) ~= isvalid(j,k,i) || isvalid(i,j,k) ~= isvalid(j,i,k) )
                error = 1;
                return;
            end
            if( isvalid(k,i,j) ~= isvalid(k,j,i) || isvalid(k,i,j) ~= isvalid(i,j,k) )
                error = 1;
                return;
            end
        end
    end
end

count = 0;
valid_triples = zeros(nGENO, nGENO, nGENO);
for i = 1:nGENO
    for j = i+1:nGENO
        for k = j+1:nGENO
            if( isvalid(i,j,k) )
                count = count + 1;
                valid_triples(i,j,k) = count;
                valid_triples(i,k,j) = count;
                valid_triples(k,i,j) = count;
                valid_triples(k,j,i) = count;
                valid_triples(j,k,i) = count;
                valid_triples(j,i,k) = count;
            end
        end
    end
end


triples.ambig_triples = valid_triples;
triples.total_triples = count;
triples.deter_triples = determined;

if( debug_mode == 1 )
    display(['          triple generating costs time: ', num2str(cputime - time), ' seconds']);  
end

end






















































