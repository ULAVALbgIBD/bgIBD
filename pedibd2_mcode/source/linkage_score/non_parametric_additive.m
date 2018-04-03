

function [z_score, p_value, error] = non_parametric_additive(alleles_all, family)

global debug_mode;
error = 0;
z_score = [];
p_value = [];


[nIND, fields] = size(family);
if( nIND <= 0 || fields < 12 )
    error = 1;
    disp('error in family structure');
    return;
end
genotyped = (family(:,7) == 1);


if( ndims(alleles_all) ~= 3 )
    error = 1;
    disp('error in global IBD');
    return;
else
    [nSEG, d2, d3] = size(alleles_all);
    if( nSEG <= 0 || d2 ~= nIND || d3 ~= 2 )
        error = 1;
        disp('error in global IBD');
        return;
    end
end


% linkage score only consider genotyped individuals
raw_score(1:nSEG) = 0;
for i = 1:nSEG
    alleles = reshape(alleles_all(i,1:nIND,1:2), [nIND,2]);
    [r,c] = size(alleles);
    if( r ~= nIND || c ~= 2 )
        error = 1;
        disp('error in global IBD');
        return;
    end
    if( any( alleles(genotyped,1:2) == 0 ) )
        error = 1;
        disp('error in global IBD, not all genotyped individuals assigned');
        return;
    end
    raw_score(i) = additive1family(family, alleles);
end

% p value is based on permutation
perm_times = min(100000, 10 * 2^(2*nIND));
perm_score(1:perm_times) = 0;
for i = 1:perm_times
    perm_score(i) = additive1family_randomized(family);
end

% empirical p value
p_value = ones(nSEG,1);
for i = 1:nSEG
    gt = nnz(perm_score > raw_score(i));
    eq = nnz(perm_score == raw_score(i));
    lt = nnz(perm_score < raw_score(i));
    p_value(i) = (gt + eq./2)./perm_times;
end
p_value(p_value <= 0) = 1./(2 * perm_times);
p_value(p_value >= 1) = 1 - 1./(2 * perm_times);


mu = mean(perm_score);
delta = std(perm_score);

z_score = zeros(nSEG, 1);
if( delta == 0 )
    for i = 1:nSEG
        if( raw_score(i) ~= mu )
            if( debug_mode )
                error = 1;
                disp('error in computing linkage score, no variance detected');
                return;
            end
        end
        z_score(i) = 0;
    end
else
    for i = 1:nSEG
        z_score(i) = norminv(1 - p_value(i), 0, 1);
    end
end




end







