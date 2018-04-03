
function [imputed_alleles error] = imputation(alleles_all, range, kinship2ex)

error = 0;
imputed_alleles = [];
global debug_mode;

%%

if( isempty(range) )
    error = 1;
    disp('error in family structures');
    return;
end
genotyped = range.family_range;
nGENO = length(genotyped);
nIND = length(range.pedigree_range_full);
if( nIND <= 0 )
    error = 1;
    disp('error in family structures');
    return;
end
if( nGENO > nIND || any(genotyped > nIND) )
    error = 1;
    disp('error in family structures');
    return;
end
family = range.structure;
[r c] = size(family);
if( r ~= nIND || c < 12 )
    error = 1;
    disp('error in family structures');
    return;
end
[r c] = size(alleles_all);
if( r ~= nIND || c ~= 2 )
    error = 1;
    disp('error in homologous chromosome segregation');
    return;
end
if( nGENO <= 1 )
    imputed_alleles = alleles_all;
    return;
end
[d1 d2 d3 d4] = size(kinship2ex);
if( d1 <= 0 || d1 ~= nIND || d1 ~= d2 || d3 ~= 2 || d4 ~= 2 )
    error = 1;
    disp('error in kinship');
    return;
end


isGENO(1:nIND) = (family(1:nIND,7) == 1);
for i = 1:nIND
    if( isGENO(i) )
        if( any(alleles_all(i,1:2) == 0) )
            error = 1;
            disp('genotyped individual not assigned');
            return;
        end
    end
end

%%

[alleles_all phased error] = phase_alleles(alleles_all, family, kinship2ex);
if( error == 1 )
    disp('error in determining parental sources');
    return;
end

% children voting for each allele
unGENOTYPED(1:nIND,1:2*nIND+1) = 0;
allele_map(1:2*nIND+1) = 0;
allele_map(1:nIND) = -nIND:-1;
allele_map(nIND+2:2*nIND+1) = 1:nIND;
for i = 1:nIND
    if( all(alleles_all(i,1:2) ~= 0, 2) )
        if( ~phased(i) )
            continue;
        end
        p = alleles_all(i,1);
        m = alleles_all(i,2);
        if( p == 0 || m == 0 )
            error = 1;
            disp('error in IBD partition');
            return;
        end
        if( p > nIND || p < -nIND || m > nIND || m < -nIND )
            error = 1;
            disp('error in IBD partition');
            return;
        end
        
        father = family(i,3);
        mother = family(i,4);
        if( father ~= 0 && ~isGENO(father) )
            unGENOTYPED(father, p + nIND + 1) = unGENOTYPED(father, p + nIND + 1) + 1;
        end
        if( mother ~= 0 && ~isGENO(mother) )
            unGENOTYPED(mother, m + nIND + 1) = unGENOTYPED(mother, m + nIND + 1) + 1;
        end
    end
end


imputed_alleles = alleles_all;
for i = 1:nIND
    oneline = unGENOTYPED(i,1:2*nIND+1);
    if( any(oneline > 0) )
        if( any(imputed_alleles(i,1:2) > 0) )
            % reassign in new iteration
            imputed_alleles(i,1:2) = 0;
        end
        if( nnz(oneline) == 2 )
            imputed_alleles(i,1:2) = allele_map(oneline > 0);
        end
        if( nnz(oneline) == 1 )
            imputed_alleles(i,1) = allele_map(oneline > 0);
        end
        if( nnz(oneline) > 2 )
            if( debug_mode == 1 )
                disp(['           ', num2str(i), ': conflict in imputing alleles: more than 2 alleles']);
            end
            [~, index] = sort(oneline, 'descend');
            imputed_alleles(i,1:2) = allele_map(index(1:2));
        end
    end
end

imputed_alleles = remapFOUNDER(imputed_alleles);

end



















