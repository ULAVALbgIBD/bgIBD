
function [raw_score error] = additive1family(family, alleles)

    global fac;
    
    error = 0;
    raw_score = 0;
    
    [rows, cols] = size(family);
    if( rows <= 0 || cols < 12 )
        error = 1;
        disp('error in family structure');
        return;
    end
    
    [nIND, ncol] = size(alleles);
    if( nIND ~= rows || ncol ~= 2 )
        error = 1;
        disp('error in global IBD');
        return;
    end
    
    affected(1:nIND) = 0;
    for i = 1:nIND
        if( family(i,6) == 2 )
            affected(i) = 1;
        end
    end
    
    genotyped(1:nIND) = 0;
    for i = 1:nIND
        if( family(i,7) == 1 )
            genotyped(i) = 1;
        end
    end
    
    % change type
    % otherwise all following operations overflow at 127
    alleles = int32(alleles);
    
    min_allele = min(min(alleles));
    max_allele = max(max(alleles));
    range = max_allele - min_allele + 1;
    offset = - min_allele + 1;
    
    if( isempty(fac) || length(fac) < range )
        for i = 1:range
            fac(i) = factorial(double(i));
        end
    end
    
    allele_count = zeros(range,1);
    for i = 1:nIND
        % only include affected, assigned individuals
        if( affected(i) == 1 && genotyped(i) == 1 )
            p = alleles(i,1);
            m = alleles(i,2);
            if( p ~= 0 )
                p = p + offset;
                allele_count(p) = allele_count(p) + 1;
            end
            if( m ~= 0 )
                m = m + offset;
                allele_count(m) = allele_count(m) + 1;
            end          
        end
    end
    
    raw_score = 1;
    for i = 1:range
        if( allele_count(i) > 0 )
            raw_score = raw_score * fac(allele_count(i));
        end
    end
    
    % if no genotyped affected individuals
    % raw_score will be 1
    
end
















