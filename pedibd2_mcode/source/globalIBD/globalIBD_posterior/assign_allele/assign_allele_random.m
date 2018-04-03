function [ assignment error ] = assign_allele_random(map)

error = 0;
assignment = [];

if( isempty(map) )
    return;
end

num_geno = length(map.id_genotyped);
num_all = length(map.id);
paternal_alleles = map.paternal_alleles;
maternal_alleles = map.maternal_alleles;
relevance = map.relevance;

if( num_all <= 0 || num_geno > num_all )
    error = 1;
    return;
end

[r, c] = size(paternal_alleles);
if( r ~= num_all || c ~= 2 )
    error = 1;
    return;
end
[r, c] = size(maternal_alleles);
if( r ~= num_all || c ~= 2 )
    error = 1;
    return;
end
len = length(relevance);
if( len ~= num_all )
    error = 1;
    return;
end

for i = 1:num_all
    if( any( abs(paternal_alleles(i,1:2)) > i ) )
        error = 1;
        return;
    end
    if( any( abs(maternal_alleles(i,1:2)) > i ) )
        error = 1;
        return;
    end
end

assignment = zeros(num_all,2);

for i = 1:num_all
    if( relevance(i) == 1 )
        assignment(i,1) = -i;
        assignment(i,2) = i;
        p = paternal_alleles(i,1);
        m = maternal_alleles(i,2);
        if( p > 0 )
            assignment(i,1) = assignment(p,2);
        end
        if( p == 0 )
            error = 1;
            return;
        end
        if( p < 0 )
            assignment(i,1) = assignment(-p,1);
        end
        if( m > 0 )
            assignment(i,2) = assignment(m,2);
        end
        if( m == 0 )
            error = 1;
            return;
        end
        if( m < 0 )
            assignment(i,2) = assignment(-m,1);
        end
    end
end

end






































