function [assignment, typing_error, error] = assign_locus_vote(constraint, genotype)

assignment = [];
error = 0;
typing_error = 0;

[num, cols] = size(constraint);
if( num <= 0 || cols ~= 2 )
    error = 1;
    return;
end

allele_vote(1:num,1:2,1:2) = 0;
assignment(1:num,1:2) = 0;


for i = 1:num
    for j = 1:2
        source = constraint(i,j);
        if( source > 0 )
            if( source > num )
                error = 1;
                return;
            end
            if( constraint(source,2) ~= source )
                error = 1;
                return;
            end
        end
        if( source < 0 )
            if( source < -num )
                error = 1;
                return;
            end
            if( constraint(-source,1) ~= source )
                error = 1;
                return;
            end
        end
    end
end

for i = 1:num
    p = genotype(i,1);
    m = genotype(i,2);
    if( p < 0 || m < 0 || p > 2 || m > 2 )
        disp('error in genotypes, not bi-allelic');
    end    
    if( p == 0 || m == 0 )
        continue;
    end    
    for j = 1:2
        anc = constraint(i,j);
        if( anc == 0 )
            continue;
        end
        if( anc > 0 )
            allele_vote(anc,2,p) = allele_vote(anc,2,p) + 1;
            allele_vote(anc,2,m) = allele_vote(anc,2,m) + 1;
        end
        if( anc < 0 )
            allele_vote(-anc,1,p) = allele_vote(-anc,1,p) + 1;
            allele_vote(-anc,1,m) = allele_vote(-anc,1,m) + 1;
        end
    end
end

for i = 1:num
    max_vote(1:2) = 0;
    vote(1:2,1:2) = 0;
    locgeno(1:2) = genotype(i,1:2);
    for j = 1:2
        anc = constraint(i,j);
        if( anc == 0 )
            continue;
        end
        if( anc > 0 )
            vote(j,1:2) = allele_vote(anc,2,1:2);
            if( allele_vote(anc,2,1) > allele_vote(anc,2,2) )
                max_vote(j) = 1;
            end
            if( allele_vote(anc,2,2) > allele_vote(anc,2,1) )
                max_vote(j) = 2;
            end
        end
        if( anc < 0 )
            vote(j,1:2) = allele_vote(-anc,1,1:2);
            if( allele_vote(-anc,1,1) > allele_vote(-anc,1,2) )
                max_vote(j) = 1;
            end
            if( allele_vote(-anc,1,2) > allele_vote(-anc,1,1) )
                max_vote(j) = 2;
            end
        end
    end
    if( ~all(locgeno(1:2)) )
        assignment(i,1:2) = max_vote(1:2);
    end
    if( all(locgeno(1:2)) )
        if( vote(1, locgeno(1)) + vote(2, locgeno(2)) >= vote(1, locgeno(2)) + vote(2, locgeno(1)) )
            assignment(i,1:2) = locgeno(1:2);
        else
            assignment(i,1:2) = locgeno(2:-1:1);
        end
    end
    if( all(locgeno(1:2)) )
        if( all(max_vote(1:2)) )
            if( (max_vote(1) == locgeno(1) && max_vote(2) == locgeno(2)) )
            else
                if( max_vote(1) == locgeno(2) && max_vote(2) == locgeno(1) )
                else
                    typing_error = typing_error + 1;
                end
            end
        end
        if( nnz(max_vote(1:2)) == 1 )
            if( max_vote(1) ~= 0 && ~any(locgeno(1:2) == max_vote(1)) )
                typing_error = typing_error + 1;
            end
            if( max_vote(2) ~= 0 && ~any(locgeno(1:2) == max_vote(2)) )
                typing_error = typing_error + 1;
            end
        end
    end
end


end



















