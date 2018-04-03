
%be careful of two options below, additive and nonadditive

function [raw_score, error] = additive1family_randomized(family)

error = 0;
raw_score = 0;

[rows, cols] = size(family);

if( rows <= 0 || cols < 12 )
    error = 1;
    disp('error in family structure');
    return;
end

nind = rows;

for i = 1:nind
    id = family(i,2);
    father = family(i,3);
    mother = family(i,4);
    if( id ~= i || father >= id || mother >= id )
        error = 1;
        disp('error in family structure');
        return;
    end
end

alleles(1:nind,1:2) = 0;
for i = 1:nind
    id = family(i,2);
    father = family(i,3);
    mother = family(i,4);
    if( father == 0 )
        alleles(id,1) = -i;
    else
        if( rand() > 0.5 )
            alleles(id,1) = alleles(father,1);
        else
            alleles(id,1) = alleles(father,2);
        end
    end
    if( mother == 0 )
        alleles(id,2) = i;
    else
        if( rand() > 0.5 )
            alleles(id,2) = alleles(mother,1);
        else
            alleles(id,2) = alleles(mother,2);
        end
    end
end

[raw_score, error] = additive1family(family, alleles);


end


