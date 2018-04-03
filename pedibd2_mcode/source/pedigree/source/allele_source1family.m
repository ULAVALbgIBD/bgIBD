function [relevance, error] = allele_source1family(family)

    relevance = [];
    error = 0;

    if( isempty(family) )
        error = 1;
        disp('empty family');
        return;
    end

    [rows, cols] = size(family);
    if( rows <= 0 || cols < 12 )
        error = 1;
        disp('error in family structure');
        return;
    end
    
    nind = rows;
    
    ids = unique(family(:,2:4));
    if( max(ids) > nind || min(ids) < 0 )
        error = 1;
        disp('error in family structure');
        return;
    end
    
    % nearest genotyped ancestor
    nga(1:nind,1:nind) = 0;
    % founder ancestor
    fa(1:nind,1:nind) = 0;
    % all ancestor
    aa(1:nind,1:nind) = 0;
    
    genotyped(1:nind) = 0;
    founder(1:nind) = 0;
    for i = 1:nind
        id = family(i,2);
        father = family(i,3);
        mother = family(i,4);
        if( id ~= i )
            error = 1;
            disp('error in family indexing');
            return;
        end
        if( father ~= 0 )
            if( father >= id )
                error = 1;
                disp('error in family indexing');
                return;
            end
        end
        if( mother ~= 0 )
            if( mother >= id )
                error = 1;
                disp('error in family indexing');
                return;
            end
        end
        if( father == 0 && mother == 0 )
            founder(id) = 1;
        end
        if( family(i,7) == 1 )
            genotyped(id) = 1;
        end
    end
    
    for i = 1:nind
        id = family(i,2);
        father = family(i,3);
        mother = family(i,4);
        if( father ~= 0 )
            aa(id,father) = 1;
            aa(id,:) = aa(id,:) | aa(father,:);
        end
        if( mother ~= 0 )
            aa(id,mother) = 1;
            aa(id,:) = aa(id,:) | aa(mother,:);
        end
    end
    
    % included in the transmission enumeration
    relevance(1:nind) = 0;
    for i = 1:nind
        if( genotyped(i) == 1 )
            relevance(i) = 1;
            relevance(aa(i,:)>0) = 1;
        end
    end
    
    
    
end