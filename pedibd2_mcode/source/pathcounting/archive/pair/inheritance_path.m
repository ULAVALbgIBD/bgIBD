
%make sure parents are relative indexed by id

function vector = inheritance_path( ind1, ind2, parents, ancestor, max_path )
    if( ind1 ~= 0 && ind2 ~= 0 )
        f1 = parents(ind1,2);
        m1 = parents(ind1,3);
        f2 = parents(ind2,2);
        m2 = parents(ind2,3);
        if( ind1 == ind2 )
            vector(1:max_path) = 0;
            vector(1) = 2;
            if( f1 ~= 0 && m1 ~= 0 )
                temp = inheritance_path(f1, m1, parents, ancestor, max_path);
                vector(1+2:max_path) = vector(1+2:max_path) + temp(1:max_path-2);
            end
        else
            if( ancestor(ind1,ind2) )
                temp = inheritance_path(f1, ind2, parents, ancestor, max_path) + inheritance_path(m1, ind2, parents, ancestor, max_path);
            else
                temp = inheritance_path(f2, ind1, parents, ancestor, max_path) + inheritance_path(m2, ind1, parents, ancestor, max_path);
            end
            vector(1:max_path) = 0;
            vector(2:max_path) = vector(2:max_path) + temp(1:max_path-1);        
        end
    else
        vector(1:max_path) = 0;
    end
end

