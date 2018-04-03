function [ consistency ] = set_assign( index, assign, offset )

% assign is either 1 or 2

global ancestral;

if( index < 0 )
    if( ancestral(-index,5) == 0 )
        if( offset == 0 )
            if( assign == 1 )
                ancestral(-index,5) = 1;
            else
                ancestral(-index,5) = 2;
            end
        else
            if( assign == 1 )
                ancestral(-index,5) = 2;
            else
                ancestral(-index,5) = 1;
            end
        end
        consistency = 1;
    else
        if( offset == 0 )
            if( ancestral(-index,5) ~= assign )
                consistency = 0;
            else
                consistency = 1;
            end
        else
            if( ancestral(-index,5) == assign )
                consistency = 0;
            else
                consistency = 1;
            end
        end
    end
else
    % index > 0
    if( ancestral(index,6) == 0 )
        if( offset == 0 )
            if( assign == 1 )
                ancestral(index,6) = 1;
            else
                ancestral(index,6) = 2;
            end
        else
            if( assign == 1 )
                ancestral(index,6) = 2;
            else
                ancestral(index,6) = 1;
            end
        end
        consistency = 1;
    else
        if( offset == 0 )
            if( ancestral(index,6) ~= assign )
                consistency = 0;
            else
                consistency = 1;
            end
        else
            if( ancestral(index,6) == assign )
                consistency = 0;
            else
                consistency = 1;
            end
        end
    end
end

end

