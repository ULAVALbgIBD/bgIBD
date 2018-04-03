
% an inheritance path is a path
% that starts from an ancestor goes down to 
% two descendants

% not all paths in graph are valid

% starting point i
% currently at j
% family relationship P
% direction upORdown

function traverse(i, j, P, upORdown)

    global freeALLELE;
    nIND = size(freeALLELE,1);

    if( upORdown )
        for k = 1:nIND
            % if k is the parent of j
            % go to ancestor
            if( P(j,k) )
                if( assignALLELE(i,j,k) )
                    % depth first traverse upwards to k
                    traverse(i, k, P, true);
                end
            end
        end
    else
        for k = 1:nIND
            % if k is the child of j
            % go to descendants
            if( P(k,j) )         
                if( assignALLELE(i,j,k) )
                % depth first traverse downwards to k
                    traverse(i, k, P, false);
                end
            end
        end
    end

    % change direction, to down
    if( upORdown )
        traverse(i, j, P, false);
    end



end

function changed = assignALLELE(i,j,k)

    global freeALLELE;
    edge = freeALLELE(i,k,1:2);
    
    
    if( freeALLELE(i,k,1) == 0 )
        if( freeALLELE(i,j,2) == -1 )
            freeALLELE(i,k,1:2) = [j,k];
        else
            freeALLELE(i,k,1:2) = freeALLELE(i,j,1:2);
        end
    else
        if( freeALLELE(i,j,2) == -1 )
            freeALLELE(i,k,1:2) = -1;
        else
            if( any(freeALLELE(i,k,1:2) ~= freeALLELE(i,j,1:2)) )
                freeALLELE(i,k,1:2) = -1;
            end
        end
    end 
    changed = ( any( edge ~= freeALLELE(i,k,1:2) ) );
    
end











