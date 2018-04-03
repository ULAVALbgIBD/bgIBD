
function states = identity_states(num)

% all subset combinations of size of num
% each set can choose to join a previous set, or emerge as a new one
% bell number

% the sequence is strictly partial ordered

states = [];
current(1:num) = 0;
i = 1;
j = 1;
while(1)
    if( j == num + 1 )
        states(i,1:num) = current;
        i = i + 1;
        j = j - 1;
        continue;
    end
    if( j == 0 )
        break;
    end
    
    forward = 0;
    for k = 1:j-1
        if( current(k) > current(j) )
            % join each of the previous sets
            current(j) = current(k);
            j = j + 1;
            forward = 1;
            break;
        end
    end
    if( forward == 0 && j > current(j) )
        current(j) = j;
        j = j + 1;
        forward = 1;
    end
    
    if( forward == 1 )
        continue;
    else
        current(j) = 0;
        j = j - 1;
    end
end


