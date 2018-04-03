function code = identity_coding(vec)

    code(1:length(vec)) = 0;
    for i = 1:length(vec)
        for j = 1:i
            if( vec(i) == vec(j) )
                code(i) = j;
                break;
            end
        end
    end
    
    % uniform the coding to alphabet 1:length(vect)

end