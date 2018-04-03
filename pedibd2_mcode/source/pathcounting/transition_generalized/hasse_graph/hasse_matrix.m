function [ is transition ] = hasse_matrix( num )

is = identity_states(num);
hasse = hasse_diagram(num);

transition(1:length(is),1:length(is)) = 0;

pro = hasse(1:end,1:num);
change = hasse(1:end,num+1);
suc = hasse(1:end,num+2:num+num+1);

for i = 1:length(hasse(:,1))
    for j = 1:length(is)
        success = 1;
        for k = 1:num
            if( is(j,k) ~= pro(i,k) )
                success = 0;
                break;
            end
        end
        if( success == 1 )
            a = j;
            break;
        end
    end
    for j = 1:length(is)
        success = 1;
        for k = 1:num
            if( is(j,k) ~= suc(i,k) )
                success = 0;
                break;
            end
        end
        if( success == 1 )
            b = j;
            break;
        end
    end  
    transition(a,b) = change(i,1);
end

end

