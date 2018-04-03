function [ states_prob, tran_list ] = condense_inheritance( inheritance )

count = 0;

num = inheritance(1,1);

[is transition] = hasse_matrix(num);
ip = inheritance(:,2:end);

states_prob(length(is(:,1))) = 0;

tran_list = [];
% in case no transition, d1234 only

for i = 1:length(ip(:,1))
    code = ip(i,:);
    id = code(2:num+1);
    gf = code(num+2:end);
    prob = code(1);
    for j = 1:length(gf)
        prob = prob * 0.5^(gf(j));
    end    
    for j = 1:length(is(:,1))
        success = 1;
        for k = 1:num
            if( is(j,k) ~= id(k) )
                success = 0;
                break;
            end
        end
        if( success == 1 )
            id_number = j;
            break;
        end
    end
    for j = 1:length(transition(id_number,:))
        if( transition(id_number,j) ~= 0 )
            ind = -transition(id_number,j);
            set1 = 0;
            set2 = 0;
            for k = 1:num
                if( is(j,k) == is(id_number,ind) )
                    set1 = set1 + 2^(k-1);
                end
                if( is(j,k) == is(j,ind) )
                    set2 = set2 + 2^(k-1);
                end
            end
            count = count + 1;
            tran_list(count,1) = prob;
            tran_list(count,2:3) = [id_number,j];
            tran_list(count,4) = gf(set1) + gf(set2);
        end
    end
    states_prob(id_number) = states_prob(id_number) + prob;
end


end

