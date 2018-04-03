function [id_prob, tran_matrix_all] = transition_background(num, markerlist)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[is transition] = hasse_matrix(num);

mg = 0.0113; % background marginal probability, LD level
m = 200;     % number of meioses apart

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

unit = 10^6; % base pair unit
cM = 0.01;

for i = 1:length(markerlist)    
    if( i == 1 )
        r(i) = cM*(markerlist(i,2))/unit;
    else
        r(i) = cM*(markerlist(i,2)-markerlist(i-1,2))/unit;
    end 
end

r(1) = mean(r(2:end));

num_is = length(is(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

count = 0;
for i = 1:num
    for j = i+1:num
        count = count + 1;
        map(count, 1:2) = [i,j];
    end
end

state(1:2^count,1:num) = 0;
plain(1:2^count,1:count) = 0;

for i = 1:2^count
    for j = 1:num
        state(i,j) = j;
    end
    prob = 1;
    for k = 1:count
        if( bitget(i-1,k) == 1 )
            state(i,map(k,2)) = state(i,map(k,1));
            plain(i,k) = 2;
            prob = prob * mg;
        else
            plain(i,k) = 1;
            prob = prob * (1-mg);
        end
    end
    state_prob(i) = prob;
end


for i = 1:length(is(:,1))
    id2state{i} = [];
    id_prob(i) = 0;
end

for i = 1:2^count
    for j = 1:num_is
        suc = 1;
        for k = 1:num
            if( state(i,k) ~= is(j,k) )
                suc = 0;
                break;
            end
        end
        if( suc == 1 )
            state2id(i) = j;
            id2state{j} = union(id2state{j}, i);
            id_prob(j) = id_prob(j) + state_prob(i);
            break;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for l = 1:length(r)
    
    recom = 1 - ((1+exp(-2*r(l)))/2)^m; %probability of recombination on m meioses 
    non_recom = ((1+exp(-2*r(l)))/2)^m;
    
    t(2,2) = non_recom;
    t(2,1) = recom;
    t(1,2) = mg*recom/(1-mg);
    t(1,1) = 1-t(1,2);

    tran_matrix(1:num_is,1:num_is) = 0;
    for i = 1:num_is
        for j = 1:num_is
            for k_i = 1:length(id2state{i})
                code1 = id2state{i}(k_i);
                for k_j = 1:length(id2state{j})
                    code2 = id2state{j}(k_j);
                    temp = 1;
                    for s = 1:count
                        temp = temp * t(plain(code1,s),plain(code2,s));
                    end
                    tran_matrix(i,j) = tran_matrix(i,j) + temp * state_prob(code1);
                end
            end
            tran_matrix(i,j) = tran_matrix(i,j)/id_prob(i);
        end
    end
    
    tran_matrix_all{l} = tran_matrix;
    
end

end
