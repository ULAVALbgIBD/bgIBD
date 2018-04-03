function [tran_matrix] = transition_background_aux(prob, transition, m, r)


% change accumulated in one column or one row

% assume is strictly ordered

% r is genetic distance
% probability of no recombination exp(-r)
% recombination fraction, odd number of recombination
% odd: [1 - exp(-2r)]/2 ~ r - o(r^2) if r is small
    % 1 - exp(-r) + o(r^2)
% even >= 2: [1 + exp(-2r)]/2 - exp(-r)
% 0: exp(-r)

% recom = 1 - (exp(-r))^m;
recom = 1 - ((1+exp(-2*r))/2)^m; %probability of recombination on m meioses 
non_recom = ((1+exp(-2*r))/2)^m;

num = length(transition(:,1));

tran_matrix(1:num,1:num) = 0;

for i = 1:num
    if( prob(i) == 0 )
        tran_matrix(i,i) = 1;
        continue;
    end
    outbound = 0;
    for j = 1:i
        outbound = outbound + tran_matrix(j,i);
    end
    for j = i+1:num
        if( prob(j) ~= 0 && transition(i,j) ~= 0 )
            tran_matrix(i,j) = recom * (1-outbound);
            tran_matrix(j,i) = recom * prob(i) / prob(j);
        end
    end
    tran_matrix(i,i) = 1 - sum(tran_matrix(i,[1:i-1,i+1:num]));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% bypass situations for big r

if( recom > 1/num )
    for i = 1:num
        if( prob(i) ~= 0 )
            for j = 1:num
                tran_matrix(i,j) = prob(j);
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




