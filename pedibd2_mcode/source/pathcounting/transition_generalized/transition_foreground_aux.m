function [tran_matrix] = transition_foreground_aux(margin_prob, tran_list, r)


% assume is strictly ordered

% r is GENETIC DISTANCE
% probability of no recombination exp(-r)
% recombination fraction, odd number of recombination
% odd: [1 - exp(-2r)]/2 ~ r - o(r^2) if r is small
    % 1 - exp(-r) + o(r^2)
% even >= 2: [1 + exp(-2r)]/2 - exp(-r)
% 0: exp(-r)

% recom = 1 - (exp(-r))^m;
no_recom = ((1+exp(-2*r))/2); % probability of no recombination on 1 meiosis 

max_len = 0;

tran_matrix(1:length(margin_prob),1:length(margin_prob)) = 0;

if( ~isempty(tran_list) )
    for i = 1:length(tran_list(:,1))
        p = tran_list(i,1);
        pro = tran_list(i,2);
        suc = tran_list(i,3);
        len = tran_list(i,4);
        tran_matrix(pro,suc) = tran_matrix(pro,suc) + p * (1-no_recom^len);
        if( max_len < len ) 
            max_len = len;
        end
    end
end

for i = 1:length(tran_matrix(:,1))
    if( margin_prob(i) ~= 0 )
        tran_matrix(i,:) = tran_matrix(i,:)./margin_prob(i); %normalize
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% modify to bypass situations when r is large

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% reverse transition probability, assuming equilibrium
for i = 1:length(tran_matrix(:,1))
    for j = i+1:length(tran_matrix(:,1))
        if( tran_matrix(i,j) ~= 0 )
            tran_matrix(j,i) = ( tran_matrix(i,j) * margin_prob(i) ) / margin_prob(j);
        end
    end
end


for i = 1:length(tran_matrix(:,1))
    temp = sum(tran_matrix(i,[1:i-1,i+1:end]));
    tran_matrix(i,i) = 1 - temp;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% bypass situations for big r
num = length(margin_prob);

if( 1-no_recom^max_len > 1/num )
    for i = 1:num
        if( margin_prob(i) ~= 0 )
            for j = 1:num
                tran_matrix(i,j) = margin_prob(j);
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end










