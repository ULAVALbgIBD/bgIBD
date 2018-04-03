 function [est_bg, est_tran] = get_estimate(I_prob_o, T_prob_o, O_prob_o, oblist, states, alpha, beta)

% states={'Sunny','Cloudy','Rainy'}; %The 3 "hidden" states
Nst=length(states);
% obs={'Dry','Dryish','Damp','Soggy'}; % the 4 observations.

%Initial Probabilities of states
% I_prob=[0.63 0.17 0.2]';
% 
% % Transition Probabilities of states
% T_prob=[0.5 .25 .25;.375 .125 .375;.125 .675 .375]; % check 5) in readme (Needs to be transposed)
% 
% % Observation prob( Prob of observation given the state)
% O_prob=[0.6 0.2 0.15 0.05; 0.25 0.25 0.25 0.25; 0.05 0.10 0.35 0.50 ]';
% 
% % Here 1 corresponds to Dry, 2 to Dryish, 3 to Damp and 4 to soggy
% oblist=[1 3 4 1 1];% List of observations

lob=length(oblist);

sum1 = 0;
sum2 = 0;
sum3 = 0;

for t = 1:lob
    temp = alpha(:,t).*beta(:,t);
    %temp = temp./sum(temp);
    sum1 = sum1 + sum(temp([3,6,9]));
    sum2 = sum2 + sum(temp([7,8,9]));
    sum3 = sum3 + sum(temp);
end

est_bg = (sum1+sum2)/(2*sum3);

temp = zeros(Nst,Nst);

count = 0;
sum1 = 0;
    
sum0 = 0;
sum2 = 0;

sum3 = 0;
    
for t = 1:lob-1
    for i = 1:Nst
        for j = 1:Nst
            temp(i,j) = alpha(i,t).*O_prob_o{t+1}(oblist(t+1),j).*T_prob_o{t+1}(i,j).*beta(j,t+1);
        end
    end

    %normalize
    total = sum(sum(temp));
    temp = temp./total;
    
    tran11 = T_prob_o{t+1}(2,2)./sum( T_prob_o{t+1}(2,1:3) );
   % tran1b1 = T_prob_o{t+1}(2,2)./sum( T_prob_o{t+1}(2,:) );
    
    A = temp(8,8);
    B =  sum(temp(8,[2,5,8]));
    if( B ~= 0 )
        sum0 = sum0 + A;
        sum2 = sum2 + B;
        sum1 = sum1 + A/B;
        sum3 = sum3 + A*(A/B)/tran11;
        count = count + 1;
    end
    
    A = temp(7,7);
    B = sum(temp(7,[1,4,7]));
    if( B ~= 0 )
        sum0 = sum0 + A;
        sum2 = sum2 + B;
        sum1 = sum1 + A/B;
        sum3 = sum3 + A*(A/B)/tran11;
        count = count + 1;
    end
    
    A = temp(6,6);
    B = sum(temp(6,[4,5,6]));
    if( B ~= 0 )
        sum0 = sum0 + A;
        sum2 = sum2 + B;
        sum1 = sum1 + A/B;
        sum3 = sum3 + A*(A/B)/tran11;
        count = count + 1;
    end
    
    A = temp(3,3);
    B = sum(temp(3,[1,2,3]));
    if( B ~= 0 )
        sum0 = sum0 + A;
        sum2 = sum2 + B;
        sum1 = sum1 + A/B;
        sum3 = sum3 + A*(A/B)/tran11;
        count = count + 1;
    end
    
    A = temp(9,9);
    B = sum(temp(9,[7,8,9]));
    if( B ~= 0 )
        sum0 = sum0 + A;
        sum2 = sum2 + B;
        sum1 = sum1 + A/B;
        sum3 = sum3 + A*(A/B)/tran11;
        count = count + 1;
    end
    
    A = temp(9,9);
    B = sum(temp(9,[3,6,9]));
    if( B ~= 0 )
        sum0 = sum0 + A;
        sum2 = sum2 + B;
        sum1 = sum1 + A/B;
        sum3 = sum3 + A*(A/B)/tran11;
        count = count + 1;
    end   
end

est_tran = sum1/count;
est_tran = sum0/sum2;

est_tran = sum3/sum0;







