 function [est_bg, est_tran, temp1] = get_estimate2(I_prob_o, T_prob_o, O_prob_o, oblist, states, alpha, beta, sampled_markerlist)


Nst=length(states);

intervals = sampled_markerlist(:,2);

lob=length(oblist);

sum1 = 0;
sum2 = 0;
sum3 = 0;

% [1p2p, 1p2m, 1m2p, 1m2m] included states
states = [
    1, 1, 1, 1; %1 %ibd0
    2, 1, 1, 1; %2 %ibd1
    3, 1, 1, 1; %3
    1, 2, 1, 1; %4
    1, 3, 1, 1; %5
    1, 1, 2, 1; %6
    1, 1, 3, 1; %7
    1, 1, 1, 2; %8
    1, 1, 1, 3; %9
    2, 1, 1, 2; %10 %ibd2
    3, 1, 1, 2; %11
    2, 1, 1, 3; %12
    3, 1, 1, 3; %13
    1, 2, 2, 1; %14
    1, 3, 2, 1; %15
    1, 2, 3, 1; %16
    1, 3, 3, 1; %17
    ];

s_bg = [3,5,7,9,11,12,15,16];
d_bg = [13,17];

temp1(1:17,1:lob) = 0;

for t = 1:lob
    temp = alpha(:,t).*beta(:,t);
    temp1(1:17,t) = temp;
    temp1(1:17,t) = temp1(1:17,t)./sum( temp1(1:17,t) );
    %temp = temp./sum(temp);
    sum1 = sum1 + sum(temp(s_bg));
    sum2 = sum2 + sum(temp(d_bg));
    sum3 = sum3 + sum(temp);
end

est_bg = (sum1+2*sum2)/(2*sum3);

temp = zeros(Nst,Nst);

count = 0;
sum1 = 0;
    
sum0 = 0;
sum2 = 0;

sum3 = 0;
sum4 = 0;
    
for t = 1:lob-1
    for i = 1:Nst
        for j = 1:Nst
            temp(i,j) = alpha(i,t).*O_prob_o{t+1}(oblist(t+1),j).*T_prob_o{t+1}(i,j).*beta(j,t+1);
        end
    end

    %normalize
    total = sum(sum(temp));
    temp = temp./total;
    interval = intervals(t+1)-intervals(t);
    base = exp(-interval*0.01/1000000);
    
    tran11 = T_prob_o{t+1}(2,2)./sum( T_prob_o{t+1}(2,1:3) );
   % tran1b1 = T_prob_o{t+1}(2,2)./sum( T_prob_o{t+1}(2,:) );
    
    A = temp(3,3);
    B =  sum(temp(3,[1,2,3]));
    if( B ~= 0 )
        sum0 = sum0 + A;
        sum2 = sum2 + B;
        sum1 = sum1 + A/B;
        sum3 = sum3 + A*(A/B)/tran11;
        sum4 = sum4 + A*log(A/B)/log(base);
        count = count + 1;
    end
    
    A = temp(11,11);
    B =  sum(temp(11,[8,10,11]));
    if( B ~= 0 )
        sum0 = sum0 + A;
        sum2 = sum2 + B;
        sum1 = sum1 + A/B;
        sum3 = sum3 + A*(A/B)/tran11;
        sum4 = sum4 + A*log(A/B)/log(base);
        count = count + 1;
    end

    A = temp(13,13);
    B =  sum(temp(13,[9,12,13]));
    if( B ~= 0 )
        sum0 = sum0 + A;
        sum2 = sum2 + B;
        sum1 = sum1 + A/B;
        sum3 = sum3 + A*(A/B)/tran11;
        sum4 = sum4 + A*log(A/B)/log(base);
        count = count + 1;
    end
    
end

est_tran = sum1/count;
est_tran = sum0/sum2;

est_tran = sum3/sum0;
est_tran = sum4/sum0;







