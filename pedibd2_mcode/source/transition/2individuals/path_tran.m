function [sum_prio, sum_tran] = path_tran(pt, r)

sum_prio = 0;
sum_tran = 0;

% pt is the number of inheritance path of length i

for i = 1:length(pt)
    num_path = pt(i);
    path_length = i;
    prio = num_path*((1/2)^path_length);
%     tran = (1-r)^path;
    tran = (exp(-r))^path_length;
    tran = ((1+exp(-2*r))/2)^path_length;
    
    sum_tran = sum_tran + prio*tran;
    sum_prio = sum_prio + prio;
end