function [alpha error] = get_forward(I_prob_o, T_prob_o, O_prob_o, oblist, states)

error = 0;
alpha = [];


Nst=length(states);
lob=length(oblist);
if( Nst <= 0 )
    error = 1;
    disp('error in hidden states');
    return;
end
if( lob <= 0 )
    error = 1;
    disp('error in hidden states');
    return;
end
[d1 d2 d3] = size(T_prob_o);
if( d1 ~= lob || d2 ~= Nst || d3 ~= Nst )
    error = 1;
    disp('error in transition probability');
    return;
end
[d1 d2 d3] = size(O_prob_o);
if( d1 ~= lob || d2 ~= Nst || d3 <= 0 )
    error = 1;
    disp('error in emission probability');
    return;
end
n_ibs = d3;


[alpha error] = cForward(I_prob_o, T_prob_o, O_prob_o, oblist, states);
alpha = alpha';
return;


alpha = zeros(Nst, lob);

for t = 1:1
    alpha(1:Nst,t) = I_prob_o(1:Nst).*(O_prob_o(t,1:Nst,oblist(t)));
    alpha(1:Nst,t) = alpha(1:Nst,t)./sum(alpha(1:Nst,t));
end

% transition probability at locus 1 is not used

T = [];
O = [];
pre = [];
cur = [];
for t = 2:lob
    T(1:Nst, 1:Nst) = T_prob_o(t, 1:Nst, 1:Nst);
    O(1:Nst) = O_prob_o(t, 1:Nst, oblist(t));
    pre(1:Nst) = alpha(1:Nst, t-1);
    cur(1:Nst) = (pre*T).*O;
    alpha(1:Nst,t) = cur(1:Nst)./sum(cur(1:Nst));
end


end



