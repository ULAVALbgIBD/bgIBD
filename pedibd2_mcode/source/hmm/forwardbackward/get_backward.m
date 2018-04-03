function [beta error] = get_backward(I_prob_o, T_prob_o, O_prob_o, oblist, states)

beta = [];
error = 0;

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

[beta error] = cBackward(I_prob_o, T_prob_o, O_prob_o, oblist, states);
beta = beta';

return;

beta = zeros(Nst, lob);

for t = lob:lob
    beta(1:Nst,t) = 1;
    beta(1:Nst,t) = beta(1:Nst,t)./sum(beta(1:Nst,t));
end

% transition probability at locus 1 is not valid

T = [];
O = [];
pos = [];
cur = [];
for t = lob-1:-1:1
    T(1:Nst,1:Nst) = T_prob_o(t+1, 1:Nst, 1:Nst);
    O(1:Nst) = O_prob_o(t+1, 1:Nst, oblist(t+1));
    pos(1:Nst) = beta(1:Nst,t+1);
    cur(1:Nst) = T * (pos .* O)';
    beta(1:Nst,t) = cur(1:Nst)./sum(cur(1:Nst));   
end



end







