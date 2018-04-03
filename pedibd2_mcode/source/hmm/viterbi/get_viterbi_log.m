function [dec_state, del, error] = get_viterbi_log(I_prob, T_prob_l, O_prob_l, oblist, states)

dec_state = [];
del = [];
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
[d1 d2 d3] = size(T_prob_l);
if( d1 ~= lob || d2 ~= Nst || d3 ~= Nst )
    error = 1;
    disp('error in transition probability');
    return;
end
[d1 d2 d3] = size(O_prob_l);
if( d1 ~= lob || d2 ~= Nst || d3 <= 0 )
    error = 1;
    disp('error in emission probability');
    return;
end
n_ibs = d3;


del = zeros(Nst, lob);
maxlist = del;

I_prob_l = I_prob;
I_prob_l(I_prob<=0) = 10e-10;
I_prob_l = log(I_prob_l);

[dec_state, del, pre, error] = cViterbi(I_prob_l, T_prob_l, O_prob_l, oblist);

return;



del = zeros(Nst, lob);
maxlist = del;

O = [];
T = [];
stuffing(1:Nst) = 1;

for t = 1:lob
    O(1:Nst) = O_prob_l(t, 1:Nst, oblist(t));
    T(1:Nst, 1:Nst) = T_prob_l(t, 1:Nst, 1:Nst);
    if( t == 1 )
        del(1:Nst,t) = I_prob_l(1:Nst) + O(1:Nst); %Initialization
        continue;
    end
    pre(1:Nst, 1:Nst) = del(1:Nst, t-1)*stuffing + T(1:Nst, 1:Nst);
    [del(1:Nst,t) maxlist(1:Nst, t)] = max(pre(1:Nst, 1:Nst));
    del(1:Nst,t)=del(1:Nst,t) + O(1:Nst)';

end


dec_state(1:lob) = 0;
[pstar dec_state(lob)]=max(del(1:Nst,lob));

for t=lob-1:-1:1
    dec_state(t)=maxlist(dec_state(t+1),t+1);
end



end



