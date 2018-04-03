function [likelihood] = get_likelihood_log(I_prob_o, T_prob_o, O_prob_o, oblist, states, dec)

% get cumulated likelihood of hidden states
% for debug usage

Nst=length(states);


lob=length(oblist);

del(1:lob)=0;% collection of the maximum probability values at each stage


I_prob_l = log(I_prob_o);
for i = 1:length(T_prob_o)
    T_prob_l{i} = log(T_prob_o{i});
end
for i = 1:length(O_prob_o)
    O_prob_l{i} = log(O_prob_o{i});
end


for t=1:lob
    if t==1
        del(t) = I_prob_l(dec(t)) + O_prob_l{t}(oblist(t),dec(t)); %Initialization
        continue;
    end
    %Recursive Phase
    
    [del(t)] = del(t-1) + T_prob_l{t}(dec(t-1),dec(t));
    del(t)=del(t) + O_prob_l{t}(oblist(t),dec(t));

end



likelihood = del;

