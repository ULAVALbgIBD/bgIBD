function [dec_state, del] = get_viterbi_locus(I_prob, T_prob, O_prob, oblist, states)

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

del=zeros(Nst,lob);% collection of the maximum probability values at each stage
maxlist=del;
mx=zeros(1,lob);
for t=1:lob
    if t==1
        del(:,t)=I_prob.*O_prob{t}(oblist(t),:)'; %Initialization
       [p mx(t)]=max(del(:,t));
        continue;
    end
    %Recursive Phase
    
    for j=1:Nst
        [del(j,t) maxlist(j,t)]=max(del(:,t-1).*T_prob{t}(:,j));
    end
    del(:,t)=del(:,t).*O_prob{t}(oblist(t),:)';
    [p mx(t)]=max(del(:,t)); % Read 4) in readme
    if( sum(del(:,t)) < 0.01 ) %&& (t > 11060 || t < 10960) )
        del(:,t) = del(:,t) * 100;
    end
end

        % Termination and Backtrack stage
dec_state=zeros(1,lob);
% decode_state=cell(1,lob);
[pstar dec_state(lob)]=max(del(:,lob));

for t=lob-1:-1:1
    dec_state(t)=maxlist(dec_state(t+1),t+1);
end



