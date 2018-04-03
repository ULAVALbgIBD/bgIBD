function [dec_state, del] = get_viterbi_locus_log(identity, oblist)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



prob_f = identity.prob_f;
prob_b = identity.prob_b;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nst = nnz(prob_f)*nnz(prob_b);
lob=length(oblist);
del=zeros(Nst,lob);% collection of the maximum probability values at each stage
maxlist=del;
mx=zeros(1,lob);
dec_state_temp=zeros(1,lob);
dec_state=zeros(lob,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ct_combine = 0;
ct_decode = 0;


for t=1:lob  
       
    
    tran_matrix_f = identity.tran_f{t};
    tran_matrix_b = identity.tran_b{t};
       
    time = cputime;
        
    [margin_prob_log, tran_matrix_log, ix] = combine_fb(prob_f, tran_matrix_f, prob_b, tran_matrix_b);
    
    ct_combine = ct_combine + cputime - time;
    
    time = cputime;
        
    O_prob_l = log(identity.emission{t});
    O_prob_l = O_prob_l(ix(:,3),:);
    I_prob_l = margin_prob_log';
    T_prob_l = tran_matrix_log; % log transformation done inside
    
    if t==1
        del(:,t)=I_prob_l + O_prob_l(:,oblist(t)); %Initialization
        [p mx(t)]=max(del(:,t));
        continue;
    end
    %Recursive Phase
    
    
    for j=1:length(T_prob_l(1,:))
        [del(j,t) maxlist(j,t)]=max(del(:,t-1) + T_prob_l(:,j));
    end
    del(:,t)=del(:,t) + O_prob_l(:,oblist(t));
    [p mx(t)]=max(del(:,t)); % Read 4) in readme
%     if( sum(del(:,t)) < -1000 )
%         del(:,t) = del(:,t) + 1000;
%     end

    ct_decode = ct_decode + cputime - time;
    
end

ct_decode
ct_combine
% Termination and Backtrack stage

[pstar dec_state_temp(lob)]=max(del(:,lob));

for t=lob-1:-1:1
    dec_state_temp(t)=maxlist(dec_state_temp(t+1),t+1);
end

for t = 1:lob
    dec_state(t,1:2) = ix(dec_state_temp(t),1:2);
end









