

% get transition probability

%%

bg_ibd(s_id) = 0.0113; %estimated from the whole population
% bg_tran(s_id) = 0.9886;
bg_tran(s_id) = 400;    %generations apart

% bg_ibd = 0.2; %estimated from the whole population
% bg_tran = 0.3;


%%

for i = 1:1

    script_transition;

    %%%%%%%%%estimate background ibd%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    % oblist = input_oblist_genotype;
    oblist = input_oblist_ibs{s_id} + 1;
    I_prob = input_pr';
    T_prob = input_tpl;
    O_prob = [];
    for i = 1:length(input_epl)
        O_prob{i} = input_epl{i}';
    end
    hidden_states = [1:length(input_epl{1}(:,1))];

    est_alpha = get_forward(I_prob, T_prob, O_prob, oblist, hidden_states);
    est_beta = get_backward(I_prob, T_prob, O_prob, oblist, hidden_states);
    [est_bg, est_tran, est_alphabeta] = get_estimate2(I_prob, T_prob, O_prob, oblist, hidden_states, est_alpha, est_beta, sampled_markerlist);
    
    bg_ibd(s_id) = est_bg;
    bg_tran(s_id) = est_tran;
    bg_tran(s_id) = 400;    %generations apart
    
    clear i;
    clear I_prob T_prob O_prob states oblist hidden_states;

end

%%


%%





















