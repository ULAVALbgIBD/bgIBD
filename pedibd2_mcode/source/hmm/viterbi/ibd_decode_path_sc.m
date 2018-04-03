

%%


function output_dec_state_sc = ibd_decode_path_sc(paths, ob, pm)


%%%%%%%%%locus specific decoding%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% single chain decoding using genotype as emission states, 4 states

oblist = ob;
I_prob = paths.pr_sc';
T_prob = paths.tpl_sc;
O_prob = [];
for i = 1:length(paths.epl_sc)
    O_prob{i} = paths.epl_sc{i}';
    
end

hidden_states = [1:length(paths.epl_sc{1}(:,1))];

[temp, output_del] = get_viterbi_log(I_prob, T_prob, O_prob, oblist, hidden_states);



output_all_dec_state_sc = temp;
output_dec_state_sc = separate_bg_sc(temp, pm.sampled_markerlist);



end





















