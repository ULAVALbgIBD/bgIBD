

%%

function [output_dec_state error] = ibd_decode_path(paths, ob, viterbi_option)

output_dec_state = [];
error = 0;

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

oblist = ob;    %no need to map, ready states
I_prob = paths.pr;
T_prob = paths.tpl_log;
O_prob = paths.epl_log;
[d1, d2, d3] = size(T_prob);
if( d1 <= 0 || d2 <= 0 || d2 ~= d3 )
    error = 1;
    disp('error in transition states');
    return;
end
nloc = d1;
nstates = d2;
[d1, d2, d3] = size(O_prob);
if( d1 ~= nloc || d2 ~= nstates || d3 <= 0 )
    error = 1;
    disp('error in emission states');
    return;
end
len = length(I_prob);
if( len <= 0 || len ~= nstates )
    error = 1;
    disp('error in prior probability');
    return;
end
len = length(oblist);
if( len <= 0 || len ~= nloc )
    error = 1;
    disp('error in observation list');
    return;
end

hidden_states = [1:nstates];


[output_all_dec_state, output_del error] = get_viterbi_log(I_prob, T_prob, O_prob, oblist, hidden_states);
if( error ~= 0 )
    disp('error in generating viterbi decoding');
    return;
end


% 2 states
output_dec_state = separate_bg(output_all_dec_state, viterbi_option);

end





















