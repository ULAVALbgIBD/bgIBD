function ibd_option = generate_ibd_option()

% currently static, modify it to be more versatile
% use identity states for more universal usage

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


num_states = length(states(:,1));
% number of output states
foreground_map(1:num_states) = 0;
background_map(1:num_states) = 0;



for i = 1:num_states
    foreground_map(i) = nnz(states(i,:)==2)+1;
    background_map(i) = nnz(states(i,:)==3)+1;
end

ibd_option.viterbi.foreground_map = foreground_map;
ibd_option.viterbi.background_map = background_map;





end