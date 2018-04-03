

function [output_dec_state] = separate_bg(temp, ibd_option)


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

output_dec_state(1:length(temp),1:2) = 1;

output_dec_state(:,1) = ibd_option.foreground_map(temp);
output_dec_state(:,2) = ibd_option.background_map(temp);



end






