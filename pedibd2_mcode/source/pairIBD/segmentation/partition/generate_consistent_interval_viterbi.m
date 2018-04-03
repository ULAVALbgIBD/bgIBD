function [ interval ] = generate_consistent_interval_viterbi( viterbi )

% only viterbi foreground probabilities are considered

all_interval = [];

if( ndims(viterbi) == 3 )
    [npairs nloc cols] = size(viterbi);
else
    npairs = 1;
    [nloc cols] = size(viterbi);
    temp = zeros(1,nloc,cols);
    temp(1,1:nloc,1:cols) = viterbi;
    viterbi = temp;
end

if( npairs <= 0 )
    disp('error in viterbi decoding');
    return;
end

if( nloc <= 0 || cols ~= 2 )
    disp('error in viterbi decoding');
    return;
end

cur(1:npairs, 1:nloc+1) = 0;
pre(1:npairs, 1:nloc+1) = 0;

cur(1:npairs, 1:nloc) = viterbi(1:npairs,1:nloc,1);
cur(1:npairs, nloc+1) = -1;
pre(1:npairs, 2:nloc+1) = cur(1:npairs, 1:nloc);
pre(1:npairs, 1) = -1;

change(1:npairs, 1:nloc+1) = (pre ~= cur);
anychange(1:nloc+1) = any(change, 1);

points = find(anychange);
len = length(points);
interval(1:len-1, 1:3) = 0;

interval(1:len-1, 2) = points(2:len) - 1;
interval(1:len-1, 1) = points(1:len-1);

interval(1:len-1,3) = round(mean(interval(1:len-1,1:2), 2));

end




















