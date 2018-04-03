function [ interval error ] = quick_interval(viterbi)

interval = [];
error = 0;

[nloc c] = size(viterbi);

if( nloc <= 0 || c ~= 1 )
    error = 1;
    return;
end

value = unique(viterbi(1:nloc));

cur(1:nloc+1) = 0;
cur(1:nloc) = viterbi(1:nloc);
cur(nloc+1) = -1;

pre(1:nloc+1) = 0;
pre(1) = -1;
pre(2:nloc+1) = viterbi(1:nloc);

change(1:nloc+1) = ( pre(1:nloc+1) ~= cur(1:nloc+1) );

points = find(change(1:nloc+1));

len = length(points);
interval(1:len-1, 1:3) = 0;

for i = 1:len
    if( i > 1 )
        interval(i-1, 2) = points(i) - 1;
    end
    if( i < len )
        interval(i, 1) = points(i);
    end
end

interval(1:len-1,3) = round(mean(interval(1:len-1,1:2), 2));

end


