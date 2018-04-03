function smoothed = smooth(data, interval)

smoothed = [];
if( isempty(data) )
    return;
end

[nrow, ncol] = size(data);
if( nrow ~= 1 && ncol ~= 1 )
    return;
    % data must be one-demensional
end

if( length(data) < 2*interval + 1 )
    interval = floor((length(data)-1)/2);
end

smoothed(1:length(data)) = data + 0;
s = sum(data(1:2*interval+1));

smoothed(1+interval) = s/(2*interval+1);
for i = 1+interval+1:length(data)-interval
    s = s - data(i-interval-1);
    s = s + data(i+interval);
    smoothed(i) = s/(2*interval+1);
end

smoothed(1:interval) = smoothed(1+interval);
smoothed(length(data)-interval+1:end) = smoothed(length(data)-interval);



end