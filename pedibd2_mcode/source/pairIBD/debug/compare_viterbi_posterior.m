function compare_viterbi_posterior(pairs, intervals)

ml = pairs.pair_max;
vi = pairs.pair_vit;

[rows, cols] = size(intervals);

if( cols < 2 )
    disp('error in segmentation');
    return;
end
if( rows ~= length(ml) )
    disp('error in output pairs');
    return;
end
if( rows ~= length(vi) )
    disp('error in output pairs');
    return;
end

diff(1:rows) = 0;

for i = 1:rows
    if( length(ml{i}) ~= length(vi{i}) )
        disp('error in output pairs');
        return;
    end
    diff(i) = nnz(ml{i}-vi{i});
end

for i = 1:rows

    line(intervals(i,1:2), [diff(i), diff(i)], 'color', 'r', 'linewidth', 2);
    text(intervals(i,3), diff(i), num2str(i));

end

end