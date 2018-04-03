function output_segments(fid2, segments, markerlist)


temp6 = markerlist;

[nseg c] = size(segments);
if( c ~= 2 )
    return;
end

for j = 1:nseg
    a = segments(j,1);
    b = segments(j,2);
    fprintf(fid2, '%d\t%d\t%d\t%d\n', a, b, temp6(a,2), temp6(b,2));
end

fclose('all');

end
