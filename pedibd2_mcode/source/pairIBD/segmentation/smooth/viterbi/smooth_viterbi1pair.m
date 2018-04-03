function [smoothed_vit error] = smooth_viterbi1pair(oblist, viterbi, posterior, posteriorIBD1, map)

smoothed_vit = int8([]);
error = 0;

[nmarkers, c] = size(map);
if( nmarkers <= 0 || c ~= 1 )
    error = 1;
    disp('error in marker list');
    return;
end

[r, c] = size(viterbi);
if( r <= 0 || r ~= nmarkers || c ~= 2 )
    error = 1;
    disp('error in viterbi decoding');
    return;
end

segments =  generate_consistent_interval_viterbi( viterbi );
[nseg, c] = size(segments);
if( nseg <= 0 || nseg > nmarkers || c ~= 3 || segments(nseg,2) ~= nmarkers )
    error = 1;
    disp('error in segmenting chromosomes');
    return;
end

[r, c] = size(posterior);
if( r <= 0 || r ~= nmarkers || c ~= 3 )
    error = 1;
    disp('error in posterior decoding');
    return;
end

[d1, d2, d3] = size(posteriorIBD1);
if( d1 <= 0 || d1 ~= nmarkers || d2 ~= 2 || d3 ~= 2)
    error = 1;
    disp('error in posterior decoding');
    return;
end

len = length(oblist);
if( len <= 0 || len ~= nmarkers )
    error = 1;
    disp('error in observation list');
    return;
end

smoothed_vit = viterbi;

[smoothed_vit error] = resist_fluctuation(oblist, smoothed_vit, posteriorIBD1, posterior, map);

[smoothed_vit error] = remove_pikes(oblist, smoothed_vit, posteriorIBD1, posterior, map);

[smoothed_vit error] = suppress_ld(oblist, smoothed_vit, posterior, map);


end


































