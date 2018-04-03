function [pairs corrected_viterbi error] ...
    = generate_consistent_pairs(oblist, viterbi, posterior, posteriorIBD1, ...
    kinship2ex, range, parameters)

error = 0;
pairs = [];
corrected_viterbi = int8([]);

[nmarkers, cols] = size(parameters.sampled_markerlist);

if( nmarkers < 1 || cols ~= 2 )
    error = 1;
    disp('error in marker information');
    return;
end

if( isempty(range) )
    error = 1;
    disp('empty family');
    return;
end

ngeno = length(range.family_range);

if( ngeno < 2 )
    if( ~isempty(viterbi) || ~isempty(posterior) || ~isempty(posteriorIBD1) )
        error = 1;
        disp('error in family structures');
        return;
    end
    if( ~isempty(range.reverse_pairs) )
        error = 1;
        disp('error in family structures');
        return;
    end
    
    disp(['family ', num2str(range.family_id), ': < 2 genotyped individuals: segmentation skipped']);
    pairs = [];
    if( ngeno == 1 )
        pairs.pair_vit(1,1,1) = 2;       
        pairs.pair_pos(1,1,1,1:3) = [0,0,1];
        pairs.pair_posIBD1(1,1,1,1:2,1:2) = [1,0;0,1];
        pairs.pair_max(1,1,1) = 2;
    end
    pairs.intervals(1,1:3) = [1, nmarkers, (1+nmarkers)/2];
    return;
end


display(' ');
disp('segmenting the chromosome ...');

time = cputime;

[npairs, c] = size(range.pairs);
if( npairs <= 0 || c ~= 2 )
    error = 1;
    disp('error in generating affinity pairs');
    return;
end

if( ndims(viterbi) ~= 3 )
    error = 1;
    disp('error in generating affinity pairs');
    return;
end
[d1 d2 d3] = size(viterbi);
if( d1 ~= npairs || d2 ~= nmarkers || d3 ~= 2 )
    error = 1;
    disp('error in generating affinity pairs');
    return;
end

if( ndims(posterior) ~= 3 )
    error = 1;
    disp('error in generating affinity pairs');
    return;
end
[d1 d2 d3] = size(posterior);
if( d1 ~= npairs || d2 ~= nmarkers || d3 ~= 3 )
    error = 1;
    disp('error in generating affinity pairs');
    return;
end

if( ndims(posteriorIBD1) ~= 4 )
    error = 1;
    disp('error in generating affinity pairs');
    return;
end
[d1 d2 d3 d4] = size(posteriorIBD1);
if( d1 ~= npairs || d2 ~= nmarkers || d3 ~= 2 || d4 ~= 2 )
    error = 1;
    disp('error in generating affinity pairs');
    return;
end

if( ndims(oblist) ~= 2 )
    error = 1;
    disp('error in generating affinity pairs');
    return;
end
[d1 d2] = size(oblist);
if( d1 ~= npairs || d2 ~= nmarkers )
    error = 1;
    disp('error in generating affinity pairs');
    return;
end

[smoothed_viterbi error] = smooth_viterbi(oblist, viterbi, posterior, posteriorIBD1, parameters.sampled_markerlist(1:nmarkers,2));
if( error ~= 0 )
    disp('error in segmentation');
    return;
end

% rounding applied here to save memory
[corrected_viterbi error] = correct_viterbi(smoothed_viterbi, kinship2ex, range, parameters.sampled_markerlist(1:nmarkers,2));
if( error ~= 0 )
    disp('error in segmentation');
    return;
end


% output original posterior in the output file
[posterior, posteriorIBD1 error] = smooth_posterior(corrected_viterbi, posterior, posteriorIBD1, nmarkers);

if( error == 0 )
    [ pairs error ] = generate_segments( corrected_viterbi, posterior, posteriorIBD1, range.reverse_pairs );
    if( isempty(pairs) || error ~= 0 )
        disp('error in segmentation');
        error = 1;
        return;
    end
else
    disp('error in segmentation');
    error = 1;
    return;
end

    
display(['segmentation costs time: ', num2str(cputime - time), ' seconds']);  
  

end






