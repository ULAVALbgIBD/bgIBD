function [smoothed_viterbi error] = smooth_viterbi(oblist, viterbi, posterior, posteriorIBD1, map)

smoothed_viterbi = [];
error = 0;
global debug_mode;



time = cputime;

[nmarkers, c] = size(map);
if( nmarkers <= 0 || c ~= 1 )
    error = 1;
    disp('error in marker list');
    return;
end

if( ndims(viterbi) ~= 3 )
    error = 1;
    disp('error in generating affinity pairs');
    return;
end
[d1 d2 d3] = size(viterbi);
npairs = d1;
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

smoothed_viterbi = zeros(npairs, nmarkers, 2); 
% pre-assign to speed up
% flatten the representation to speed up

for i = 1:npairs
    ob = reshape(oblist(i,1:nmarkers), [nmarkers,1]);
    vit = reshape(viterbi(i,1:nmarkers,1:2), [nmarkers,2]);
    pos = reshape(posterior(i,1:nmarkers,1:3), [nmarkers,3]);
    pos1 = reshape(posteriorIBD1(i,1:nmarkers,1:2,1:2), [nmarkers,2,2]);
    if( i == 48 )
    end
    [smoothed_viterbi(i, 1:nmarkers, 1:2) error] = smooth_viterbi1pair(ob, vit, pos, pos1, map);
    if( error ~= 0 )
        disp('error in generating affinity pairs');
        return;
    end
end

if( debug_mode == 1 )
    display(['          ibd smoothing costs time: ', num2str(cputime - time), ' seconds']);  
end

end


































