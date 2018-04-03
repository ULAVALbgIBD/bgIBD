function [posterior, posteriorIBD1 error] = smooth_posterior(viterbi, posterior, posteriorIBD1, nmarkers)

global debug_mode;


% due to tremendous size of this two variables, do not duplicate for output

error = 0;

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

time = cputime;


for i = 1:npairs
    markers =  generate_consistent_interval_viterbi( viterbi(i, 1:nmarkers, 1:2) );
    [d1, d2] = size(markers);
    if( d1 <= 0 || d2 ~= 3 || markers(d1,2) ~= nmarkers )
        disp('error in segmentation');
        error = 1;
        return;
    end

    for j = 1:d1
        for k = 1:3
            posterior(i, markers(j,1):markers(j,2),k) = mean(posterior(i, markers(j,1):markers(j,2),k));
        end
        for p = 1:2
            for q = 1:2
                posteriorIBD1(i, markers(j,1):markers(j,2),p,q) = mean(posteriorIBD1(i, markers(j,1):markers(j,2),p,q));
            end
        end
    end

end

if( debug_mode == 1 )
    display(['          posterior smoothing costs time: ', num2str(cputime - time), ' seconds']);  
end

end







