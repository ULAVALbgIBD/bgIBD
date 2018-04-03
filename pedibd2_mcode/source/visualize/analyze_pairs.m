function analyze_pairs(paths, family_range, family_id)


if( isempty(paths) )
    return;
end
if( isempty(family_range) || isempty(family_range{family_id}) )
    return;
end

posterior = paths.posterior{family_id};
viterbi = paths.viterbi{family_id};

pairs = family_range{family_id}.pairs;
family = family_range{family_id}.structure;

[npairs, fields] = size(pairs);

if( fields ~= 2 || npairs <= 0 )
    return;
end

len1 = length(posterior);
if( len1 ~= npairs )
    return;
end

len2 = length(viterbi);
if( len2 ~= npairs )
    return;
end

list = [];
count = 0;
len_nonibd = 0;
len_bgibd = 0;
for i = 1:npairs
    temp_vit = [];
    vitf = (viterbi{i}(:,1)==2) | (viterbi{i}(:,1)==3);
    vitb = (viterbi{i}(:,2)==2) | (viterbi{i}(:,2)==3);
    len_nonibd = len_nonibd + nnz(vitf == 0);
    len_bgibd = len_bgibd + nnz(vitb > 0);
    vit = (viterbi{i}(:,1)==2) | (viterbi{i}(:,1)==3);
    temp_vit = vit;
    temp_vit(:,2) = 0;
    pos = posterior{i}(:,2)+posterior{i}(:,3);
    markers =  generate_consistent_interval_viterbi( temp_vit );
    [r, c] = size(markers);
    if( r <= 0 || c ~= 3 )
        return;
    end
    for j = 1:r
        if( vit(markers(j,3)) == 1 )
            count = count + 1;
            list(count,1) = markers(j,2) - markers(j,1) + 1;
            list(count,2) = mean(pos(markers(j,1):markers(j,2)));
            list(count,3) = family(pairs(i,1),8);
            list(count,4) = family(pairs(i,2),8);
        end
    end
end

scatter(list(:,1), list(:,2));
title(['average IBD interval length ', num2str(mean(list(:,1))), ' markers']);
disp(['background ibd level ', num2str(len_bgibd/len_nonibd)]);

end
















