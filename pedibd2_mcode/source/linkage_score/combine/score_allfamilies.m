function [combined_score, error] = score_allfamilies(assignment, families)

% correctness check will fail in debug mode
% if only the first family is processed

% consider each p value as a uniform distribution, 
% build a final statistic

error = 0;
combined_score = [];

[combined_score, error] = combine_pvalue(assignment);
if( error ~= 0 )
    disp('error in calculating p value');
    return;
end

z_score = combined_score.all_z_score;
p_value = combined_score.all_p_value;

nfam = length(families);
[nrows, ncols] = size(z_score);
if( ncols ~= nfam || ncols < 1 )
    error = 1;
    disp('not all families processed');
    return;
end
nint = nrows;
[r, c] = size(p_value);
if( c ~= nfam || r < 1 || r ~= nint )
    error = 1;
    disp('not all families processed');
    return;
end


nGENOaffected(1:nfam) = 0;
for i = 1:nfam
    if( isempty(families) )
        error = 1;
        disp('error in family structures');
        return;
    end
    family = families{i}.structure;
    [r, c] = size(family);
    if( r <= 0 || c < 12 )
        error = 1;
        disp('error in family structures');
        return;
    end
    count = 0;
    for j = 1:r
        if( family(j,7) == 1 && family(j,6) == 2 )
            count = count + 1;
        end
    end
    nGENOaffected(i) = count;
end

for j = 1:nfam
    if( sum(p_value(1:nint, j)) == nint )
        % do not count nonvariance families
        nGENOaffected(j) = 0;
    end
end

weight(1:nfam) = 0;
total = sum(nGENOaffected(1:nfam));
if( total > 0 )
    for i = 1:nfam
        weight(i) = (nGENOaffected(i)/total)^0.5;
    end
else
    for i = 1:nfam
        weight(i) = 0;
    end
end
% squared summation of weight == 1

p_combined = zeros(nrows, 1);
z_combined = zeros(nrows, 1);
u_combined = zeros(nrows, 1);
% combined value of p, uniform distribution

for i = 1:nrows
    temp_z = 0;
    temp_u = 0;
    for j = 1:nfam
        temp_z = temp_z + z_score(i,j) * weight(j);
        temp_u = temp_u + weight(j) * ...
            (p_value(i,j) - 0.5) * 2 * sqrt(3);
    end
    z_combined(i) = temp_z;
    u_combined(i) = temp_u;
end

% final statistic is standard normal
% from z value
if( total > 0 )
    for i = 1:nrows
        p_combined(i) = 1 - normcdf(z_combined(i),0,1);
    end
else
    for i = 1:nrows
        p_combined(i) = 1;
    end
end

% final statistic is standard normal
% from p value
% if( total > 0 )
%     for i = 1:nrows
%         p_combined(i) = normcdf(u_combined(i),0,1);
%     end
% else
%     for i = 1:nrows
%         p_combined(i) = 1;
%     end
% end

if( nnz(weight) == 1 )
    p_combined = prod(p_value,2);
end

combined_score.p_combined = p_combined;

end




















