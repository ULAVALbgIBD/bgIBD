


function [error] = output_linkage_score(fid, family_range, combined_score, parameters)

% output linkage score, all segments combined
error = 0;

if( isempty(combined_score) )
    error = 1;
    disp('linkage score not generated');
    return;
end

if( isempty(parameters) )
    error = 1;
    return;
end

sampled_markerlist = parameters.sampled_markerlist;
chr = parameters.chr;

[nmarkers, fields] = size( sampled_markerlist );
if( nmarkers < 1 || fields ~= 2 )
    error = 1;
    return;
end

if( chr < 0 )
    error = 1;
    return;
end

all_p_value = combined_score.all_p_value;
all_intervals = combined_score.all_intervals;
combined_p_value = combined_score.p_combined;

[nint, col] = size(all_intervals);
if( col ~= 2 )
    error = 1;
    return;
end

if( nint < 1 )
    error = 1;
    return;
end

if( all_intervals(1,1) < 1 || all_intervals(nint,2) > nmarkers )
    error = 1;
    return;
end

nfam = length(family_range);
if( nfam <= 0 )
    error = 1;
    disp('no families');
    return;
end

[r, c] = size(all_p_value);
if( r ~= nint || c ~= nfam )
    error = 1;
    return;
end

len = length(combined_p_value);
if( len ~= nint )
    error = 1;
    return;
end


fprintf(fid, '*****       non-parametric linkage       *****\n\n\n\n');

fprintf(fid, '                                     p_value  \t\t  family\n');
fprintf(fid, 'chr                    region        combined \t\t');
for i = 1:nfam
    fprintf(fid, '%8d\t', family_range{i}.family_id);
end
fprintf(fid, '\n');

map = sampled_markerlist(1:nmarkers, 2);
for i = 1:nint
    s = all_intervals(i,1);
    t = all_intervals(i,2);
    fprintf(fid, '%3d   ', chr);
    fprintf(fid, '%10d ~ %10d        ', map(s), map(t)); 
    fprintf(fid, '%.2e\t\t', combined_p_value(i));
    for j = 1:nfam
        fprintf(fid, '%.6f\t', all_p_value(i,j));
    end
    fprintf(fid, '\n');
end





end













