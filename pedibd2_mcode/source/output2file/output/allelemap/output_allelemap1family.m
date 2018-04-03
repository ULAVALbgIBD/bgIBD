function error = output_allelemap1family(filehead, family, assignment, parameters)

suc = 0;
error = 0;


[nind ~] = size(family.structure);

alleles_all = assignment.alleles_all;
if( ndims(alleles_all) ~= 3 )
    error = 1;
    disp('error in global IBD');
    return;
else
    [nseg, d2, d3] = size(alleles_all);
    if( nseg <= 0 || d2 ~= nind || d3 ~= 2 )
        error = 1;
        disp('error in global IBD');
        return;
    end
end


len1 = length(assignment.score.p_value);
if( nseg ~= len1 )
    error = 1;
    disp('error in NPL p value');
    return;
end

intervals = assignment.intervals;
[nrows, ncols] = size( intervals );
if( nrows ~= nseg )
    error = 1;
    disp('error in chromosome segmentation');
    return;
end

if( ncols < 2 )
    error = 1;
    disp('error in chromosome segmentation');
    return;
end

markerloc = parameters.sampled_markerlist(:,2);

len2 = length(markerloc);

if( max(max(intervals(:,1:2))) > len2 || min(min(intervals(:,1:2))) < 1 )
    error = 1;
    disp('error in marker location');
    return;
end

chr = parameters.chr;

valid = assignment.recombination.valid;
len3 = length(valid);
if( len3 <= 0 || len3 ~= nseg )
    error = 1;
    disp('error in recombination positions');
    return;
end


seg_count = 0;
for i = 1:nseg

    if( valid(i) ~= 1 )
        continue;
    end
    seg_count = seg_count + 1;
    alleles_all = reshape(assignment.alleles_all(i,1:nind,1:2), [nind,2]);
    p_value = assignment.score.p_value(i);
    basepair = markerloc(intervals(i,1:2));
    
    filename = [filehead, '_seg', num2str(seg_count), '.png']; 
    
    error = output_allelemap1family1seg(filename, alleles_all, family, chr, basepair, p_value);
    if( error ~= 0 )
        disp('error in generating alleleMap');
        return;
    end

end

end














