function test_visualization(assignment, family, parameters, seg)


suc = 0;
error = 0;

if( exist([pwd, '\alleleMap'], 'dir') == 0 )
    suc = mkdir('alleleMap');
else
    suc = 1;
end

if( suc ~= 1 )
    error = 1;
    disp('cannot create output folder');
    return;
end

filename = ['alleleMap/', 'myplot.png'];

alleles_all = squeeze(assignment.alleles_all(seg,:,:));
p_value = assignment.score.p_value(seg);
chr = parameters.chr;
basepair = parameters.sampled_markerlist(assignment.intervals(seg,1:2),2);


error = output_allelemap1family1seg(filename, alleles_all, family, chr, basepair, p_value);

end




