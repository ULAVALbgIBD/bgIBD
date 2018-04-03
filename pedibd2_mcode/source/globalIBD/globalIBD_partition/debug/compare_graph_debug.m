% find the maximum clique of a graph

function compare_graph_debug(assignment1, assignment2)




for seg = 1:length(assignment1.alleles)
    
    difference(seg) = nnz(allele2pair(assignment1.alleles{seg}) - allele2pair(assignment2.alleles{seg}));  
    line(assignment1.intervals(seg,1:2), [difference(seg), difference(seg)], 'color', 'b', 'linewidth', 3);
    text(assignment1.intervals(seg,3), difference(seg), num2str(seg));
    difference1vit(seg) = assignment1.diff(seg);
    difference2vit(seg) = assignment2.diff(seg);
    line(assignment1.intervals(seg,1:2), [difference1vit(seg), difference1vit(seg)], 'color', 'r', 'linewidth', 2);
end


end

