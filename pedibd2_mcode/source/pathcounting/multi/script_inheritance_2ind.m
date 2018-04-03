



function ih = script_inheritance_2ind(input_pedigree, family_id, ind1, ind2)

    ih.descent_graph = descent_graph_2ind(input_pedigree, family_id, ind1, ind2);
    ih.num = ih.descent_graph(1,1);
    [ih.is ih.transition] = hasse_matrix(ih.num);
    [ih.prob_f, ih.tran_list] = condense_inheritance(ih.descent_graph);
        
end









