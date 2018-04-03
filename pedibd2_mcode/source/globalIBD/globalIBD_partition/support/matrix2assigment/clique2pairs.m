
function output_pairs = clique2pairs(clique_assignment)
    for i = 1:length(clique_assignment(1,:))
        for j = 1:length(clique_assignment(1,:))
            output_pairs(i,j) = sum((clique_assignment(:,i) & clique_assignment(:,j)));
        end
    end
end