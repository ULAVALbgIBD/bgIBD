function clique_grow = aggressive_grow(pair, clique)
    
clique_grow = clique;
for i = 1:length(pair)
    if( sum(pair(i,clique_grow)>=1) == length(clique_grow) )
        clique_grow = union(clique_grow, i);
    end
end
    
end