function output_genF = descent_graph_combine(genFp, genFm)



% genF format
% content format
% first column is coefficient
% 2~slots+1, is the allele involved
% 2+slots~, is the power of each term
    
% consider using unique to simply

genF = [genFp;genFm];
[alleleANDpower, ~, ix] = unique(genF(:,2:end), 'rows');
% alleleANDpower(ix) = genF(:,2:end);
coef = accumarray(ix, genF(:,1), ...
    [size(alleleANDpower,1), 1], @sum, 0);
output_genF = [coef,alleleANDpower];















