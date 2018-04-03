

function all_data = transposeGENOTYPE(genotype)

[nLOC, nIND, c] = size(genotype);

all_data = zeros(nIND, nLOC * c);
for i = 1:c
    all_data(1:nIND, i:c:end) = genotype(1:nLOC, 1:nIND, i)';
end


end