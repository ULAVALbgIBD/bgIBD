
function [triples, triples_view, error] ...
    = valid_triples(family_range, kinship, kinship2ex, ...
    family_genotype, viterbi, intervals)

error = 0;
triples_view = [];
triples = [];

if( isempty(intervals) )
    error = 1;
    return;
end
[nSEG, c] = size(intervals);
if( nSEG <= 0 || c ~= 3 )
    error = 1;
    disp('error in segmentation');
    return;
end
nIND = length(family_range.pedigree_range_full);
if( nIND <= 0 )
    error = 1;
    disp('error in family structures');
    return;
end
genotyped = family_range.family_range;
if( any(genotyped > nIND) )
    error = 1;
    disp('error in family structures');
    return;
end
nGENO = length(genotyped);

if( nGENO < 2 )
    return;
end

display(' ');
display(['family ',num2str(family_range.family_id), ': removing third chromosome sharing ...']);
time = cputime;

kinship2 = kinship2ex(genotyped,genotyped,1:2,1:2);
[triples, error] = generate_triples(family_range, kinship, kinship2);
if( error ~= 0 )
    disp('error in separating');
    return;
end
[triples_view, error] ...
    = smooth_triples(family_range, triples, family_genotype, viterbi, intervals);
if( error ~= 0 )
    disp('error in separating');
    return;
end

if( nGENO > 0 && any(size(triples_view.triple_viewpair) ~= [nSEG,nGENO,nGENO,nGENO]) )
    error = 1;
    disp('error in triplets');
    return;
end


[error] = verify_triples(triples_view, nSEG, nGENO);
if( error ~= 0 )
    disp('error in separating');
    return;
end



display(['removing costs: ', num2str(cputime - time), ' seconds']);

end


