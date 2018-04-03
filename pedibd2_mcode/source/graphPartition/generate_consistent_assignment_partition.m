function [assignment error] ...
    = generate_consistent_assignment_partition ...
    (input_pairs, ...
    input_triples, ...
    kinship2, ...
    singleLINEAL, ...
    verify, debug, range)

%%

assignment = [];
error = 0;

if( isempty(input_pairs) || isempty(kinship2) || isempty(input_triples) )
    error = 1;
    disp('DNA-phoresis for 3 and more genotyped individuals');
    return;
end

intervals = input_pairs.intervals;

[nSEG, ncols] = size(intervals);
if( nSEG <= 0 || ncols ~= 3 )
    error = 1;
    return;
end

[nIND, ~] = size(range.structure);
nGENO = length(range.family_range);


pair_pos = input_pairs.pair_pos;
pair_max = input_pairs.pair_max;
pair_posIBD1 = input_pairs.pair_posIBD1;
pair_vit = input_pairs.pair_vit;

if( ndims(pair_pos) ~= 4 || any(size(pair_pos) ~= [nSEG, nGENO, nGENO, 3]) )
    error = 1;
    return;
end

if( ndims(pair_max) ~= 3 || any(size(pair_max) ~= [nSEG, nGENO, nGENO]) )
    error = 1;
    return;
end

if( ndims(pair_vit) ~= 3 || any(size(pair_vit) ~= [nSEG, nGENO, nGENO]) )
    error = 1;
    return;
end

if( ndims(pair_posIBD1) ~= 5 || any(size(pair_posIBD1) ~= [nSEG, nGENO, nGENO, 2, 2]) )
    error = 1;
    return;
end

if( nGENO > 0 && any(size(input_triples.triple_viewpair) ~= [nSEG,nGENO,nGENO,nGENO]) )
    error = 1;
    return;
end

%%

alleles_all = zeros(nSEG, nIND, 2);

seg_triples.num_triples = input_triples.num_triples;
seg_triples.ambig_triples = input_triples.ambig_triples;
for seg = 1:nSEG
    disp(['          processing segment: ', num2str(seg), ' ...']);
    % false positive, cliques not included in the correct ones
    
    seg_triples.triple_viewpair ...
        = reshape(input_triples.triple_viewpair(seg,1:nGENO,1:nGENO,1:nGENO), ...
        [nGENO, nGENO, nGENO]);
    
    [alleles error] ...
        = partition_1segment(reshape(pair_vit(seg, 1:nGENO, 1:nGENO), [nGENO, nGENO]), ...
        reshape(pair_max(seg, 1:nGENO, 1:nGENO), [nGENO,nGENO]), ...
        kinship2, ...
        singleLINEAL, ...
        reshape(pair_posIBD1(seg, 1:nGENO, 1:nGENO, 1:2, 1:2), [nGENO,nGENO,2,2]), ...
        reshape(pair_pos(seg, 1:nGENO, 1:nGENO, 1:3), [nGENO,nGENO,3]), ...
        seg_triples, ...
        range, verify, debug);
    
    if( error ~= 0 )
        disp('error in global IBD');
        return;
    end
    
    if( isempty(alleles) )
        disp(['segment: ', num2str(seg), ' assignment error']);
    else
        alleles_all(seg, 1:nIND, 1:2) = alleles;
    end
end

assignment.alleles_all = alleles_all;
assignment.intervals = intervals;



end











