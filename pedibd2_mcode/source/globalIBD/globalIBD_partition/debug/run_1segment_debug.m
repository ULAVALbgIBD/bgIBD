
function run_1segment_debug(input_pairs, input_triples, seg, assignment, range, kinship2ex, singleLINEAL)


input_pairs.intervals(seg,[1,2])
if( ~isempty(assignment) )
    verify = squeeze(assignment.alleles(seg,:,:));
else
    verify = [];
end
debug = 1;
seg_triples.num_triples = input_triples.num_triples;
seg_triples.ambig_triples = input_triples.ambig_triples;
seg_triples.triple_viewpair = squeeze(input_triples.triple_viewpair(seg,:,:,:));

[alleles_all error] = ...
    partition_1segment( ...
    squeeze(input_pairs.pair_vit(seg,:,:)), ...
    squeeze(input_pairs.pair_max(seg,:,:)), ...
    kinship2ex, ...
    singleLINEAL, ...
    squeeze(input_pairs.pair_posIBD1(seg,:,:,:,:)), ...
    squeeze(input_pairs.pair_pos(seg,:,:,:)), ...
    seg_triples, ...
    range, verify, debug);

if( error ~= 0 )
    disp('error in DNAphoresis');
    return;
end

[alleles_all error] = imputation(alleles_all, range, kinship2ex);
alleles = alleles_all(range.family_range,1:2);
inferred_pairs = allele2pair(alleles);

intervals = input_pairs.intervals;

title(['segment', num2str(seg), num2cell(intervals(seg,[1,2]))]);

set(gcf, 'Position',[300 300 300 + 1200 300 + 400], 'Visible', 'on');

if( ~isempty(verify) )
    pairs = allele2pair(verify);
    disp(['likelihood to viterbi difference: ', num2str(nnz(squeeze(input_pairs.pair_vit(seg,:,:)) - pairs))]);
    disp(['graph to viterbi difference: ', num2str(nnz(squeeze(input_pairs.pair_vit(seg,:,:)) - inferred_pairs))]);
    disp(['graph to likelihood difference: ', num2str(nnz(pairs - inferred_pairs))]);
end

end