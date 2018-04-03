function [triples_seg, error] = smooth_triples(family_range, triples, family_genotype, vit, intervals)


triples_seg = [];
error = 0;

global debug_mode;

if( isempty(family_range) )
    error = 1;
    return;
end
nind = length(family_range.pedigree_range_full);
if( nind <= 0 )
    error = 1;
    disp('error in family structures');
    return;
end
family = family_range.structure;
[r, c] = size(family);
if( r <= 0 || r ~= nind || c < 12 )
    error = 1;
    disp('error in family structures');
    return;
end
genotyped = family_range.family_range;
ngeno = length(genotyped);
if( any(genotyped > nind) || any(genotyped < 0 ) )
    error = 1;
    disp('error in family structures');
    return;
end
[nmarkers, d2, d3] = size(family_genotype);
if( d2 <= 0 || d2 ~= nind )
    error = 1;
    disp('error in genotype data');
    return;
end
if( nmarkers <= 0 || d3 ~= 2 )
    error = 1;
    disp('error in genotype data');
    return
end

[nseg, c] = size(intervals);
if( nseg <= 0 || c < 2 )
    error = 1;
    disp('error in chromosome segmentation');
    return;
end
for i = 1:nseg
    a = intervals(i,1);
    b = intervals(i,2);
    if( a <= 0 || a > b || b > nmarkers )
        error = 1;
        disp('error in chromosome segmentation');
        return;
    end
end

if( ngeno < 1 )
    return;
end

ambig_triples = triples.ambig_triples;
num_triples = triples.total_triples;
deter_triples = triples.deter_triples;

[d1, d2, d3] = size(ambig_triples);
if( d1 ~= d2 || d2 ~= d3 || d1 ~= ngeno )
    error = 1;
    disp('error in triplets');
    return;
end
if( ~isempty(ambig_triples) && max(max(max(ambig_triples))) > num_triples )
    error = 1;
    disp('error in triplets');
    return;
end
[d1, d2, d3] = size(deter_triples);
if( d1 ~= d2 || d2 ~= d3 || d1 ~= ngeno )
    error = 1;
    disp('error in triplets');
    return;
end
for i = 1:ngeno
    for j = 1:ngeno
        for k = 1:ngeno
            if( deter_triples(i,j,k) ~= 0 && ambig_triples(i,j,k) > 0 )
                error = 1;
                disp('error in triplets');
                return;
            end
        end
    end
end

if( ngeno < 1 )
    all_triplesview = [];
else
    all_triplesview = zeros(nseg, ngeno, ngeno, ngeno);
end


triples_seg.triple_viewpair = all_triplesview;
triples_seg.ambig_triples = ambig_triples;
triples_seg.deter_triples = deter_triples;    
triples_seg.num_triples = num_triples;


if( ngeno < 3 )   
    return;
end

pairs = family_range.pairs;
[npairs, d2] = size(pairs);
if( npairs <= 0 || d2 ~= 2 )
    error = 1;
    return;
end
reverse_pairs = family_range.reverse_pairs;
[r, c] = size(reverse_pairs);
if( r ~= ngeno || c ~= ngeno )
    error = 1;
    return;
end
if( ndims(vit) ~= 3 )
    error = 1;
    disp('error in viterbi decoding');
    return;
end    
[d1 d2 d3] = size(vit);
if( d1 ~= npairs || d2 ~= nmarkers  || d3 ~= 2 )
    error = 1;
    disp('error in viterbi decoding');
    return;
end

if( num_triples <= 0 )
    return;
end

time = cputime;

triple_viewpair = zeros(num_triples,nmarkers);

for i = 1:ngeno
    for j = i+1:ngeno
        index_pair = reverse_pairs(i,j);
        if( index_pair <= 0 || index_pair > npairs )
            error = 1;
            disp('error in family structures');
            return;
        end
        affinity1_pairs = (vit(index_pair,1:nmarkers,1) == 2);
        if( ~any(affinity1_pairs) )
            continue;
        end
        pair_ibs(1:nmarkers) = 0;
        ind1 = pairs(index_pair,1);
        ind2 = pairs(index_pair,2);
        if( ind1 ~= genotyped(i) || ind2 ~= genotyped(j) )
            error = 1;
            disp('error in family structures');
            return;
        end
        if( ~any(ambig_triples(:,i,j)) )
            continue;
        end
        [pair_ibs error] = count_ibs(reshape(family_genotype(1:nmarkers,ind1,1:2),[nmarkers,2]), reshape(family_genotype(1:nmarkers,ind2,1:2), [nmarkers,2]));
        if( error ~= 0 )
            disp('error in genotype data');
            return;
        end
        for k = 1:ngeno 
            triple_index = ambig_triples(k,i,j);            
            if( triple_index <= 0 )
                continue;
            end
            pivot = genotyped(k);
            if( pivot <= 0 || pivot > nind )
                error = 1;
                disp('error in family structures');
                return;
            end
            index1 = reverse_pairs(k,i);
            index2 = reverse_pairs(k,j);
            if( index1 <= 0 || index1 > npairs || index2 <= 0 || index2 > npairs )
                error = 1;
                disp('error in family structures');
                return;
            end
            view1 = (vit(index1,1:nmarkers,1) == 2);
            view2 = (vit(index2,1:nmarkers,1) == 2);
            triple_ibd = view1 & view2 & affinity1_pairs;
            if( ~any( triple_ibd ) )
                continue;
            end
            
            [pivot_hetero error] = count_hetero(reshape(family_genotype(1:nmarkers, pivot, 1:2), [nmarkers,2]));
            if( error ~= 0 )
                disp('error in genotype data');
                return;
            end
            pivot_hetero(pair_ibs == 3) = false;    %exclude missing
            pair_identical = false(nmarkers, 1);
            pair_identical(pivot_hetero & pair_ibs == 2) = true;
            
  
            temp_vit = vit([index_pair, index1, index2], 1:nmarkers, 1:2);          
            regions = generate_consistent_interval_viterbi(temp_vit);
            [nr, c] = size(regions);
            if( nr <= 0 || c ~= 3 )
                error = 1;
                disp('error in segmentation');
                return;
            end
            if( regions(nr,2) ~= nmarkers )
                error = 1;
                disp('error in segmentation');
                return;
            end
            for t = 1:nr
                a = regions(t,1);
                b = regions(t,2);
                if( b < a || a <= 0 || b > nmarkers )
                    error = 1;
                    disp('error in segmentation');
                    return;
                end
                if( any(triple_ibd(a:b)) )
                    if( ~all(triple_ibd(a:b)) )
                        error = 1;
                        disp('error in segmentation');
                        return;
                    end
                end
            end
            
            temp_viewpair(1:nmarkers) = -1;
            for t = 1:nr
                a = regions(t,1);
                b = regions(t,2);
                if( any(triple_ibd(a:b)) )
                    if( all(triple_ibd(a:b)) )
                        numhetero = nnz(pivot_hetero(a:b));
                        numidentical = nnz(pair_identical(a:b));
                        if( numhetero > 0 )
                            % if not enough heterozygous marker at pivot
                            % not enough confidence to call
                            if( numhetero >= 10 )
                                temp_viewpair(a:b) = numidentical/numhetero;
                            end
                        end
                    else
                        error = 1;
                        disp('error in segmentation');
                        return;
                    end
                end
            end
            triple_viewpair(triple_index, 1:nmarkers) = temp_viewpair(1:nmarkers);           
        end       
    end
end

for i = 1:nseg
    
    a = intervals(i,1);
    b = intervals(i,2);
    
    temp_triplesview = zeros(num_triples, 1);
    value = mean(triple_viewpair(1:num_triples,a:b),2);
    temp_triplesview(value > 0.3) = 1;
    temp_triplesview(value >= 0 & value <= 0.3) = -1;
    temp_triplesview(value == -1) = 0;
    
    all_triplesview = zeros(ngeno, ngeno, ngeno);

    all_triplesview(ambig_triples > 0) = temp_triplesview(ambig_triples(ambig_triples > 0));
    triples_seg.triple_viewpair(i, 1:ngeno, 1:ngeno, 1:ngeno) = all_triplesview;
end

if( debug_mode == 1 )
    display(['          triple smoothing costs time: ', num2str(cputime - time), ' seconds']);  
end

end
















