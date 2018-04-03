function plot_posterior_wrap(input, output, parameters, f, id1, id2, limit)

genotyped = input.family_range{f}.family_range;
family = input.family_range{f}.structure;
temp = family(genotyped,8);
in1 = find(temp == id1);
in2 = find(temp == id2);
if( length(in1) ~= 1 || length(in2) ~= 1 || id1 == id2 )
    return;
end

s_id = input.family_range{f}.reverse_pairs(in1, in2);

coor = parameters.sampled_markerlist(:,1:2);
if( length(limit) ~= 2 )
    limit = [1,length(parameters.sampled_markerlist(:,1))];
else
    if( limit(1) > limit(2) || limit(1) < 1 || limit(2) > length(parameters.sampled_markerlist(:,1)) )
        return;
    end
end
option = 1;



if( option == 0 )
    posterior = [];
    posteriorIBD1 = [];
else
    posterior = squeeze(output.paths.posterior{f}(s_id,:,:));
    posteriorIBD1 = squeeze(output.paths.posteriorIBD1{f}(s_id,:,:,:));
end
viterbi = squeeze(output.paths.viterbi{f}(s_id,:,:));

print_id1 = input.family_range{f}.family_range(in1);
print_id2 = input.family_range{f}.family_range(in2);
print_id1 = id1;
print_id2 = id2;

genotype = input.family_genotype{f};
r1 = genotyped(in1);
r2 = genotyped(in2);
fa1 = family(r1,3);
ma1 = family(r1,4);
fa2 = family(r2,3);
ma2 = family(r2,4);
nmarkers = length(coor(:,1));
geno(1:4, 1:nmarkers, 1:2) = 0;
geno(1, 1:nmarkers, 1:2) = genotype(1:nmarkers, r1, 1:2);
geno(2, 1:nmarkers, 1:2) = genotype(1:nmarkers, r2, 1:2);
geno(3, 1:nmarkers, 1:2) = 0;
geno(4, 1:nmarkers, 1:2) = 0;
% plot parents if individuals are full siblings
if( fa1 == fa2 && fa1 ~= 0 && ma1 == ma2 && ma1 ~= 0 )
    geno(3, 1:nmarkers, 1:2) = genotype(1:nmarkers, fa1, 1:2);
    geno(4, 1:nmarkers, 1:2) = genotype(1:nmarkers, fa2, 1:2);
end
[oblist emission_option error] = generate_oblist(genotype, r1, r2, parameters.sampled_markerlist, 2);

plot_posterior(geno, oblist, posterior, posteriorIBD1, viterbi, coor, print_id1, print_id2, limit, option);

end



%  generate the observation list on the fly
%     sampled_markerlist = parameters.sampled_markerlist;
%     family_genotype = input.family_genotype{f};
%     input_pair = input.family_range{f}.pairs(s_id,1:2);
%     [oblist emission_option error] = generate_oblist(family_genotype, input_pair(1), input_pair(2), sampled_markerlist);

