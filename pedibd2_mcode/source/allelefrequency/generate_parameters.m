
function [parameters error] = generate_parameters(founder_source, map, genotype)


error = 0;
parameters = [];

physical_map = map.physical_map;
scaled_map = map.scaled_map;

[nrows, ncols] = size(scaled_map);
if( nrows <= 0 || ncols ~= 2 )
    error = 1;
    disp('error in processing marker list');
    return;
end

num_markers = nrows;

[nrows, ncols] = size(physical_map);
if( nrows <= 0 || ncols ~= 2 || nrows ~= num_markers )
    error = 1;
    disp('error in processing marker list');
    return;
end

% select all markers
chr_code = unique(scaled_map(1:num_markers,1));
if( length(chr_code) ~= 1 )
    error = 1;
    disp('error in chromsome number');
    return;
end

if( any(physical_map(1:num_markers,1) ~= chr_code) )
    error = 1;
    disp('error in chromosome number');
    return;
end

[nrows, nind, c] = size( genotype );

if( num_markers ~= nrows || nind <= 0 || c ~= 2 )
    error = 1;
    disp('genotype_file and marker_file number of SNPs do not match');
    return;
end

if( nind ~= length(founder_source) )
    error = 1;
    disp('error in processing genotype data');
    return;
end

markerlist(1:num_markers,1) = 1:num_markers;
markerlist(1:num_markers,2) = physical_map(1:num_markers,2);

parameters.chr = chr_code;
parameters.physical_map = markerlist;
if( isfield(map, 'genetic_map') )
    parameters.genetic_map = map.genetic_map;
end


[parameters.allele_freq error] = allele_frq(genotype, founder_source, num_markers);
if( error ~= 0 )
    disp('error in generating allele frequency');
    return;
end


sampled_markerlist = parameters.physical_map(1:1:num_markers, 1:2);
sampled_allele_count = parameters.allele_freq.all_af(sampled_markerlist(1:num_markers,1), 1:3);

sampled_af = minor_allele_frq(sampled_allele_count, 50);

parameters.sampled_markerlist = sampled_markerlist;
parameters.sampled_af = sampled_af;

parameters.typing_error = 0.001;
parameters.missing_rate = 0.01;

end















