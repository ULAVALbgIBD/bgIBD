




% go to the directory of the program


% move to source code directory
cur_directory = pwd;
cd([pwd, '/source/paths']);

% compile c sources
mex 'cForward.c';
mex 'cBackward.c';
mex 'cViterbi.c';


% add souce code paths
% and all its sub-folders
addpath([cur_directory, '/source/paths']);

% move to input file directory
input_directory = [cur_directory, '/sample_input'];
cd(input_directory);

% run the haplotyping function

PedIBD2('pedigree_file.txt', 'genotype_file.txt', 'marker_file.txt');




