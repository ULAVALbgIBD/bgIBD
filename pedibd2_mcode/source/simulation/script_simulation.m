
%%

input_recfrac(1:length(parameters.sampled_markerlist(:,2))-1) = 0;
for i = 2:length(parameters.sampled_markerlist(:,2))
    temp = parameters.sampled_markerlist(i,2) - parameters.sampled_markerlist(i-1,2);
    temp = temp/10^6;
    temp = temp/10^2;
    input_recfrac(i-1) = [1 + exp(-2*temp)]/2; %probability of even recombinations
end

simulation.recfrac = input_recfrac;

clear i temp input_recfrac chr;


%%
% families
% [1,6,11,15,18,23]
% family 6 both parents
% family 7 father
% 

%%%%% ucsd data !!!

% see UCSD script for more details of how to generate these data

simulation.founder_source = allFOUNDERsource > 0 & allERROR < 0.02;
simulation.haplotypes = allHAPLOTYPE;
%%%%%


%%

% generate allele frequency to feed to merlin

simulation.allele_freq = ...
    allele_frq(simulation.haplotypes, ...
    simulation.founder_source, ...
    size(parameters.sampled_markerlist,1));

%%

