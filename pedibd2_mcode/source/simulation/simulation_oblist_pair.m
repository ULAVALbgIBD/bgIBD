function [ random_ibd_result, random_ibs_result ] = simulation_oblist_pair( s_id, input_pairlist, pedigree_all_missing, random_ibd, random_genotype )

temp = input_pairlist(s_id, 1:2);
sp = pedigree_all_missing(temp,2);

% sp is the inner id

random_ibd_result = ibs(random_ibd(sp(1), 1:end), random_ibd(sp(2), 1:end));

% for i = 1:length(sampled_markerlist(:,1))
%     temp(i) = random_ibd(sp(2), (i-1)*2 + 2);
% end

% plot(temp);

% plot(random_ibd_result);
% hold on;

random_ibs_result = ibs(random_genotype(sp(1), 1:end), random_genotype(sp(2), 1:end));

end

