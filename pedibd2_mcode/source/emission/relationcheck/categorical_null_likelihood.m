function [p_value, incon, error] = categorical_null_likelihood(emission, observation, viterbi_foreground, emission_option, ibd_viterbi)

error = 0;
% multinomial test
% approximate it using chi-square test

%emission is only correct for line 1 and line 2
states = [
    1, 1, 1, 1; %ibd0
    2, 1, 1, 1; %ibd1
    3, 1, 1, 1;
    1, 2, 1, 1;
    1, 3, 1, 1;
    1, 1, 2, 1;
    1, 1, 3, 1;
    1, 1, 1, 2;
    1, 1, 1, 3;
    2, 1, 1, 2; %ibd2
    3, 1, 1, 2;
    2, 1, 1, 3;
    3, 1, 1, 3;
    1, 2, 2, 1;
    1, 3, 2, 1;
    1, 2, 3, 1;
    1, 3, 3, 1;
    ];


len1 = length(emission);
len2 = length(observation);
len3 = length(viterbi_foreground);
p_value = 1;
incon = [];

if( len1 ~= len2 || len2 ~= len3 )
    error = 1;
    disp('error in check data consistency');
    return;
end

if( len1 <= 0 )
    error = 1;
    disp('no genotype data for evaluating relationship');
    return;
end

num_ob = len1;

ob_missing_code = emission_option.ob_missing_code;
pair2ob = emission_option.pair2ob;
ob_states = unique(pair2ob);

if( any((viterbi_foreground>2)) )
    return;
    % only works for ibd0 and ibd1 check
end
viterbi_states = unique(ibd_viterbi.foreground_map);

likelihood_all(viterbi_states) = 0;
likelihood_target = 0;
ratio(viterbi_states) = 0;
ratio2max = 0;

ibs0 = 0;
% compare the likelihood of different hidden states
% only consider bare foreground states of ibd0 and ibd1
for i = 1:num_ob
    ob = observation(i);
    hid = viterbi_foreground(i);
    em(viterbi_states,1:length(emission{i}(1,:))) = 0;
    for j = 1:length(viterbi_states)
        em(viterbi_states(j),:) = mean(emission{i}(ibd_viterbi.foreground_map == viterbi_states(j) & ibd_viterbi.background_map == 1,:), 1);
        %pick one representative from emission probability
        %do not consider background emission
        %background ibd has a very small prio, not suitable for likelihood
    end
    likelihood_all(viterbi_states) = likelihood_all(viterbi_states) + log(em(:,ob))';
    likelihood_target = likelihood_target + log(em(hid,ob));
    ratio = ratio + (-2)*log(em(hid,ob)./em(:,ob))';
    ratio2max = ratio2max + (-2)*log(em(hid,ob)./max(em(:,ob)));
    if( ob == 1 )
        ibs0 = ibs0 + 1;
    end
end

% only test ibd1, viterbi_foreground
if( hid == 2 )
    incon = ibs0/num_ob;
    %disp(['estimated genotyping error ', num2str(4*(ibs0/num_ob))]);
else
    incon = [];
end

% output typing error 

% likelihood ratio test
df = length(ob_states) - 1;
p_value = 1 - chi2cdf(max(ratio), df);



end






