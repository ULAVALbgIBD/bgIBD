function p_value = categorical_null(emission, observation, ibd, emission_option)

global temp1;
global temp2;
if( isempty(temp2) )
    temp2 = 0;
end


% multinomial test
% approximate it using chi-square test

len1 = length(emission);
len2 = length(observation);
len3 = length(ibd);
p_value = 1;

if( len1 ~= len2 || len2 ~= len3 )
    disp('error in check data consistency');
    return;
end

if( len1 <= 0 )
    disp('no genotype data for evaluating relationship');
    return;
end


len = len1;

ob_missing_code = emission_option.ob_missing_code;
pair2ob = emission_option.pair2ob;
ob_states = unique(pair2ob);

expected(ob_states) = 0;
occurence(ob_states) = 0;


for i = 1:len
    ob = observation(i);
    hid = ibd(i);
    em = emission{i};
    expected(ob_states) = expected(ob_states) + em(hid,ob_states);
    occurence(ob) = occurence(ob) + 1;
end

temp2 = temp2 + 1;
temp1(temp2, 1:length(ob_states)) = occurence;
temp1(temp2, length(ob_states)+1:2*length(ob_states)) = expected;


%do not count missing
temp = occurence(ob_missing_code);
occurence = occurence .* (len/(len-temp));
temp = expected(ob_missing_code);
expected = expected .* (len/(len-temp));



chisquare = 0;
for i = 1:length(ob_states)
    ob = ob_states(i);
    if( ob ~= ob_missing_code )
        chisquare = chisquare + ((occurence(ob)-expected(ob))^2)/expected(ob);
    end
end

df = length(ob_states) - 1 - 1;
p_value = 1 - chi2cdf(chisquare, df);

% allele_frequency = 0;
% for j = 1:100
%     p = j/200;
%     if( ibd(1) == 2 )
%         expected(1) = 0;
%         expected(2) = 2*p*(1-p);
%         expected(3) = p^3 + (1-p)^3;
%         expected(4) = p*(1-p)^2 + (1-p)*p^2;
%     else
%         expected(1) = 2 * p^2 * (1-p)^2;
%         expected(2) = p^2*p*(1-p)*4 + (1-p)^2*(1-p)*p*4;
%         expected(3) = p^4 + (1-p)^4;
%         expected(4) = (2*p*(1-p))^2;
%     end
%     expected(1:4) = expected(1:4) * len;
%     for i = 1:length(ob_states)
%         ob = ob_states(i);
%         if( ob ~= ob_missing_code )
%             chisquare = chisquare + ((occurence(ob)-expected(ob))^2)/expected(ob);
%         end
%     end
%     df = length(ob_states) - 1 - 1;
%     p_value_temp = 1 - chi2cdf(chisquare, df);   
%     if( p_value_temp > p_value )
%         p_value = p_value_temp;
%         allele_frequency = p;
%     end
% end
% 
% disp(['suggested allele frequency', num2str(allele_frequency)]);

temp1(temp2, 2*length(ob_states)+1) = p_value;
temp1(temp2, 2*length(ob_states)+2) = ibd(1);



end






