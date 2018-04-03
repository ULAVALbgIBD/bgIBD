function p_value = single_logistic_regression(Y,X)

for i = 1:length(X(:,1))
    i
    %[B,dev,stats] = mnrfit([X(i,:);X(i,:)>=2]',Y);
    [B,dev,stats] = mnrfit(X(i,:),Y);
    p_value(i) = stats.p(2);
end

end