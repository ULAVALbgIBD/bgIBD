function pairs_cor = ibd2_reconcile(pairs_vit, pairs_pos)

ibd1 = pairs_vit == 1;
ibd2 = pairs_vit == 2;
ibd0 = pairs_vit == 0;
pairs_cor = pairs_vit;

num = length(ibd2(1,:));

% compare the posterior probability to remove errors


% remove ibd2 transistivity errors
% make all groups different
equivalence(1:num) = 1:num;

for i = 1:length( ibd2(:,1) )
    for j = 1:length( ibd2(i,:) )
        if( ibd2(i,j) == 1 )
            cor_pairs(i,j) = 2;
            if( equivalence(j) ~= equivalence(i) )
                for k = 1:num
                    if( equivalence(k) == equivalence(j) )
                        equivalence(k) = equivalence(i);
                    end
                end
            end
        end
    end
end


for i = 1:num
    temp = find(equivalence == i);
    if( ~isempty(temp) )
        pairs_cor(temp,temp) = 2;
    end
end

end