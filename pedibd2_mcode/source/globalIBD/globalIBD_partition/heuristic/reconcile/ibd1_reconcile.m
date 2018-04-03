function cor_pairs = ibd1_reconcile(pairs_vit, pairs_pos, range)

ibd1 = pairs_vit == 1;
ibd2 = pairs_vit == 2;
ibd0 = pairs_vit == 0;
cor_pairs = pairs_vit;
cor_pairs(:,:) = 0;


% compare the posterior probability to remove errors


% remove ibd2 transistivity errors

% previously removed


% reconcile ibd1 relationships

neighbor(1:length(ibd2(:,1))) = 1:length(ibd2(:,1));

for i = 1:length(ibd2(:,1))
    for j = i+1:length(ibd2(i,:))
        if( ibd2(i,j) == 1 )
            neighbor(j) = neighbor(i);
        end
    end
end

% ibd2 errors are in general rare, treatment is not important

for i = 1:length(ibd2(:,1))
    for j = 1:length(ibd2(:,1))
        n1 = find(neighbor == neighbor(i));
        n2 = find(neighbor == neighbor(j));
        % compare ibd0 and ibd1 relationships across all ibd2 neighbors
        temp1 = pairs_pos(n1,n2,1);
        temp2 = pairs_pos(n1,n2,2);
        temp3 = pairs_pos(n1,n2,3);
        temp4(1) = prod(prod(temp1));
        temp4(2) = prod(prod(temp2));
        temp4(3) = prod(prod(temp3)); %ibd2 situation not considered again
        [x, ix] = max(temp4);
        cor_pairs(i,j) = ix-1;
        cor_pairs(i,j) = round(mean(mean(pairs_vit(n1,n2))));
        cor_pairs(i,j) = max(max(pairs_vit(n1,n2)));
        if( cor_pairs(i,j) ~= pairs_vit(i,j) )
%             range([n1,n2])
%             pairs_vit([n1,n2],[n1,n2])
%             cor_pairs([n1,n2],[n1,n2])
        end
    end
end

end