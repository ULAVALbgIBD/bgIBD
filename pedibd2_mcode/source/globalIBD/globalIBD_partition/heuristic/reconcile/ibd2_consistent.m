function cor_pairs = ibd2_consistent(pairs)

ibd1 = pairs == 1;
ibd2 = pairs == 2;

cor_pairs = 0;



% check ibd2 transistivity
for i = 1:length(ibd2(:,1))
    for j = i+1:length(ibd2(i,:))
        if( ibd2(i,j) == 1 )
            if(nnz(xor(ibd2(i,:),ibd2(j,:))) > 0)
                nnz(pairs(i,:)~=pairs(j,:));
                cor_pairs = cor_pairs + 1;
%                 i
%                 pairs(i,:)
%                 j
%                 pairs(j,:)
            end
        end
    end
end


% check ibd1 consistency within ibd2 groups
for i = 1:length(ibd2(:,1))
    for j = i+1:length(ibd2(i,:))
        if( ibd2(i,j) == 1 )
            if(nnz(xor(ibd1(i,:),ibd1(j,:))) > 0)
               cor_pairs = cor_pairs + 1;
%                 i
%                 ibd1(i,:)
%                 j
%                 ibd1(j,:)
            end
        end
    end
end

end