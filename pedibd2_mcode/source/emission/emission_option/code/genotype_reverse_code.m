function reverse_pair = genotype_reverse_code(forward_pair)



reverse_pair(1:max(max(forward_pair)), 1:2) = 0;

for i = 1:length(forward_pair(:,1))
    for j = 1:length(forward_pair(i,:))
        if( forward_pair(i,j) <= 0 )
            disp('genotype_pair coding not compact');
        else
            reverse_pair(forward_pair(i,j),1:2) = [i,j];
        end
    end
end



end