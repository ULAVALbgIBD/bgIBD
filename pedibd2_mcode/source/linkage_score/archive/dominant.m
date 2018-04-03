function output = dominant(assignment, family_range, pedigree);

disease_status = pedigree(family_range.family_range_outer, 6);

for chr = 1:length(assignment)
    temp = assignment{chr}.alleles;
    result = [];
    for i = 1:length(temp)
        temp1 = [];
        for j = 1:length(disease_status)
            if( disease_status(j) == 1 )
                temp1 = union(temp1,temp{i}(j,1:2)); 
            end
        end
        count = 0;
        result(i,1:2) = assignment{chr}.intervals(i,1:2);
        result(i,3) = 1;
        temp3 = union(temp{i}(:,1),temp{i}(:,2));
        for j = 1:length(disease_status)
            if( disease_status(j) == 2 )
                count = count + 1;
                temp2{count} = [];
                temp2{count} = union(temp2{count}, temp{i}(j,1:2));
                temp3 = intersect(temp3, temp{i}(j,1:2));
                temp2{count} = setdiff(temp2{count}, temp1);
                temp3 = setdiff(temp3, temp1);
                if( isempty(temp2{count}) )
                    result(i,3) = 0;
                end
                if( isempty(temp3) )
                    result(i,3) = 0;
                end
            end
        end
    end
    output{chr}.intervals = result;
end


end

