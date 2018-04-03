



function output_genF = descent_graph_2ind(input_pedigree, family_id, ind1, ind2)

    id1 = input_pedigree(ind1, 2);
    id2 = input_pedigree(ind2, 2);
    temp = find(input_pedigree(:,1) == family_id);
    family = input_pedigree(temp,3:4);
    
    parental(1:2*length(temp),1:2) = 0;
    ancestor(1:2*length(temp),1:2*length(temp)) = 0;
    related(1:2*length(temp),1:2*length(temp)) = 0;
    
    for i = 1:length(family(:,1))
        if( family(i,1) ~= 0 )
            %paternal allele source
            parental(i*2-1,1:2) = [family(i,1)*2-1, family(i,1)*2];
            ancestor(family(i,1)*2-1,i*2-1) = 1;
            ancestor(family(i,1)*2,i*2-1) = 1;
        else
            parental(i*2-1,1:2) = 0;
        end
        if( family(i,2) ~= 0 )
            %maternal allele source
            parental(i*2,1:2) = [family(i,2)*2-1, family(i,2)*2];
            ancestor(family(i,2)*2-1,i*2) = 1;
            ancestor(family(i,2)*2,i*2) = 1;
        else
            parental(i*2,1:2) = 0;
        end
    end
    
    for k = 1:length(ancestor)
        for i = 1:length(ancestor)
            for j = 1:length(ancestor)
                if( ancestor(i,j) == 1 )
                    ancestor(i,:) = ancestor(i,:) | ancestor(j,:);
                end
            end
        end
    end
    
    for i = 1:length(related)
        for j = 1:length(related)
            if( sum(ancestor(:,i) & ancestor(:,j)) > 0 || ancestor(i,j) )
                related(i,j) = 1;
            end
        end
    end
    
    vec = [id1*2-1, id1*2, id2*2-1, id2*2];
    
    genF = descent_graph( parental, ancestor, related, vec );
    
    output_genF(1:length(genF(:,1)),2:length(genF(1,:))+1) = genF;
    output_genF(:,1) = length(vec);
    
end









