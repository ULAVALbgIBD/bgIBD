function plot_relationship(member_list, pedigree, relationship)

interperson = 30;
interallele = 8;
intergeneration = -10;
allelesize = 5;



hold on;
for i = 1:length(member_list)
    plot(i*interperson, pedigree(member_list(i),9) * intergeneration, 'o', 'MarkerSize', allelesize);
    plot(i*interperson+interallele, pedigree(member_list(i),9) * intergeneration, 'o', 'MarkerSize', allelesize);
    if( i <= length(relationship) )
        text(i*interperson,pedigree(member_list(i),9) * intergeneration - 1 ,['(', num2str(relationship(i,1)), ',', num2str(relationship(i,2)), ')']);
    end
 
    if( pedigree(member_list(i),6) == 2 )
        %text(i*interperson,pedigree(member_list(i),9) * intergeneration + 1 , num2str(pedigree(member_list(i),8)), 'color', 'red');
        text(i*interperson,pedigree(member_list(i),9) * intergeneration + 1 , num2str(pedigree(member_list(i),2)), 'color', 'red');
    else
        %text(i*interperson,pedigree(member_list(i),9) * intergeneration + 1 , num2str(pedigree(member_list(i),8)), 'color', 'blue');
        text(i*interperson,pedigree(member_list(i),9) * intergeneration + 1 , num2str(pedigree(member_list(i),2)), 'color', 'blue');
    end
end

for i = 1:length(member_list)
    
    this = member_list(i);
    f = pedigree(this,1);
    temp1 = find(pedigree(:,1) == f);
 
    for j = 3:4
        inner_parent = pedigree(this, j);
        if( inner_parent ~= 0 )
            parent = temp1(inner_parent);
            temp = find( member_list == parent );
            if( ~isempty(temp) )
                if( relationship(i,1) == 0 || relationship(temp,1) == 0 )
                    continue;
                end                
                for m = 1:2
                    if( m == 2 && relationship(i,1) == relationship(i,2) )
                        line([i*interperson,i*interperson+interallele*(1)], [pedigree(member_list(i),9) * intergeneration, pedigree(member_list(i),9) * intergeneration]);
                        
                    end
                    for n = 1:2
                        if( ( j == 3 && m == 1 ) || ( j == 4 && m == 2 ) )
                            if( n == 2 && relationship(temp,1) == relationship(temp,2) )
                                line([temp*interperson,temp*interperson+interallele*(1)], [pedigree(member_list(temp),9) * intergeneration, pedigree(member_list(temp),9) * intergeneration]);
                                continue;
                            end                        
                            if( relationship(i,m) == relationship(temp,n) )
                                line([i*interperson+interallele*(m-1),temp*interperson+interallele*(n-1)], [pedigree(member_list(i),9) * intergeneration, pedigree(member_list(temp),9) * intergeneration]);
                            end
                        end
                    end
                end
            end
        end
    end
end

axis off;

hold off;

end