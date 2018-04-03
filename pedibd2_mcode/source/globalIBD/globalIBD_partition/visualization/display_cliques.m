
function display_cliques(clique_assignment, pairs, range)

list = [];

for i = 1:length(clique_assignment(:,1))-1
    temp = find(clique_assignment(i,:) == 1);
    temp = setdiff(temp, list);
    list = [list, temp];
end

inferred_pairs = clique2pairs(clique_assignment);

nnz(inferred_pairs-pairs);
num = length(list);
pairs = pairs(list,list);
inferred_pairs = inferred_pairs(list,list);

ibd1 = [];
ibd2 = [];
difference = [];

for i = 1:num
    for j = 1:num
        if( inferred_pairs(i,j) == 1 )
            ibd1 = [ibd1;[i,j]];
        end
        if( inferred_pairs(i,j) == 2 )
            ibd2 = [ibd2;[i,j]];
        end
        if( inferred_pairs(i,j) ~= pairs(i,j) )
            difference = [difference;[i,j]];
        end
    end
end



hold on;
plot(ibd1(:,1), ibd1(:,2), '+', 'MarkerSize', 5);
plot(ibd2(:,1), ibd2(:,2), 'o', 'MarkerSize', 5);
if( ~isempty(difference) )
    plot(difference(:,1), difference(:,2), 'x');
end
hold off;

set(gca,'XTick',1:num);
set(gca,'XTickLabel',(num2cell(range(list))));
set(gca,'YTick',1:num);
set(gca,'YTickLabel',(num2cell(range((list)))));



end