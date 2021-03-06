function plot_1pair(pairs, range)

figure;
if( isempty(range) )
    return;
end
genotyped = range.family_range;
nind = length(genotyped);
if( nind <= 0 )
    return;
end

family = range.structure;
[r, c] = size(family);
if( r < nind || c < 12 )
    return;
end

printid = family(genotyped,2);

[r, c] = size(pairs);
if( r <= 0 || r ~= c || r~= nind )
    return;
end


ibd12 = (pairs >= 1);
[index, neighbor12] = biggest_degree_node(ibd12, 1:nind, 1:nind);

neighbor12 = setdiff(neighbor12, find(pairs(index,:)==2));

neighborhood = ibd12(neighbor12,neighbor12);
subplot(1,2,1);
hold on;
spy(neighborhood, 'o');
spy(neighborhood==0, 'x');
set(gca,'XTick',1:length(neighbor12));
set(gca,'XTickLabel',num2str(printid(neighbor12)));
set(gca,'YTick',1:length(neighbor12));
set(gca,'YTickLabel',num2str(printid(neighbor12)));
xlabel('');

neighborhood = neighborhood * 2 - 1;
[V, D] = eig(neighborhood);
[d, i] = max(diag(D));
vector = V(:,i);
vector = vector./max(abs(vector));


subplot(1,2,2);
hold on;
for i = 1:length(vector)
    if( vector(i) > 0 )
        text(i,vector(i) + 0.1, num2str(printid(neighbor12(i))),'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');
    else
        text(i,vector(i) - 0.1, num2str(printid(neighbor12(i))),'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');
    end
end
bar(vector);
ylim([-1.1,1.1]);
axis off;

% set(gca,'XTick',1:length(vector));
% set(gca,'XTickLabel',num2str(printid(neighbor12)));

end












