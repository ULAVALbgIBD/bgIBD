function error = plot_inheritance(assignment, family, list, coor, limit)

error = 0;

nloc = length(coor);
if( nloc <= 0 )
    error = 1;
    return;
end

[nind, c] = size(family);
if( nind <= 0 || c < 12 )
    error = 1;
    return;
end

alleles = assignment.alleles_all;
intervals = assignment.intervals;

[nseg c] = size(intervals);
if( nseg <= 0 || c < 2 )
    error = 1;
    return;
end

if( any(intervals(1:nseg,1) > nloc) || any(intervals(1:nseg,2) > nloc) )
    error = 1;
    return;
end
newint = coor(intervals(1:nseg,1:2));

len = length(alleles);
if( len ~= nseg )
    error = 1;
    return;
end

for i = 1:nseg
    oneseg = alleles{i};
    [r c] = size(oneseg);
    if( r ~= nind || c ~= 2 )
        error = 1;
        return;
    end    
end

num = length(list);
if( num <= 0 )
    return;
end

if( any(list > nind) || any(list < 1) )
    error = 1;
    return;
end

categories = [];
for i = 1:nseg
    categories = union(categories, unique(alleles{i}(list(1:num),1:2)));
end

ncol = ceil((length(categories))^(1/3));

colmap = zeros(ncol^3, 3);

temp = 0;
for i = 1:ncol
    for j = 1:ncol
        for k = 1:ncol
            temp = temp + 1;
            colmap(temp, 1) = (k-1)/(ncol-1);
            colmap(temp, 2) = (j-1)/(ncol-1);
            colmap(temp, 3) = (i-1)/(ncol-1);
        end
    end
end

mid = mean(coor);

for i = 1:num
    
    text(mid, i*4 - 2 + 1, num2str(family(list(i),8)), 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'FontUnits', 'normalized', 'FontSize', 0.05);
    
    for j = 1:2
        y = i*4 - 2 + j;
        for k = 1:nseg
            
            if( ~isempty(limit) )
                if( intervals(k,2) < limit(1) || intervals(k,1) > limit(2) )
                    continue;
                end
            end
            
            curseg = alleles{k}(list(i),j);
            index = find(categories == curseg, 1);
            color = colmap(index, 1:3);

             
            line([newint(k,1), newint(k,2)], [y,y], 'Color', color, 'LineWidth', 4);
            text(mean(newint(k,1:2)), y+0.2, num2str(curseg), 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'FontUnits', 'normalized', 'FontSize', 0.02);
        end
    end

end

for i = 1:nseg
    text( mean(newint(i,1:2)), 0, num2str(i), 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');
end

ylim([0,num*4+4]);




















