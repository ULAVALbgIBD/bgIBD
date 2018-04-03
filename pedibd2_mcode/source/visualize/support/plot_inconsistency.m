function plot_inconsistency(merged_assignment, consistency, mendel_error, markerlist, genotype)

% plot inconsistency rate for each segment, merged segment

% plot rate of non-missing loci for each merged segment

global disp_mode;
disp_mode = 1;
if( ~isempty(markerlist) )
    if( disp_mode == 0 ) 
        coor = markerlist(:,1);
    else
        coor = markerlist(:,2);
    end
else
    coor = 1:length(consistency);
end

bit = consistency < 17;

intervals = merged_assignment.intervals;

[len c] = size(intervals);

figure;
for i = 1:len
    a = intervals(i,1);
    b = intervals(i,2);
    pro = nnz(bit(a:b))./(b-a+1);
    line(coor(intervals(i,1:2)), [pro,pro], 'linewidth', 3);
end

title('ratio of inconsistency')
ylim([0,1]);


if( ~isempty(genotype) )
    
    figure;
    hold on;
    
    bit_int = genotype(:,:,1) ~= genotype(:,:,2);
    count = sum(bit_int, 2);
    bit_non = genotype(:,:,1) ~= 0 & genotype(:,:,2) ~= 0;
    count_non = sum(bit_non, 2);
    
    count = count./count_non;
    count(isnan(count)) = 0;    %   heterozyosity rate
    count = count_non;          %   non-missing rate
    
    % hete/non-missing rate in inconsistent loci
    for i = 1:len
        a = intervals(i,1);
        b = intervals(i,2);
        disp_bit = bit;
        disp_bit(1:a-1) = false;
        disp_bit(b+1:end) = false;
        pro = mean(count(disp_bit));
        wid = 1 + (10*nnz(disp_bit)/(b-a+1));
        line(coor(intervals(i,1:2)), [pro,pro], 'linewidth', wid, 'color', 'r');
    end
    % hete/non-missing rate in consistent loci
    for i = 1:len
        a = intervals(i,1);
        b = intervals(i,2);
        disp_bit = ~bit;
        disp_bit(1:a-1) = false;
        disp_bit(b+1:end) = false;
        pro = mean(count(disp_bit));
        wid = 1 + (10*nnz(disp_bit)/(b-a+1));
        line(coor(intervals(i,1:2)), [pro,pro], 'linewidth', wid, 'color', 'b');
    end
    
end

end
