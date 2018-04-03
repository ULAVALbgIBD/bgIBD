function markersize_pix = view_mendelian_map(family_display, phase, total_pixel_width)

    % markersize is fixed in pixels, fonts are scalable

    % the axes is within boundary [0,num+1]
    % font size of alleles is scaled inversely to num
    % title size is scaled to 1/20
    % marker size is scaled inversely to num * maximum individuals in a cell
    
    if( isempty(family_display) )
        
    end
 
    assignment = family_display.alleles;
    nind = family_display.num_ind;
    num = length(family_display.unique_alleles);
    
    %%%% generate positions of each individual in each cell
    
    if( phase == 0 )
        for i = 1:nind
            p = assignment(i,1);
            m = assignment(i,2);
            % make sure 2 is smaller than 1
            if( abs(p) > abs(m) )
                assignment(i,1:2) = [p,m];
            end
            if( abs(m) > abs(p) )
                assignment(i,1:2) = [m,p];
            end
            if( abs(m) == abs(p) )
                if( m < 0 )
                    assignment(i,1:2) = [p,m];
                else
                    assignment(i,1:2) = [m,p];
                end
            end
        end
    end
    [x x_offset y y_offset ratio] = generate_individual_positions(assignment, family_display.unique_alleles);
    % ratio is the marker distance with respect to cell width
    
        
    center(1:num) = 1:num;
    x_pos = center(x) + x_offset;
    y_pos = center(y) + y_offset;
    x_pos(nind+1:nind*2) = y_pos(1:nind);
    y_pos(nind+1:nind*2) = x_pos(1:nind);
    % expand to fill symmetric cells
    
    %%%%% ratio is normalized proportion of a cell
    % pixel width of each cell, total_pixel_width/(2*num)
    
    hold on;
    threshold = 20; % threshold for smallest cell numbers range for good visual;
    cellnum = num + 1; % actual cell numbers, including marginal blank cells;
    if( cellnum < threshold )
        cellnum = threshold;
    end
    cellsize_pix = 0.5*total_pixel_width/cellnum;
    markersize_pix = ratio * 0.5 * cellsize_pix;
    

    %%%% put individuals into cells
    
    
    plot_individuals(x_pos, y_pos, markersize_pix, family_display, 1);

    
    %%%% draw cell dividers
    
    if( phase == 1 )
        for i = 1:num+1
            line([i-0.5,i-0.5], [-0.5,num+0.5+1]);
            line([-0.5,num+0.5+1], [i-0.5,i-0.5]);
        end
    else
        for i = 1:num+1
            line([i-0.5,i-0.5], [0.5,i+0.5]);
            line([i-1-0.5,num+0.5], [i-0.5,i-0.5]);
        end        
    end
    
    % draw allele names
    allele_location_x(1:num) = 1:num;
    allele_location_y(1:num) = 0;
    allele_location_x(1:num) = num+1;
    allele_location_y(1:num) = 1:num;
    
    
    fontsize = 0.5/cellnum;

    for i = 1:num        
        text(i,0,family_display.allele_name{i}, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontUnits', 'normalized', 'FontSize', fontsize * 0.9, 'Rotation', 37.5);
        text(num+1,i,family_display.allele_name{i}, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right', 'FontUnits', 'normalized', 'FontSize', fontsize * 0.9, 'Rotation', -37.5);
        text(i,i,family_display.allele_name{i}, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'FontUnits', 'normalized', 'FontSize', fontsize * 0.7, 'Rotation', 0);
    end
    
    % Mendelian cross
    scatter(num+1, 0, cellsize_pix^2, 'X');
    
    title_size = 0.02;
    title('ALLELE MAP', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', title_size, 'FontUnits', 'normalized');
    
    % change to 45 angle view
    view(225,90);
    
    
    % filled cell boundaries
    leftmost = -0.5;
    lowermost = -0.5;
    rightmost = num+1+0.5;
    uppermost = num+1+0.5;
    
    % wrapping blank cell boundaries

    x_true_center = (rightmost - leftmost)*2/3 + leftmost;
    x_target_center = (cellnum+1+0.5 + 0.5) *2/3 + -0.5;
    x_difference = x_true_center - x_target_center;
    y_true_center = (uppermost - lowermost)*1/3 + lowermost;
    y_target_center = (cellnum+1+0.5 + 0.5) *1/3 + -0.5;
    y_difference = y_true_center - y_target_center;
    leftmost = leftmost + x_difference;
    rightmost = leftmost + cellnum + 2;
    lowermost = lowermost + y_difference;
    uppermost = lowermost + cellnum + 2;

    % move view center
    xlim([leftmost,rightmost]);
    ylim([lowermost,uppermost]);
    
    % set axis off, remove ticks
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    axis off;
    set(gcf,'Color', 'w');
    hold off;
    
    
end













