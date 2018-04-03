function [ x x_offset y y_offset ratio] = generate_individual_positions( assignment, alleles )



    
    num = length(alleles);
    n_ind = length(assignment(:,1));

    cell_width = 1;
    % use normalized scale
    
    
    x(1:n_ind) = 0;
    y(1:n_ind) = 0;
    x_offset(1:n_ind) = 0;
    y_offset(1:n_ind) = 0;
    

    if( n_ind <= 0 )
        return;
    end
    
    % caculate the occurence of each allele in the family
    % get the cell identity of each individual    
    count(1:num,1:num) = 0;
    order(1:n_ind) = 0;
    for i = 1:n_ind
        p = assignment(i,1);
        q = assignment(i,2);
        x(i) = find(alleles == p);
        y(i) = find(alleles == q);

        count(x(i),y(i)) = count(x(i),y(i)) + 1;
        count(y(i),x(i)) = count(y(i),x(i)) + 1;
        order(i) = count(x(i),y(i));
    end
    

    % use square cells
    width = cell_width;
    len = cell_width;
    
    max_volume = max(max(count));
    % maximum number of markers in a cell
    % generate marker distance based on the fullest cell
    [dist] = generate_marker_dist(width, len, max_volume);
    ratio = dist/cell_width;

    % get the position of each individual in a cell
    for i = 1:n_ind        
        [x_offset(i), y_offset(i)] = get_cell_position(width, len, count(x(i),y(i)), order(i), dist);
    end


end

