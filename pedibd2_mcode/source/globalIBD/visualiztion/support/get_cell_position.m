function [x, y] = get_cell_position(width, len, num_marker, order, dist)

% markersize and distance between markers


x = 0;
y = 0;


if( order > num_marker )
    return;
end

% number of markers per row and per column
[x_num y_num] = generate_layout(width, len, num_marker);

row = 1 + floor((order-1)/x_num);
col = 1 + mod((order-1),x_num);

% distance from the first marker center
dist_x = (row-1) * (dist);
dist_y = (col-1) * (dist);

% check whether the cell can hold all markers
% leave half distance margin on each end
min_width = (x_num - 1 + 0.5 + 0.5) * dist;
min_len = (y_num - 1 + 0.5 + 0.5) * dist;


if( min_width > width || min_len > len )
    disp('cell can not hold all markers');
    return;
end

% coordinates with respect to the cell center

% first marker to last marker, center to center
total_dist_x = (dist) * (x_num - 1);
total_dist_y = (dist) * (y_num - 1);

% x should col order
% y use row order
x = - total_dist_x/2 + (col-1) * (dist);
y = - total_dist_y/2 + (row-1) * (dist);



end
