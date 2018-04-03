function [marker_dist] = generate_marker_dist(width, len, num_marker)


% dist is the distance between marker centers
% distance should be at least 2* marker size
% dist = ratio * marker_size;


[x_num y_num] = generate_layout(width, len, num_marker); 

w_res = width/(x_num - 1 + 0.5 + 0.5);
l_res = len/(y_num - 1 + 0.5 + 0.5);

% leave half distance on each size of edges

marker_dist = min(w_res, l_res);

ratio = 2;
marker_size = marker_dist/ratio;

end