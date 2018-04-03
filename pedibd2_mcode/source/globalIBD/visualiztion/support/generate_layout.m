
% number of individuals per row and per column to fit in num individuals

function [x_num y_num] = generate_layout(width, len, num_marker)

    x_num = ceil(sqrt(num_marker*(len/width)));
    y_num = ceil(num_marker/x_num);
%     y_num = ceil(sqrt(num_marker*(width/len)));

end