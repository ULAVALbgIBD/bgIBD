function [result, rows] = assign_clique(result, clique)

[rows cols] = size(result);
result(rows, clique) = 1;
result(rows+1,:) = 0;

end