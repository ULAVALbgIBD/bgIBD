function [margin_prob, tran_matrix] = transition_identity(input_ip, markerlist)


% change accumulated in one column or one row
% ip is inheritance path

[fp, ft] = transition_foreground(input_ip, markerlist);
[bp, bt] = transition_background(input_ip, markerlist);


for i = 1:length(markerlist(:,1))
    [margin_prob, tran_matrix{i}] = combine_fb(fp, ft{i}, bp, bt{i});
end

end










