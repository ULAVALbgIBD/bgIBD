function freq = get_freq_locus(data_locus, member_weight)

freq(1:3) = 0;
n_ind = length(member_weight);
[n_rows, n_cols] = size(data_locus);
if( n_rows ~= n_ind || n_cols ~= 2 )
    disp('error in generating allele frequency');
end
for i = 1:length(member_weight)
    if( member_weight(i) ~= 0 )
        seq = data_locus(i, 1 : 2);
        freq = freq + member_weight(i) .* get_freq(seq);
    end
end

