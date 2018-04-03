function ld = get_ld_locus(data, member_list, locus1, locus2)

seq(1:length(member_list),1:4) = 0;

seq(:,1:2) = data(member_list, 6 + locus1 * 2 - 1 : 6 + locus1 * 2);
seq(:,3:4) = data(member_list, 6 + locus2 * 2 - 1 : 6 + locus2 * 2);

ld = seq;

