function hetero = get_hetero_locus(data, locus, member_list)

hetero(1:3) = 0;
for i = 1:length(member_list)
    seq = data(member_list(i), 6 + locus * 2 - 1 : 6 + locus * 2);
    if( seq(1) ~= 0 && seq(2) ~= 0 )
        if( seq(1) == seq(2) )
            hetero(1) = hetero(1) + 1;
        else
            hetero(2) = hetero(2) + 1;
        end
    else
        hetero(3) = hetero(3) + 1;
    end
end

