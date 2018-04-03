function [ p t ix ] = combine_fb( fp_i, ft_i, bp_i, bt_i )

% foreground and background probability
s = length(fp_i);



fnz = find( fp_i ~= 0 );
bnz = find( bp_i ~= 0 );

ix(1:length(fnz)*length(bnz), 1:2) = 0;

fp_log = log(fp_i);
ft_log = log(ft_i);
bp_log = log(bp_i);
bt_log = log(bt_i);

p(1:length(fnz)*length(bnz)) = 0;
t(1:length(fnz)*length(bnz), 1:length(fnz)*length(bnz)) = 0;

for i1 = 1:length(fnz)
    for j1 = 1:length(bnz)
        code1 = (i1-1) * s + j1-1 + 1;
        p(code1) = fp_log(fnz(i1)) + bp_log(bnz(j1));
        ix(code1,1:2) = [fnz(i1), bnz(j1)];
        ix(code1,3) = (fnz(i1)-1) * s + bnz(j1)-1 + 1;
        for i2 = 1:length(fnz)
            for j2 = 1:length(bnz)
               code2 = (i2-1) * s + j2-1 + 1;
               t(code1,code2) = ft_log(fnz(i1),fnz(i2)) + bt_log(bnz(j1),bnz(j2));
            end
        end
    end
end

end

