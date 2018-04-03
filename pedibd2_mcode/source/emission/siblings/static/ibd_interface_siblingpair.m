function allele2ibd = ibd_interface_siblingpair()

allele2ibd = [];

for i1 = 1:2
    for i2 = 1:2
        for j1 = 1:2
            for j2 = 1:2
                allele2ibd(i1,i2,j1,j2) = ibd_map_siblingpair(i1,i2,j1,j2);
            end
        end
    end
end



end