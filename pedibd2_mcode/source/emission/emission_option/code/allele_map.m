function code = allele_map()


g11 = 1;
g12 = 2;
g21 = 2;
g22 = 3;

code(1:2,1:2) = 0;
code(1,1) = g11;
code(1,2) = g12;
code(2,1) = g21;
code(2,2) = g22;

end