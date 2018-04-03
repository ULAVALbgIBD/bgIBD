function [pair2ob missing] = pair2obMap(eo, option)


% option 1: genotype
% option 2: ibs
% option 3: ibs, differatiate ibd2



genotype2pair = eo.genotype2pair;
allele2genotype = eo.allele2genotype;
genotype2pair = eo.genotype2pair;
g0 = eo.pair_missing_code;

pair2ob(1:max(max(genotype2pair))) = 0;
% initial assignment is valid

g11 = allele2genotype(1,1);
g12 = allele2genotype(1,2);
g21 = allele2genotype(2,1);
g22 = allele2genotype(2,2);



if( option == 1 )
    
    pair2ob(1:max(max(genotype2pair))) = 1:max(max(genotype2pair));
    missing = g0;
    
end

if( option == 2 )


    ibs0 = 1;
    ibs1 = 2;
    ibs2 = 3;
    missing = 4;

    g11 = 1;
    g12 = 2;
    g22 = 3;


    pair2ob(genotype2pair(g11,g11)) = ibs2;
    pair2ob(genotype2pair(g11,g12)) = ibs1;
    pair2ob(genotype2pair(g11,g22)) = ibs0;
    pair2ob(genotype2pair(g12,g11)) = ibs1;
    pair2ob(genotype2pair(g12,g12)) = ibs2;
    pair2ob(genotype2pair(g12,g22)) = ibs1;
    pair2ob(genotype2pair(g22,g11)) = ibs0;
    pair2ob(genotype2pair(g22,g12)) = ibs1;
    pair2ob(genotype2pair(g22,g22)) = ibs2;

end

if( option == 3 )
    
    ibs0 = 1;
    ibs1 = 2;
    ibs2_homo = 3;
    ibs2_hete = 4;
    missing = 5;
    
    g11 = 1;
    g12 = 2;
    g22 = 3;

    pair2ob(genotype2pair(g11,g11)) = ibs2_homo;
    pair2ob(genotype2pair(g11,g12)) = ibs1;
    pair2ob(genotype2pair(g11,g22)) = ibs0;
    pair2ob(genotype2pair(g12,g11)) = ibs1;
    pair2ob(genotype2pair(g12,g12)) = ibs2_hete;
    pair2ob(genotype2pair(g12,g22)) = ibs1;
    pair2ob(genotype2pair(g22,g11)) = ibs0;
    pair2ob(genotype2pair(g22,g12)) = ibs1;
    pair2ob(genotype2pair(g22,g22)) = ibs2_homo;
       
end

pair2ob(g0) = missing;


end










