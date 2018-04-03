function de_ibd_s = ibd_expected_sc(ind, epl_sc)



% 0: ibs0
% 1: ibs1
% 2: ibs2
% 3: missing genotype

g11 = 1;
g12 = 2;
g22 = 3;
g00 = 4;

ibd0 = 1;
ibd1 = 2;

Nmarkers = length(ind);
ce(1:2,1:4) = 0; 
% cummulated emmission

for i = 1:Nmarkers/2
    ce = ce + epl_sc{i};
end

ce = ce';

N(1:4) = 0;

count = 0;
t = 0; %consecutive homozygosity count
l(1:Nmarkers) = 0;

for i = 1:2:Nmarkers
    if( ind(i) == 0 || ind(i + 1) == 0 )
        N(g00) = N(g00) + 1;
    end
    if( ind(i) ~= ind(i+1) )
        N(g12) = N(g12) + 1;
    end 
    if( ind(i) == ind(i+1) )
        if( t == 0 )
            count = count + 1;
            t = 1;
        end
        l(count) = l(count) + 1;
        if( ind(i) == 1 )
            N(g11) = N(g11) + 1;
        else
            N(g22) = N(g22) + 1;
        end
    else
        t = 0;
    end
end


N = N';

de_ibd_s = ((ce'*ce)^(-1))*(ce'*N);


