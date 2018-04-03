function [ output_fams ] = generate_random_genotype( ibd, data, founder )


[nLOC, nIND, ~] = size(ibd);

founder = find(founder);
fn = length(founder);
founder = founder(randperm(fn)); %randomize the order

output_fams = zeros(nLOC, nIND, 2);

countp = 0;
countm = 0;

% fill in founders first
for j = 1:nIND    
    p = 1;
    m = 2;
    if( ibd(1,j,1) == j )% is a founder
        % generate a random number between 1 to fn
        temp = floor(rand()*(fn-1)) + 1;
        temp = founder(temp);
        countp = countp + 1;
        temp = founder(countp);
        output_fams(1:nLOC,j,p) = data(1:nLOC,temp,p);
    end
    if( ibd(1,j,2) == -j )
        temp = floor(rand()*(fn-1)) + 1;
        temp = founder(temp);
        countm = countm + 1;
        temp = founder(countm);
        output_fams(1:nLOC,j,m) = data(1:nLOC,temp,m);
    end
end

% fill in descendants
for i = 1:nLOC
    p = 1;
    m = 2;
    for j = 1:nIND       
       ind = ibd(i,j,p);
       if( ind > 0 )
           allele = output_fams(i,ind,p);
       else
           allele = output_fams(i,-ind,m);
       end
       output_fams(i,j,p) = allele;

       ind = ibd(i,j,m);
       if( ind > 0 )
           allele = output_fams(i,ind,p);
       else
           allele = output_fams(i,-ind,m);
       end
       output_fams(i,j,m) = allele;      
    end
end












