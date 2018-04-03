function [ output_fams ] = generate_random_ibd( family, nMARKER, recfrac )

nIND = length(family(:,1));

ibd = zeros(nMARKER, nIND, 2);

% assign first locus for all individuals
% paternal positive
for j = 1:nIND
    if( family(j,3) == 0 )
        ibd(1,j,1) = j;
    else
        if( rand() > 0.5 )
            ibd(1,j,1) = family(j,3);
        else
            ibd(1,j,1) = (-1)*family(j,3);
        end
    end
    if( family(j,4) == 0 )
        ibd(1,j,2) = -j;
    else
        if( rand() > 0.5 )
            ibd(1,j,2) = family(j,4);
        else
            ibd(1,j,2) = (-1)*family(j,4);
        end
    end
end

for i = 2:nMARKER
   p = 1;
   m = 2;
   for j = 1:nIND
       if( family(j,3) == 0 )
           ibd(i,j,p) = ibd(i-1,j,p);
       else
           if( rand() > recfrac(i-1) )
               % different from previous
               ibd(i,j,p) = -ibd(i-1,j,p);
           else
               % same as previous
               ibd(i,j,p) = ibd(i-1,j,p);
           end
       end
       if( family(j,4) == 0 )
           ibd(i,j,m) = ibd(i-1,j,m);
       else
           if( rand() > recfrac(i-1) )
               ibd(i,j,m) = -ibd(i-1,j,m);
           else
               ibd(i,j,m) = ibd(i-1,j,m);
           end
       end
   end
end

output_fams = ibd;










