function [ output_fams ] = condense_random_ibd( ibd, nLOC )

nIND = size(ibd,2);

% trace to founders
% assuming pre-ordered
for i = 1:nLOC
   p = 1;
   m = 2;
   for j = 1:nIND
       
       ind = ibd(i,j,p);
       if( ind > 0 )
           allele = ibd(i,ind,p);
       else
           allele = ibd(i,-ind,m);
       end
       while( abs(allele) ~= abs(ind) )
           ibd(i,j,p) = allele;
           ind = ibd(i,j,p);
           if( ind > 0 )
               allele = ibd(i,ind,p);
           else
               allele = ibd(i,-ind,m);
           end           
       end
       
       ind = ibd(i,j,m);
       if( ind > 0 )
           allele = ibd(i,ind,p);
       else
           allele = ibd(i,-ind,m);
       end
       while( abs(allele) ~= abs(ind) )
           ibd(i,j,m) = allele;
           ind = ibd(i,j,m);
           if( ind > 0 )
               allele = ibd(i,ind,p);
           else
               allele = ibd(i,-ind,m);
           end
       end
   end
end

output_fams = ibd;










