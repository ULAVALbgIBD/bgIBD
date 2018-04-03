function [ output_fams ] = generate_random_missing_error( random_genotype, mr, er )

[nLOC, nIND, ~] = size(random_genotype);


output_fams = random_genotype;


for i = 1:nLOC
   p = 1;
   m = 2;
   for j = 1:nIND
       
       if( rand() < er/nIND )
           if( rand() < 0.5 )
               output_fams(i,j,p) = 1;
           else
               output_fams(i,j,p) = 2;
           end
           if( rand() < 0.5 )
               output_fams(i,j,m) = 1;
           else
               output_fams(i,j,m) = 2;
           end
       end
       
       if( rand() < mr )
           
           output_fams(i,j,p) = 0;
           output_fams(i,j,m) = 0;
           
       end
       
   end
end












