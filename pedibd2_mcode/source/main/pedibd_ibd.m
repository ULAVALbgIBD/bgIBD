function pedibd_ibd(pedigree_file, genotype_file, marker_file)

    nIn = nargin;
        
    if( nIn ~= 3 )
        disp('require 3 input files, pedigree file, genotype file, marker file');
        return;
    end

    option = false(4,1);
    option(1:3) = false;
    option(1) = true;
    
    global release;
    release = 1;
    
    INfileWRAP(pedigree_file, genotype_file, marker_file, [], option);
    
    clear global release;
    
end