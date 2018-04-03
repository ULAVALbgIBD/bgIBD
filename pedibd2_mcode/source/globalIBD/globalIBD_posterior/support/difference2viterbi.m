    
function diff = difference2viterbi(assignment, vit)
    
    diff = [];
    [nseg1 ngeno1 c] = size(assignment);
    [nseg2 ngeno2 ngeno3] = size(vit);
    if( isempty(vit) )
        return;
    end
    if( nseg1 ~= nseg2 || ngeno1 ~= ngeno2 || ngeno2 ~= ngeno3 || c ~= 2 )
        return;
    end
    nseg = nseg1;
    ngeno = ngeno1;
    diff = zeros(nseg1,1);
    for i = 1:nseg
        pair = allele2pair(reshape(assignment(i,1:ngeno,1:2), [ngeno,2]));
        diff(i) = nnz(pair-reshape(vit(i,1:ngeno,1:ngeno), [ngeno,ngeno]));
    end
    
    
end