function remap = remapFOUNDER(assignment)



[nIND, nFIELDS] = size(assignment);
    
% re-map alleles to founder allele format
% FIRST ENCOUNTERED individuals will be used as founder codes
% not necessarily founder, just for set-leader purposes
tagged(1:nIND,1:2) = false;
remap = zeros(nIND,2);
for i = 1:nIND
    p = assignment(i,1);
    m = assignment(i,2);
    if( p ~= 0 )
        reP = (assignment(1:nIND,1:2) == p);
        if( all(all(~tagged(reP))) )
            remap(reP) = -i;
            tagged(reP) = true;
        end
    end
    if( m ~= 0 )
        reM = (assignment(1:nIND,1:2) == m);
        if( all(all(~tagged(reM))) )
            remap(reM) = i;
            tagged(reM) = true;
        end
    end
end



end


