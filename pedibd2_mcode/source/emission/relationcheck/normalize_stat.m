
function value = normalize_stat(dist, ob)

    value = 1;
    
    error_code = 'relationship check'
    
    len = length(dist);
    if( len <= 0 )
        disp(['emission probability not match data', error_code]);
        return;
    end
    if( ob > len )
        disp(['emission probability not match data', error_code]);
        return;
    end
    if( abs(sum(dist) - 1) > 0.00001 )
        disp(['emission probability not match data', error_code]);
        return;
    end
    
    mu = 0;
    std = 0;
    va = 0;
 
    for i = 1:len
        mu = mu + i * dist(i);
    end

    for i = 1:len
        va = va + ((i-mu)^2) * dist(i);
    end
    if( va == 0 && ob ~= mu )
        disp('error in emission probability');
        return;
    end
    std = sqrt(va);
    
    if( ob == mu )
        value = 0;
    else
        value = (ob - mu)/std;
    end

end