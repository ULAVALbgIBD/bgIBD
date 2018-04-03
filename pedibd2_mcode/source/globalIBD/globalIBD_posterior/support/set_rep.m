function [array, root_rep] = set_rep(ina, ix)

array = ina;

if( ix > 0 )
    rep = array(ix);
else
    rep = -array(-ix);
end

if( abs(array(abs(rep))) ~= abs(rep) )
    [array, root_rep] = set_rep(array, rep);
else
    root_rep = rep;
end


if( array(ix) > 0 )
    array(ix) = root_rep;
else
    array(-ix) = root_rep;
end

edn