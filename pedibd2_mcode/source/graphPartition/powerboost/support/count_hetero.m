function [het error] = count_hetero(seq)

het = [];
error = 0;

if( ndims(seq) ~= 2 )
    error = 1;
    return;
end

[r c] = size(seq);
if( r <= 0 || c ~= 2 )
    error = 1;
    disp('dimension error');
    return;
end
len = r;

het = false(len,1);

s1 = seq(1:len,1);
s2 = seq(1:len,2);

missing = ( s1 == 0 | s2 == 0 );
het(1:len) = ( s1 ~= s2 );
het(missing) = false;


end
