function [error] = check_integer(data, range)

error = 0;
if( isempty(data) )
    error = 1;
    return;
end

if( ~isa(data, 'numeric') )
    error = 1;
    return;
end

if( any(any(isnan(data))) )
    error = 1;
    disp('data contain empty cells');
    return;
end

if( any(any(fix(data) ~= data)) )
    error = 1;
    return;
end

if( any(any((data<range(1) | data>range(2)))) )
    error = 1;
    return;
end
        

end