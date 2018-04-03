

function [error] = check_genotype(all_data)

error = 0;
[nrows ncols] = size(all_data);

if( nrows <= 0 )
    error = 1;
    disp('empty genotype data');
    return;
end

if( ncols <= 6 )
    error = 1;
    disp('empty genotype data');
    return;
end

if( mod(ncols-6,2) ~= 0 )
    error = 1;
    disp('odd genotype columns');
    return;
end

error = check_integer(all_data, [0,inf]);
if( error ~= 0 )
    error = 1;
    disp('wrong genotype data format, all fields must be non-negative integers');
    return;
end

error = check_integer(all_data(:,7:end), [0,2]);
if( error ~= 0 )
    error = 1;
    disp('wrong genotype data format, alleles must be coded 1 and 2, and 0 for missing');
    return;
end
    
    
end