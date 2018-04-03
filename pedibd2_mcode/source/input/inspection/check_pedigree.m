

function [error] = check_pedigree(pedigree_info)

error = 0;
[nrows ncols] = size(pedigree_info);

if( nrows <= 0 )
    disp('empty pedigree data');
    error = 1;
    return;
end

if( ncols < 7 )
    disp('wrong pedigree data format, require at least 7 columns');
    error = 1;
    return;
end

error = check_integer(pedigree_info, [0,inf]);
if( error ~= 0 )
    disp('wrong pedigree data, all fields must be non-negative integers');
    error = 1;
    return;
end

error = check_integer(pedigree_info(:,5), [0,2]);
if( error ~= 0 )
    disp('wrong pedigree data, sex field should be coded 1 and 2, 0 for unknown');
    error = 1;
    return;
end

error = check_integer(pedigree_info(:,6), [0,2]);
if( error ~= 0 )
    disp('wrong pedigree data, disease status should be coded 1 and 2, 0 for unknown');
    error = 1;
    return;
end

error = check_integer(pedigree_info(:,7), [0,1]);
if( error ~= 0 )
    disp('wrong pedigree data, column 7 should be coded 0 and 1 for ungenotyped and genotyped individuals');
    error = 1;
    return;
end

end