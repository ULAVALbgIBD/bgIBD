

% skip data_script
function dir_all(directory)

global allfiles;
global fcount;

list = dir(directory);
len = length(list);
for i = 1:len
    if( ~strcmp(list(i).name, '.') && ~strcmp(list(i).name, '..') )
        fullname = [directory, '\', list(i).name];
        if( strcmp(list(i).name, 'data_script') )
            % exclude certain folders
            continue;
        end
        if( strcmp(list(i).name, 'export') )
            % exclude certain folders
            continue;
        end        
        if( list(i).isdir == 1 )
            dir_all(fullname);
        else
            fcount = fcount + 1;
            allfiles{fcount} = fullname;
        end
    end
end


end