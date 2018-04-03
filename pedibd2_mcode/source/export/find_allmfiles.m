function allmfiles = find_allmfiles(directory)
global allfiles;
global fcount;
fcount = 0;
allfiles = [];

dir_all(directory);

mcount = 0;
for i = 1:fcount
    if( ~isempty(findstr(allfiles{i}, '.m')) && isempty(findstr(allfiles{i}, '.mex')) )        
        mcount = mcount + 1;
        allmfiles{mcount} = allfiles{i};
    end
end

clear global allfiles;
clear global fcount;

end
