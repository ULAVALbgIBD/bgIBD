function allcfiles = find_allcfiles(directory)
global allfiles;
global fcount;
fcount = 0;
allfiles = [];

dir_all(directory);

mcount = 0;
for i = 1:fcount
    if( ~isempty(findstr(allfiles{i}, '.c')) )
        mcount = mcount + 1;
        allcfiles{mcount} = allfiles{i};
    end
end

clear global allfiles;
clear global fcount;

end
