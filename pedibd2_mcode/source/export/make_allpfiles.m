
function make_allpfiles(source, target)

temp = pwd;
cd(target);

allmfiles = find_allmfiles(source);
for i = 1:length(allmfiles)
    if( isempty(strfind(allmfiles{i}, 'make_allpfiles.m')) )
        disp(allmfiles{i});
        pcode(allmfiles{i});
    end
end

allcfiles = find_allcfiles(source);
for i = 1:length(allcfiles)
    disp(allcfiles{i});
    copyfile(allcfiles{i});
end


cd(temp);

end