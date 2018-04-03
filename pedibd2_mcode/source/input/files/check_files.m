function error = check_files(pedigree_file, genotype_file, marker_file, hotspots_file)

error = 0;

code = exist(pedigree_file, 'file');
if( code ~= 2 )
    error = 1;
    disp('cannot load pedigree file, no such file exists');
    disp(pedigree_file);
    return;
end

code = exist(genotype_file, 'file');
if( code ~= 2 )
    error = 1;
    disp('cannot load genotype file, no such file exists');
    disp(genotype_file);
    return;
end

code = exist(marker_file, 'file');
if( code ~= 2 )
    error = 1;
    disp('cannot load marker file, no such file exists');
    disp(marker_file);
    return;
end

if( ~isempty(hotspots_file) )
    code = exist(hotspots_file, 'file');
    if( code ~= 2 )
        error = 1;
        disp('cannot load hotspots file, no such file exists');
        disp(hotspots_file);
        return;
    end
end

end