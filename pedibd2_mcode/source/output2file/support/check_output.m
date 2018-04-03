function error = check_output(option)

error = 0;

fid = fopen([pwd, '/dbFamily.log'], 'w');
if( fid < 0 )
    error = 1;
    disp(['access denied to folder: ', pwd]);
    disp('cannot create output file');
    return;
else
    print_header(fid, option);
    fclose(fid);
end

delete([pwd, '/dbFamily.log']);

end

