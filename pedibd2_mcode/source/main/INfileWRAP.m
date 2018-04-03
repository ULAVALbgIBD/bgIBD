


function [error] = INfileWRAP(pedigree_file, genotype_file, marker_file, hotspots_file, option)

    error = 0;

    % make sure all files are existing for fast processing
    if( exist('cForward', 'file') ~= 3 )
        disp('cForward.c not compiled');
        disp('Compile all c source codes first using mex:');
        disp('cForward.c, cBackward.c, cViterbi.c');
        error = 1;
        return;
    end
    if( exist('cBackward', 'file') ~= 3 )
        disp('cBackward.c not compiled');
        disp('Compile all c source codes first using mex:');
        disp('cForward.c, cBackward.c, cViterbi.c');        
        error = 1;
        return;
    end
    if( exist('cViterbi', 'file') ~= 3)
        disp('cViterbi.c not compiled');
        disp('Compile all c source codes first using mex:');
        disp('cForward.c, cBackward.c, cViterbi.c');        
        error = 1;
        return;
    end 
    
    nIn = nargin;   
    if( nIn ~= 5 )
        error = 1;
        disp('internal error');
        return;        
    end

    if( isempty(option) || ~islogical(option) )
        error = 1;
        disp('internal error');
        return;
    end
    
    disp(' ');
    disp('loading files ...');
    disp(' ');
    
    [error] = check_files(pedigree_file, genotype_file, marker_file, hotspots_file);    
    if( error ~= 0 )
        return;
    end
    
    % Nov 23, 2014
    % modify this file to adapt to string styple IDs
    % create a wrapper for all family and individual IDs
    [pedigree_info, IDmapper, error] = importPED(pedigree_file);
    if( error )
        disp('error in pedigree file');
        return;
    end
    
    try
        marker_list = importdata(marker_file);
    catch ME
        disp(ME.message);
        disp('error in marker file: all fields must be integers');
        return;
    end
    
    [genotype_data, error] = importGEN(genotype_file, IDmapper);
    if( error )
        disp('error in genotype file');
        return;
    end
    
    if( ~isempty(hotspots_file) )
        try
            hotspots = importdata(hotspots_file);
        catch ME
            disp(ME.message);
            disp('error in hotspots file: all fields must be integers');
            return;
        end
    else
        hotspots = [];
    end
    
    
    if( option(3) == 1 )
        disp('launch PedIBD ...');
    else
        disp('launch PedIBD ...');
    end
    disp(' ');
    
    time = cputime;
    
    [error] = check_output(option);
    if( error ~= 0 )
        disp('fail to access the output folder');
        return;
    end

    [input, output, ~, parameters, error] = ...
        dbFamily_secure(pedigree_info, ...
        genotype_data, ...
        marker_list, ...
        hotspots, option);
    
    if( error ~= 0 )        
        disp('error in processing data: skipped');
        return;
    end
    
    [header] = generate_header(pedigree_file, genotype_file, marker_file);
    output_result(header, input, output, parameters, IDmapper, option);    
  
    fclose('all');
    
    display(['total time: ', num2str(cputime - time), ' seconds']);
    
end

function [pedigree_data, IDmapper, error] = importPED(pedigree_file)

    error = false;
    fileID = fopen( ...
        pedigree_file, 'r');
    firstLINE = fgetl(fileID); 
    token = strsplit(firstLINE, {'\t', ' '});
    numFIELD = length(token);
    fclose('all');
    if( numFIELD ~= 7 )
        error = true;
        disp('Wrong pedigree file format, require 7 fields');
        return;
    end
    fileID = fopen( ...
        pedigree_file, 'r');    
    format = [repmat('%s\t', [1 4]), repmat('%f\t', [1 3])];
    data = textscan(fileID, format, ...
        'HeaderLines', 0, 'Delimiter', '\t ');
    fclose('all');
    FAMILYid = unique(data{1});
    INDid = unique(vertcat(data{2:4}));
    pedigree_data(:,5:7) = horzcat(data{5:7});  
    INDid = INDid(~ismember(INDid,'0'));
    INDorder = (1:length(INDid))';
    INDorder(end+1) = 0;
    INDid{end+1} = '0';
    [bit, list] = ismember(data{1}, FAMILYid);
    if( ~all(bit) )
        error = true;
        disp('Error in family ID');
        return;
    end
    pedigree_data(bit,1) = list;
    for i = 1:3
        [bit, list] = ismember(data{i+1}, INDid);
        if( ~all(bit) )
            error = true;
            disp('Error in individual ID');
            return;
        end
        pedigree_data(bit,i+1) = INDorder(list(bit));
    end
    IDmapper.FAMILYid = FAMILYid;
    IDmapper.INDid = INDid;
    IDmapper.INDorder = INDorder;

end


function [genotype_data, error] = importGEN(genotype_file, IDmapper)

    error = false;
    fileID = fopen( ...
        genotype_file, 'r');
    firstLINE = fgetl(fileID); 
    token = strsplit(firstLINE, {'\t', ' '});
    token = token(~ismember(token, ''));
    numFIELD = length(token);
    fclose('all');
    if( numFIELD < 6 )
        error = true;
        disp('Wrong genotype file format, require >= 6 fields');
        return;
    end
    fileID = fopen( ...
        genotype_file, 'r');    
    format = [repmat('%s\t', [1 4]), repmat('%f\t', [1 numFIELD - 4])];
    data = textscan(fileID, format, ...
        'HeaderLines', 0, 'Delimiter', '\t ');
    fclose('all');
    genotype_data(:,5:numFIELD) = horzcat(data{5:numFIELD});  

    FAMILYid = IDmapper.FAMILYid;
    INDid = IDmapper.INDid;
    INDorder = IDmapper.INDorder;
    [bit, list] = ismember(data{1}, FAMILYid);
    if( ~all(bit) )
        error = true;
        disp('Error in family ID, doesn''t match pedigree file');
        return;
    end
    genotype_data(bit,1) = list;
    for i = 1:3
        [bit, list] = ismember(data{i+1}, INDid);
        if( ~all(bit) )
            error = true;
            disp('Error in individual ID, doesn''t match pedigree file');
            return;
        end
        genotype_data(bit,i+1) = INDorder(list(bit));
    end

end



