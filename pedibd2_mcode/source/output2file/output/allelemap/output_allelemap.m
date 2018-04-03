function error = output_allelemap(families, ...
    assignment, parameters, option)



suc = 0;
error = 0;

directory = [pwd, '\alleleMap'];
if( exist(directory, 'dir') == 0 )
    suc = mkdir(directory);
else
    suc = 1;
end

if( suc ~= 1 )
    error = 1;
    disp('cannot create output folder');
    return;
end

disp(['alleleMap depository: ', directory]);
disp(' ');

fid = fopen([directory, '/readme.txt'], 'w');
if( fid < 0 )
    error = 1;
    return;
else
    print_header(fid, option);
    fprintf(fid, '\n*****      ');
    fprintf(fid, 'alleleMap deposited for debugging purposes ...');
    fprintf(fid, '\n');
    fprintf(fid, '\n');
    fprintf(fid, '\n*****      ');
    fprintf(fid, 'alleleMaps are internal data structrues of Ped-IBD');
    fprintf(fid, '\n*****      ');
    fprintf(fid, 'coding inheritance patterns of a family for each non-recombinant chromosomal segment');
    fprintf(fid, '\n*****      ');
    fprintf(fid, '*p and *m are paternal and maternal alleles of an individual');
    fprintf(fid, '\n*****      ');
    fprintf(fid, 'alleleMap specifies the distribution of distinct alleles in a family');
    fprintf(fid, '\n');
    fclose(fid);
end


nfam = length(families);
nconf = length(assignment);
if( nfam <= 0 || nfam ~= nconf )
    error = 1;
    disp('error in global IBD configuration');
    return;
end

for i = 1:nfam
    if( isempty(families{i}) || isempty(assignment{i}) )
        error = 1;
        disp('error in generating alleleMap');
        return;
    end
    family_id = families{i}.family_id;
    filehead = [directory, '/map', num2str(i), '_fam', num2str(family_id)];
    output_allelemap1family(filehead, ...
        families{i}, assignment{i}, parameters);
end


end














