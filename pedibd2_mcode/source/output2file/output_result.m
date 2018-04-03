function [error] = output_result(header, input, output, parameters, IDmapper, option)

    % make sure that all results are generated for each family
    error = 0;
    
    disp(' ');
    disp('exporting results ...');
    disp(' ');
    disp(['output folder: ', pwd]);
    disp(' ');
    
    
    if( option(1) )
        fid = fopen([header, '_IDmapper.txt'], 'w');
        print_header(fid, option);
        fprintf(fid, '# Family and Individual IDs are converted for easy display\n\n');
        FAMILYid = IDmapper.FAMILYid;
        INDid = IDmapper.INDid;
        INDorder = IDmapper.INDorder;
        fprintf(fid, '%s\t%s\n', 'orginalFAMILYid', 'newFAMILYid');
        for i = 1:length(FAMILYid)
            fprintf(fid, '%s\t%d\n', FAMILYid{i}, i);
        end
        fprintf(fid, '\n%s\t%s\n', 'orginalINDIVIDUALid', 'newINDIVIDUALid');
        for i = 1:(length(INDid)-1)
            fprintf(fid, '%s\t%d\n', INDid{i}, INDorder(i));
        end
        fclose(fid);
    end
    
    if( option(1) && ~isempty(input.kinship) )
        fid = fopen([header, '_kin.txt'], 'w');
        if( fid < 0 )
            disp('cannot create output file, access denied');
            disp([header, '_kin.txt']);
        else
            print_header(fid, option);
            error = output_kinship(fid, input.family_range, input.kinship);
            fclose(fid);
        end        
    end
    
    if( option(1) && ~isempty(output.paths.posterior) )
        fid = fopen([header, '_ibd.txt'], 'w');
        if( fid < 0 )
            disp('cannot create output file, access denied');
            disp([header, '_ibd.txt']);            
        else
            print_header(fid, option);
            error = output_posterior(fid, input.family_range, output.paths.posterior, parameters);
            fclose(fid);
        end        
    end
    
    if( option(2) && ~isempty(output.haplotype) )
        fid = fopen([header, '_hap.txt'], 'w');
        if( fid < 0 )
            disp('cannot create output file, access denied');
            disp([header, '_hap.txt']);            
        else
            print_header(fid, option);
            error = output_haplotype(fid, input.family_range, output.haplotype, parameters.sampled_markerlist);
            fclose(fid);
        end       
    end
    
    if( option(3) && ~isempty(output.combined_score) )
        fid = fopen([header, '_lnk.txt'], 'w');
        if( fid < 0 )
            disp('cannot create output file, access denied');
            disp([header, '_lnk.txt']);            
        else
            print_header(fid, option);
            error = output_linkage_score(fid, input.family_range, output.combined_score, parameters);
            fclose(fid);
        end        
    end
   
    % added function to ped-ibd
    if( option(2) && ~isempty(output.merged_assignment) )
        fid = fopen([header, '_rmb.txt'], 'w');
        if( fid < 0 )
            disp('cannot create output file, access denied');
            disp([header, '_rmb.txt']);            
        else
            print_header(fid, option);
            error = output_recombination(fid, input.family_range, output.merged_assignment, parameters);
            fclose(fid);
        end        
    end
    
    if( option(2) && ~isempty(output.merged_assignment) )
        fid = fopen([header, '_inh.txt'], 'w');
        if( fid < 0 )
            disp('cannot create output file, access denied');
            disp([header, '_inh.txt']);            
        else
            print_header(fid, option);
            error = output_inheritanceH(fid, input.family_range, output.merged_assignment, parameters);
            fclose(fid);
        end           
    end
    
    if( option(4) && ~isempty(output.merged_assignment) )
        error = output_allelemap(input.family_range, ...
            output.merged_assignment, parameters, option);
    end    
    
    fclose('all');

end






















