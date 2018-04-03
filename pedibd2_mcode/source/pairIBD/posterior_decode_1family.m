function [oblist_ibs posterior posteriorIBD1 viterbi error] = posterior_decode_1family(family_range, all_inheritance, parameters, family_genotype)

        
        error = 0;
        
        oblist_ibs = [];
        posterior = [];
        posteriorIBD1 = [];
        viterbi = [];
        
        global debug_mode;
        
        if( isempty(family_range) )
            disp('error in family structures');
            error = 1;
            return;
        end
        if( length( family_range.pedigree_range_full ) < 1 )
            disp('empty family')
            error = 1;
            return;
        end
        input_list = family_range.pairs;
            
        [nrows, ncols] = size(input_list);        
        
        if( length(family_range.family_range) < 2 )
            % skip if empty input_list
            if( ~isempty(input_list) || nrows > 0 )
                error = 1;
                disp('error in family structures');
                return;
            end
            display(['family ',num2str(family_range.family_id), ': < 2 genotyped individuals: ibd skipped']);
            return;
        end
        
        if( nrows < 1 || ncols ~= 2 )
            display(['family ',num2str(family_range.family_id), ': error in processing the family']);
            error = 1;
            return;
        end
        
        display(' ');
        display(['family ',num2str(family_range.family_id), ': infer IBD chromosomal regions, total ', num2str(nrows), ' genotyped relative pairs']);
        time = cputime;
        
        npairs = nrows;
        sampled_markerlist = parameters.sampled_markerlist;
        [num_markers c] = size(sampled_markerlist);
        if( num_markers <= 0 || c ~= 2 )
            error = 1;
            disp('error in marker map');
            return;
        end        
        
        % initiate memory space to be integer, for less memory use
        viterbi = int8(zeros(npairs, num_markers, 2));
        posterior = single(zeros(npairs, num_markers, 3));
        posteriorIBD1 = single(zeros(npairs, num_markers, 2, 2));
        oblist_ibs = int8(zeros(npairs, num_markers));

        % memory function is not available on some machines
        if( debug_mode == 1 )
            [user sys] = memory;
            ratio = user.MemUsedMATLAB/sys.PhysicalMemory.Total;
            if( ratio > 1 )
                disp('insufficient memory');
                error = 1;
                memory
                return;
            end
        end
        
        for s_id = 1:npairs
            [vit, post, post1, ob, ~, ~, error] ...
                = posterior_decode_1pair(family_range.structure, ...
                all_inheritance, ...
                parameters, ...
                family_genotype, ...
                input_list(s_id,1:2));    
            if( error ~= 0 )
                disp('error in decoding IBD');
                return
            end
            viterbi(s_id, 1:num_markers, 1:2) = vit;
            posterior(s_id, 1:num_markers, 1:3) = post;
            posteriorIBD1(s_id, 1:num_markers, 1:2, 1:2) = post1;
            oblist_ibs(s_id, 1:num_markers) = ob;
        end
        
        [error] = check_posterior1family(family_range, posterior, sampled_markerlist);
        if( error ~= 0 )
            error = 1;
            return;
        end

        display(['generating pairwise IBD costs: ', num2str(cputime - time), ' seconds']);

end








