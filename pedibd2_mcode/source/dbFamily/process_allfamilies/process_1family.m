function [output, state_code, error] = process_1family(family_range, family_genotype, parameters, option)

        error = 0;
        state_code(1:10) = 0;
        output = [];
        global debug_mode;
   
        % remove typing errors
        if( option(1) && error == 0 )
            [family_genotype mendel_error error] ...
                = remove_typing_error(family_genotype, family_range.structure);
            if( error ~= 0 )
                disp('error in family structures');
                return;
            end
            if( error == 0 )
                state_code(10) = 1;
            end
        end
        
        % descent graph & kinship
        if( option(1) && error == 0 )
            [all_inheritance, ...
                kinship, kinship2, kinship2ex, ...
                singleLINEAL, ...
                error] ...
                = generate_inheritance_1family(family_range);
            if( error == 0 )
                state_code(1) = 1;
            end
        end
        
        
        % posterior & viterbi
        if( option(1) && error == 0 && state_code(1) == 1 )           
            [oblist posterior posteriorIBD1 viterbi error] ...
                = posterior_decode_1family ...
                (family_range, all_inheritance, parameters, family_genotype);
            if( error == 0 )
                state_code(2) = 1;
            end
        end
        

        
        % segmentation
        if( option(2) && error == 0 && state_code(2) == 1 )
            [pairs, triples_view, error] ...
                = cut_chromosome ...
                (family_range, kinship, kinship2ex, ...
                family_genotype, viterbi, posterior, posteriorIBD1, ...
                oblist, parameters);
            if( error == 0 )
                state_code(3) = 1;
            end
        end
        
        
        % global IBD
        if( option(2) && error == 0 && state_code(3) == 1 )        
            [assignment error] ...
                = generate_consistent_assignment_2methods ...
                (pairs, ...
                triples_view, ...
                family_range, ...
                kinship2ex, singleLINEAL, option );
            if( error == 0 )
                state_code(4) = 1;
            end
        end
        
        % merge assignment 
        if( option(2) && error == 0 && state_code(4) == 1 )
            [merged_assignment error] = merge_segments_secure(assignment, family_range);
            if( error == 0 )
                state_code(5) = 1;
            end
        end
        
        % relabel founder alleles to reduce recombination
        if( option(2) && error == 0 && state_code(5) == 1 )
            [merged_assignment.alleles_all merged_assignment.recombination error]...
                = reassign_segments(merged_assignment, family_range, parameters);
            if( error == 0 )
                state_code(6) = 1;
            end
        end
        
        % linkage score
        % considering moving this function to a separate one
        if( option(3) && error == 0 && state_code(6) == 1 )            
            % calculating linkage score
            [merged_assignment.score error] ...
                = non_parametric_linkage(merged_assignment, family_range);
            if( error == 0 )
                state_code(7) = 1;
            end
        end
        
        if( option(2) && error == 0 && state_code(6) == 1 )
            [merged_assignment.alleles_all error] ...
                = impute_missing(merged_assignment, family_range, kinship2ex);
            if( error == 0 )
                state_code(9) = 1;
            end
        end
        
        % haplotyping
        if( option(2) && error == 0 && state_code(6) == 1 )           
            [inheritance, haplotype, ...
                consistency, merged_assignment.inconsistency, ...
                error] ...
                = assign_haplotype( ...
                merged_assignment, ...
                family_range, ...
                parameters.sampled_markerlist, ...
                ...
                family_genotype, ...
                family_range.structure(:,8), ...
                parameters.sampled_markerlist);
            
            if( error == 0 )
                state_code(8) = 1;
            end
        end        
      
        
        if( error == 0 || debug_mode == 1 )
            if( state_code(10) == 1 )
                output.mendel_error = mendel_error;
                output.family_genotype = family_genotype;
            else
                output.mendel_error = [];
                output.family_genotype = [];
            end
            if( state_code(1) == 1 )
                output.all_inheritance = all_inheritance;
                output.kinship = kinship;
                output.kinship2 = kinship2;
                output.kinship2ex = kinship2ex;
                output.singleLINEAL = singleLINEAL;
            else
                output.all_inheritance = [];
                output.kinship = [];
                output.kinship2 = [];
                output.kinship2ex = [];
                output.singleLINEAL = [];
            end
            if( state_code(2) == 1 )
                output.oblist = oblist;
                output.posterior = posterior;
                output.posteriorIBD1 = posteriorIBD1;
                output.viterbi = viterbi;
            else
                output.oblist = [];
                output.posterior = [];
                output.posteriorIBD1 = [];
                output.viterbi = [];
            end
            if( state_code(3) == 1 )
                output.pairs = pairs;
                output.triples = triples_view;
            else
                output.pairs = [];
                output.triples = [];
            end
            if( state_code(4) == 1 )
                output.assignment = assignment;
            else
                output.assignment = [];
            end
            if( state_code(5) == 1 )
                output.merged_assignment = merged_assignment;
            else
                output.merged_assignment = [];
            end
            if( state_code(8) == 1 )
                output.inheritance = inheritance;
                output.haplotype = haplotype;
                output.consistency = consistency;
            else
                output.inheritance = [];
                output.haplotype = [];
                output.consistency = [];
            end
        end
        

end
















