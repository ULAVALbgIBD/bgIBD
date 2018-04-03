

function [viterbi, posterior, posteriorIBD1, oblist, emission, inconsistency, error] = posterior_decode_1pair(family, all_inheritance, parameters, family_genotype, input_pair)

        error = 0;
        emission = [];
        inconsistency = [];
        global debug_mode;
    
        display(['          family ',num2str(family(input_pair(1),12)), ': processing relative pair: ',num2str(family(input_pair(1),8)), ', ', num2str(family(input_pair(2),8))]);
        time = cputime;
        paths.allpaths = linkage_compute_allele_verify(all_inheritance, family, input_pair(1), input_pair(2));
        sampled_markerlist = parameters.sampled_markerlist;
        [num_markers c] = size(sampled_markerlist);
        if( num_markers <= 0 || c ~= 2 )
            error = 1;
            disp('error in marker map');
            return;
        end
        [paths.pr paths.tpl paths.tpl_log status status_bypass error] = transition(paths.allpaths, parameters); 
        if( error ~= 0 )
            disp('error in generating transition probability');
            return;
        end
        if( debug_mode == 1 )
            disp(['                    transition time: ',num2str(cputime - time)]);
        end
        time = cputime;        
        [oblist emission_option error] = generate_oblist(family_genotype, input_pair(1), input_pair(2), sampled_markerlist, 2);
        if( error ~= 0 )
            disp('error in generating IBS status');
            return;
        end
        if( debug_mode == 1 )
            disp(['                    observation time: ',num2str(cputime - time)]);
        end
        time = cputime;
        if( status_bypass == 0 )
            [paths.epl paths.epl_log, error] = emission_pair(parameters, family_genotype, input_pair(1), input_pair(2), family, emission_option);
            if( error ~= 0 )
                return;
            end
            emission = paths.epl;
        end
        if( debug_mode == 1 )
            disp(['                    emission time: ',num2str(cputime - time)]);
        end
        
        time = cputime;
        [viterbi, posterior, posteriorIBD1, ~, error] = posterior_path_bypass(paths, oblist, status);   
        if( error ~= 0 )
            disp('error in decoding IBD');
            return;
        end
        if( debug_mode == 1 )
            disp(['                    decoding time: ',num2str(cputime - time)]);
        end
        
end


