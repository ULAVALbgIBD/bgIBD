function [p_value inconsistency error] = verify_relationship(oblist, viterbi, emission_option, ibd_option, family, input_pair, family_genotype, parameters)


            [paths.epl paths.epl_sc error] = emission_pair(parameters, family_genotype, input_pair(1), input_pair(2), family, emission_option);
            if( error ~= 0 )
                return;
            end
            emission = paths.epl;            
            [p_value inconsistency error] = categorical_null_likelihood(emission, oblist, viterbi(:,1), emission_option, ibd_option.viterbi);
            if( error ~= 0 )
                return;
            end
            if( debug_mode == 1 || status_bypass == 2 )
                if( p_value < 1 )
                    if( p_value > 0 )
                        if( p_value < 10^-4 )
                            disp(['relationship mis-specification p_value: ', num2str(p_value)]);
                        end
                    else
                        disp(['relationship mis-specification p_value: ', num2str(10^-20)]);
                    end
                end
            end



end

