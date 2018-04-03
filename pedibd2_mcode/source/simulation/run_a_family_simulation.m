

function [simulation] = run_a_family_simulation(family_range, pedigree, random_pedigree, parameters)

    input_pairlist = family_range.pairs;
    sampled_markerlist = parameters.sampled_markerlist;
    pedigree_all_missing = pedigree.structure;

    simulation.paths.align_pair = [];
    
    [input.paths.epl input.paths.epl_sc] = emission(parameters);
    
    for s_id = 1:length(input_pairlist(:,1))

        s_id
        
        f = input_pairlist(s_id,3);


        input_allpaths = linkage_compute_allele(pedigree_all_missing, input_pairlist(s_id,3), input_pairlist(s_id,1), input_pairlist(s_id,2));
        [input.paths.pr input.paths.tpl status] = transition(input_allpaths, sampled_markerlist);
        

        for i = 1:1

            [random_ibd_result, random_ibs_result] = simulation_oblist_pair(s_id, input_pairlist, pedigree_all_missing, random_pedigree{f,2}, random_pedigree{f,3}); 
            ref = find_breakpoints(random_ibd_result', sampled_markerlist);
            ref = ref';
            ref = [ref-1,ref];
            simulation.ref_breakpoint{i}{s_id} = ref;
            simulation.ibs{i}{s_id} = random_ibs_result;
            % extract input from simulated data


            simulation.paths.posterior{i}{s_id} = posterior_path_bypass(input.paths, simulation.ibs{i}{s_id}, status);
            [temp1 temp2] = max(simulation.paths.posterior{i}{s_id}');
            simulation.paths.output_state{i}{s_id} = [temp2',temp1'];
            
            simulation.paths.output_breakpoint{i}{s_id} = find_breakpoints(simulation.paths.output_state{i}{s_id}, sampled_markerlist);

            [align, simulation.paths.false_pair(i,1:4)] = recombination_align(simulation.paths.output_breakpoint{i}{s_id}, simulation.ref_breakpoint{i}{s_id}, sampled_markerlist);   
            if( ~isempty(align(:,1)) )
                align(:,5) = i;
                align(:,6) = s_id;
                simulation.paths.align_pair = [simulation.paths.align_pair; align];
            end
        end
        
    end
    
    for i = 1:1
    
        simulation.paths.assignment = generate_consistent_assignment_likelihood( simulation.paths.posterior{i}, family_range, pedigree, 0 );
    
    end

end






