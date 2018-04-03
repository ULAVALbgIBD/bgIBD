function [input output error] = ...
    process_allfamilies(input, parameters, option)

    global debug_mode;
    error = 0;

    output = [];
    
    
    for family_id = 1:length(input.family_range)       
        
        
        [output1family, ~, error] = ...
            process_1family_secure ...
            (input.family_range{family_id}, ...
            input.family_genotype{family_id}, ...
            parameters, option);               
     
        if( error ~= 0 && debug_mode == 0 )
            disp('error in processing families');
            return;
        end          
        
        input.all_inheritance{family_id} = output1family.all_inheritance;
        input.kinship{family_id} = output1family.kinship;
        input.kinship2{family_id} = output1family.kinship2;
        input.kinship2ex{family_id} = output1family.kinship2ex;
        input.singleLINEAL{family_id} = output1family.singleLINEAL;
        input.oblist{family_id} = output1family.oblist;
        input.family_genotype{family_id} = output1family.family_genotype;

        output.paths.posterior{family_id} = output1family.posterior;
        output.paths.posteriorIBD1{family_id} = output1family.posteriorIBD1;
        output.paths.viterbi{family_id} = output1family.viterbi;

        output.assignment{family_id} = output1family.assignment;
        output.merged_assignment{family_id} = output1family.merged_assignment;
        output.inheritance{family_id} = output1family.inheritance;
        output.haplotype{family_id} = output1family.haplotype;
        output.consistency{family_id} = output1family.consistency;  
        output.mendel_error{family_id} = output1family.mendel_error;
        
        disp(' ');
        if( debug_mode == 1 )
            % in case of debug_mode, run only one family, export pairs
            output.paths.pairs = output1family.pairs;
            output.paths.triples = output1family.triples;
            return;
        end
        
    end
    

    
end



















