
function [viterbi, posterior, posteriorIBD1, ibd_option, error] = posterior_path_bypass(paths, oblist, status)

    error = 0;
    global debug_mode;
    
    ibd_option = generate_ibd_option();

    sum_all = sum(sum(status));
    sum1_ = sum(status(1,:));
    sum2_ = sum(status(2,:));
    sum_1 = sum(status(:,1));
    sum_2 = sum(status(:,2));


    if( sum_all == 0 || sum1_ == 1 || sum2_ == 1 || sum_1 == 1 || sum_2 == 1 )              
  
        % bypass the decoding process for special pairs to save time
        % parent-child pairs must be bypassed
        % decoding by inheritance path is not correct
        
        post_temp = [];
        vit_temp = [];
        post_temp(1:length(oblist), 1:3) = 0;
        vit_temp(1:length(oblist),1:2) = 0; %column two is background decoding
        postIBD1_temp(1:length(oblist), 1:2, 1:2) = 0;
        
        % not related
        if( sum_all == 0 )
            post_temp(1:length(oblist), 1) = 1;            
            vit_temp(:,1) = 1;                        
        end        
        %parent-child relation
        if(  sum1_ == 1 || sum2_ == 1 || sum_1 == 1 || sum_2 == 1 )
            
            post_temp(1:length(oblist), 2) = 1;    
            vit_temp(:,1) = 2;                       
            if( sum1_ == 1 )
                postIBD1_temp(:,1,1) = 0.5;
                postIBD1_temp(:,1,2) = 0.5;
            end
            if( sum2_ == 1 )
                postIBD1_temp(:,2,1) = 0.5;
                postIBD1_temp(:,2,2) = 0.5;
            end
            if( sum_1 == 1 )
                postIBD1_temp(:,1,1) = 0.5;
                postIBD1_temp(:,2,1) = 0.5;
            end
            if( sum_2 == 1 )
                postIBD1_temp(:,1,2) = 0.5;
                postIBD1_temp(:,2,2) = 0.5;
            end
            % assume this only happens between parent-child pairs
            % pairwise relationships are not able to distringuish two
            % alleles
            % therefore, assign them equal probability
               
        end
        posterior = post_temp;
        posteriorIBD1 = postIBD1_temp;
        viterbi = vit_temp;
    else
        time = cputime;       
        [posterior, posteriorIBD1, error] = posterior_path(paths, oblist);
        if( error ~= 0 )
            disp('error in generating posterior probability');
            return;
        end    
        if( debug_mode == 1 )
            disp(['                         posterior decoding time: ',num2str(cputime - time)]);
        end 
        time = cputime; 
        [viterbi error] = ibd_decode_path(paths, oblist, ibd_option.viterbi);
        if( error ~= 0 )
            disp('error in generating viterbi decoding');
            return;
        end
        if( debug_mode == 1 )
            disp(['                         viterbi decoding time: ',num2str(cputime - time)]);
        end         
    end
    
    
    
end