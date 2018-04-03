function descent_graph_interpret(input_genF)

    num = input_genF(1,1);
    genF = input_genF(:,2:end);

    for i = 1:length(genF(:,1))
        
        sum = genF(i,1);
        for j = num + 2:length(genF(i,:))
            sum = sum * 0.5^(genF(i,j));
        end
        fprintf('%f:\t', sum);
        
        fprintf('%2d * id', genF(i,1));
        fprintf('%d', genF(i,2:num+1));
        fprintf(' (');        
        for j = num + 2:length(genF(i,:))
            code = j - num - 1;
            if( genF(i,j) ~= 0 )
                fprintf(' * t');
                for k = 1:num
                    if( bitget(code,k) > 0 )
                        fprintf('%d', k);
                    end
                end
                fprintf('^%d', genF(i,j));
            end            
        end
        fprintf(' )\n');
    end

    
end