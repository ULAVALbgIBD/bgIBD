
function epl = id_expand_emission(is, emi)

    num = length(is(1,:));
    states = length(is(:,1));
    
    epl(states^2, length(emi(1,:))) = 0;
    
    temp(1:states,1:states) = 0;
    
    for i = 1:states
        code1 = is(i,:);
        for j = 1:states
            code2 = is(j,:);
            code3 = code1;
            for k = 1:num
                op = k;
                opr = code2(k);
                while( code3(op) ~= code3(code3(op)) )
                    code3(op) = code3(code3(op));
                end
                while( code3(opr) ~= code3(code3(opr)) )
                    code3(opr) = code3(code3(opr));
                end
                
                rep1 = code3(op);
                rep2 = code3(opr);
                
                if( rep1 < rep2 )
                    code3(rep2) = rep1;
                end
                if( rep2 < rep1 )
                    code3(rep1) = rep2;
                end
                
                while( code3(op) ~= code3(code3(op)) )
                    code3(op) = code3(code3(op));
                end
                while( code3(opr) ~= code3(code3(opr)) )
                    code3(opr) = code3(code3(opr));
                end                
            end
            for k = 1:num
                while( code3(k) ~= code3(code3(k)) )
                    code3(k) = code3(code3(k));
                end
            end
            for a = 1:states
                success = 1;
                for b = 1:num
                    if( code3(b) ~= is(a,b) )
                        success = 0;
                        break;
                    end
                end
                if( success == 1 )                 
                    temp(i,j) = a;
                    break;
                end
            end
            for t = 1:length(emi(1,:))
                epl( states * (i-1) + j-1 + 1, t) = emi(a,t);
            end
        end
    end    
end


