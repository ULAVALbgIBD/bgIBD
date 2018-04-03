function correct = compute_likelihood(current, i, j, match, match_prob, prob_IBD1, valid_triples, triple_viewpair, max_likelihood, map)

    correct = 1;

    global assignment;
    global all_assignment;
    global mismatch;
    global likelihood;
    global paternal;
    global maternal;    
    
    all_assignment(current,1) = paternal(current,i);
    all_assignment(current,2) = maternal(current,j);    

    order = map.reverse_list(current);

    if( order > 1 )
        mismatch(order) = mismatch(order - 1);
        likelihood(order) = likelihood(order - 1);
    end
    if( order == 1 )
        mismatch(order) = 0;
        likelihood(order) = 1;
    end       

    if( order <= 0 )
        correct = 1;
        return;
    end
    
    if( order > 0 )
        % only count genotyped individuals
        assignment(order,1:2) = all_assignment(current,1:2);
        for k = 1:order-1
            ibd_num = count_ibd(assignment(order,1:2), assignment(k,1:2));
            if( ibd_num ~= match(order, k) )
                mismatch(order) = mismatch(order) + 1;
            end
            if( ibd_num == 1 )
                [a, b] = count_ibd1(assignment(order,1:2), assignment(k,1:2));
                likelihood(order) = likelihood(order) * prob_IBD1(order,k,a,b);
            else
                likelihood(order) = likelihood(order) * match_prob(order,k,ibd_num+1);
            end
            if( likelihood(order) < max_likelihood || likelihood(order) == 0 )
                correct = 0;   
                break;
            end
        end
        for ii = 1:order-1
            ibd_num = count_ibd(assignment(order,1:2), assignment(ii,1:2));
            if( ibd_num ~= 1 )
                continue;
            end
            for jj = ii+1:order-1
                if( valid_triples(order,ii,jj) <= 0 )
                    continue;
                end
                view = triple_viewpair(order,ii,jj);
                ibd_num = count_ibd(assignment(order,1:2), assignment(jj,1:2));
                if( ibd_num ~= 1 )
                    continue;
                end
                ibd_num = count_ibd(assignment(ii,1:2), assignment(jj,1:2));
                if( ibd_num ~= 1 )
                    continue;
                end
                [a1,b1] = count_ibd1(assignment(order,1:2), assignment(ii,1:2));
                [a2,b2] = count_ibd1(assignment(order,1:2), assignment(jj,1:2));
                [a3,b3] = count_ibd1(assignment(ii,1:2), assignment(jj,1:2));
                if( view > 0 )
                    if( a1 ~= a2 )
                        likelihood(order) = likelihood(order) * 0.99;
                    end
                end
                if( view < 0 )
                    if( a1 == a2 )
                        likelihood(order) = likelihood(order) * 0.99;
                    end
                end
                if( likelihood(order) < max_likelihood || likelihood(order) == 0 )
                    correct = 0;
                    break;
                end
            end
            if( likelihood(order) < max_likelihood || likelihood(order) == 0 )
                correct = 0;
                break;
            end                
        end
    end

end












