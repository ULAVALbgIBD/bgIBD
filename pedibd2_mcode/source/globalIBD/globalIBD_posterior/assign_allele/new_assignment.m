function new_assignment(current, paternal_alleles, maternal_alleles)

    global paternal;
    global maternal;
    global paternal_num;
    global maternal_num;
    global all_assignment;

    source_alleles(1,1:2) = paternal_alleles(current,1:2);
    source_alleles(2,1:2) = maternal_alleles(current,1:2);

    all_assignment(current,1:2) = [-current,current];
    i = 1;
    for j = 1:2
        if( source_alleles(i,j) < 0 )
            paternal(current,j) = all_assignment(-source_alleles(i,j), 1);
        end
        if( source_alleles(i,j) > 0 )
            paternal(current,j) = all_assignment(source_alleles(i,j), 2);
        end
    end
    if( paternal(current,1) == paternal(current,2) )
        paternal_num(current) = 1;
    else
        paternal_num(current) = 2;
    end
    i = 2;
    for j = 1:2
        if( source_alleles(i,j) < 0 )
            maternal(current,j) = all_assignment(-source_alleles(i,j), 1);
        end
        if( source_alleles(i,j) > 0 )
            maternal(current,j) = all_assignment(source_alleles(i,j), 2);
        end
    end
    if( maternal(current,1) == maternal(current,2) )
        maternal_num(current) = 1;
    else
        maternal_num(current) = 2;
    end


end
