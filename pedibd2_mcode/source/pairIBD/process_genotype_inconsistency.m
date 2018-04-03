function average = process_genotype_inconsistency(list)

global debug_mode;

average = [];

count = 0;
summation = 0;
for i = 1:length(list)
    if( ~isempty(list{i}) )
        count = count + 1;
        summation = summation + list{i};
    end
end
if( count > 0 )
    average = summation/count;
    if( debug_mode == 1 || average * 4 > 0.01 )
        disp(['estimated genotyping error: ', num2str(average*400, '%.2f'), '%']);
    end
end


end