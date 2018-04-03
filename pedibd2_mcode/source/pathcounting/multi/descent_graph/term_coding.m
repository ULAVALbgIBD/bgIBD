function code = term_coding(id)
    
    % id is a vector, indicating which slots of alleles are being put on
    % recurrence
    % mapped to a single numeric value, for easy lookup
    code = 0;
    for i = 1:length(id)
        code = code + 2^(id(i)-1);
    end
    
end