function str = insertCommas(num)
    str = num2str(num);
    FIN = length(str); 
    for i = FIN-2:-3:2
    str(i+1:end+1) = str(i:end);
    str(i) = ',';
end