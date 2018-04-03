function epl = sc_generate_emission_all_loci(all_af, ms, te, markerlist)

epl = [];

[nloc, fields] = size(markerlist);
if( nloc <= 0 || fields ~= 2 )
    disp('error in marker list');
    return;
end

[r, c] = size(all_af);
if( r <= 0 || r ~= nloc || c ~= 3 )
    disp('error in allele frequency');
    return;
end


for i = 1:nloc
    temp = sc_generate_emission(all_af(i, :), ms, te);
    temp(3,:) = temp(2,:);
    % background ibd emission == foreground ibd emission
    epl{i} = temp;
end

end


