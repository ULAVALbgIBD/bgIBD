function hotspots = extract_hotspots(hotspots_allchr, chr)

hotspots = [];
if( chr <= 0 || chr > 24 )
    return;
end

[num_spots, c] = size( hotspots_allchr );
if( num_spots <= 0 || c ~= 5 )
    return;
end

hotspots = hotspots_allchr(hotspots_allchr(1:num_spots,1) == chr, 3:4);


end
