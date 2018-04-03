
function marker_list = extract_markerlist(marker_allchr, chr)
    temp = marker_allchr(:,1) == chr;
    marker_list = marker_allchr(temp, 1:2);
end