function sampled_af = minor_allele_frq(af, smooth_interval)



% convert to minor and major allele frequency
[nmarkers, c] = size(af);
if( c ~= 3 )
    disp('error in allele frequency');
    return;
end

major(1:nmarkers) = 0;
minor(1:nmarkers) = 0;
for i = 1:nmarkers
    major(i) = max(af(i,1:2));
    minor(i) = min(af(i,1:2));
end
    
smoothed(1:nmarkers,1) = smooth(major, smooth_interval);
smoothed(1:nmarkers,2) = smooth(minor, smooth_interval);
smoothed(1:nmarkers,3) = smooth(af(1:nmarkers,3), smooth_interval);

sampled_af = smoothed;

for i = 1:nmarkers
    if( af(i,1) >= af(i,2) )
        sampled_af(i,1) = smoothed(i,1);
        sampled_af(i,2) = smoothed(i,2);
    else
        sampled_af(i,1) = smoothed(i,2);
        sampled_af(i,2) = smoothed(i,1);
    end
    sampled_af(i,3) = smoothed(i,3);
end


end