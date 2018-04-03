function af = get_afreq(freq, default)

sum = freq(1) + freq(2);
if( sum == 0 )
    af = default;
    % completely missing locus
    % no observation data, this locus does not influence decoding
else
    af(1) = freq(1)/(freq(1)+freq(2));
    af(2) = freq(2)/(freq(1)+freq(2));
    af(3) = freq(3)/(freq(1)+freq(2)+freq(3));
end

