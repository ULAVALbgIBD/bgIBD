function [ibs error] = count_ibs(seq1, seq2)

ibs = [];
error = 0;

if( ndims(seq1) ~= 2 || ndims(seq2) ~= 2 )
    error = 1;
    return;
end


[r c] = size(seq1);
if( r <= 0 || c ~= 2 )
    error = 1;
    disp('dimension error');
    return;
end
len = r;

[r c] = size(seq2);
if( r ~= len || c ~= 2 )
    error = 1;
    disp('dimension error');
    return;
end

ibs = zeros(len,1);
s11 = seq1(1:len, 1);
s12 = seq1(1:len, 2);
s21 = seq2(1:len, 1);
s22 = seq2(1:len, 2);

missing(1:len) = ( (s11 == 0) | (s12 == 0) | (s21 == 0) | (s22 == 0) );

phase1(1:len) = ( (s11 == s21) + (s12 == s22) );
phase2(1:len) = ( (s11 == s22) + (s12 == s21) );
ibs(1:len) = max(phase1, phase2);
ibs(missing) = 3;


end
