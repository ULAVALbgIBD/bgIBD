function [int index1 index2 error] = merge2intervals(int1, int2)

error = 0;
int = [];
index1 = [];
index2 = [];


[r1, c1] = size(int1);
[r2, c2] = size(int2);

if( r1 <= 0 || r2 <= 0 || c1 ~= 2 || c2 ~= 2 )
    error = 1;
    disp('error in chromsome segmentation');
    return;
end

break1(1:r1-1) = 0;
break2(1:r2-1) = 0;
first = 1;
if( first ~= int1(1,1) || first ~= int2(1,1) )
    error = 1;
    disp('error in chromosome segmentation');
    return;
end
last = int1(r1,2);
if( last ~= int2(r2,2) )
    error = 1;
    disp('error in chromosome segmentation');
    return;
end

for i = 1:r1-1
    break1(i) = int1(i,2) + 0.5;
    if( break1(i) ~= int1(i+1,1) - 0.5 )
        error = 1;
        disp('error in chromosome segmentation');
        return;
    end
end

for i = 1:r2-1
    break2(i) = int2(i,2) + 0.5;
    if( break2(i) ~= int2(i+1,1) - 0.5 )
        error = 1;
        disp('error in chromosome segmentation');
        return;
    end
end

breaks = unique([break1,break2]);
breaks = sort(breaks, 'ascend');

len = length( breaks );

if( len < r1 - 1 || len < r2 - 1 )
    error = 1;
    disp('error in combining intervals');
    return;
end

if( len > 0 )
    if( breaks(1) < first + 0.5 || breaks(len) > last - 0.5 )
        error = 1;
        disp('error in combining intervals');
        return;
    end
end

int(1:len+1,1:2) = 0;
int(1,1) = first;
int(len+1,2) = last;

for i = 1:len
    int(i,2) = breaks(i) - 0.5;
    int(i+1,1) = breaks(i) + 0.5;
end

index1(1:len+1) = 0;
index2(1:len+1) = 0;

j = 1;
for i = 1:len + 1
    if( int(i,1) >= int1(j,1) && int(i,2) <= int1(j,2) )
        index1(i) = j;
    else
        j = j + 1;
        index1(i) = j;
    end
    if( j > r1 )
        error = 1;
        disp('error in combining intervals');
        return;
    end
end

for i = 1:len + 1
    j = index1(i);
    % whether the new interval is fully contained in the orginal one
    if( int(i,1) < int1(j,1) && int(i,2) > int1(j,2) )
        error = 1;
        disp('error in combining intervals')
        return;
    end
end

j = 1;
for i = 1:len + 1
    if( int(i,1) >= int2(j,1) && int(i,2) <= int2(j,2) )
        index2(i) = j;
    else
        j = j + 1;
        index2(i) = j;
    end
    if( j > r2 )
        error = 1;
        disp('error in combining intervals');
        return;
    end
end

for i = 1:len + 1
    j = index2(i);
    % whether the new interval is fully contained in the original one
    if( int(i,1) < int2(j,1) && int(i,2) > int2(j,2) )
        error = 1;
        disp('error in combining intervals');
        return;
    end
end


end




















