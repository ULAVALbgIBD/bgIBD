
function temp5 = extrapolate(intervals, markers)
            temp5 = intervals;
            for j = 1:length(temp5(:,1))
                if( j == 1 )
                    if( temp5(j,1) ~= 1 )
                        temp5(j,1) = 1;
                    end
                end
                if( j == length(temp5(:,1)) )
                    if( temp5(j,2) ~= markers )
                        temp5(j,2) = markers;
                    end
                end
                if( j < length(temp5(:,1)) )
                    if( temp5(j+1,1) - temp5(j,2) > 1 )
                        mid = (temp5(j+1,1)+temp5(j,2))/2;
                        mid = floor(mid);
                        temp5(j,2) = mid;
                        temp5(j+1,1) = mid + 1;
                    end
                end              
            end
end