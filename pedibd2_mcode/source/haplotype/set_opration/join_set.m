function [ consistency ] = join_set( index1, index2, offset, ibd_link )


update_ancestor(index1);
update_ancestor(index2);

[anc1, off1, ans1] = find_immediate(index1);
[anc2, off2, ans2] = find_immediate(index2);

[assign1] = find_assign(anc1);
[assign2] = find_assign(anc2);


if( anc1 == anc2 )
    if( mod( off1+off2+offset, 2 ) ~= 0 )
        consistency = 0;
    else
        consistency = 1;
    end
else
    if( assign1 == 0 )
        consistency = 1;
    else
        if( assign2 == 0 )
            consistency = 1;
        else
            %both fixed
            if( mod(off1+off2+offset,2) == 0 )
                if( assign1 == assign2 )
                    consistency = 1;
                else
                    consistency = 0;
                end
            else
                if( assign1 ~= assign2 )
                    consistency = 1;
                else
                    consistency = 0;
                end
            end
        end
    end
    if( ibd_link || consistency == 1 )
        set_immediate(anc1, anc2, mod(off1+off2+offset, 2));
    end
end

update_ancestor(index1);
update_ancestor(index2);

end

