function ibd = ibd_map_siblingpair(i1,i2,j1,j2)

ibd0 = 1;
ibd_l = 2;
ibd_r = 3;
ibd_lr = 4;



if( i1 == j1 )
    if( i2 == j2 )
        ibd = ibd_lr;
    else
        ibd = ibd_l;
    end
else
    if( i2 == j2 )
        ibd = ibd_r;
    else
        ibd = ibd0;
    end
end

end