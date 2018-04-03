function [ conflict ] = update_ancestor( index )


conflict = 1;

[anc, off, a] = find_immediate(index);

if( anc ~= find_immediate(anc) )
    update_ancestor( anc );
    [new_anc, add_off, temp] = find_immediate(anc);
    set_immediate(index, new_anc, mod(add_off+off,2));
   
end

    [anc, off, temp] = find_immediate(index);
    [ans1, ans2, a] = find_immediate(anc);
    if( a ~= 0 )
        [conflict] = set_assign(index, a, off);
    end 

end

