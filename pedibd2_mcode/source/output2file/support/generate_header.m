function header = generate_header(pf, gf, mf)

p = parse_header(pf);
g = parse_header(gf);
m = parse_header(mf);

len1 = length(p);
if( len1 <= 3 )
    sp = ['P', p];
else
    sp = ['P',p(len1-2:len1)];
end
len2 = length(g);
if( len2 <= 3 )
    sg = ['G', g];
else
    sg = ['G', g(len2-2:len2)];
end
len3 = length(m);
if( len3 <= 3 )
    sm = ['M', m];
else
    sm = ['M', m(len3-2:len3)];
end

name = [sp,sg,sm];
name = p;
if( length(name) <= 0 )
    name = 'dbF';
end

header = [pwd, '/', name];

end