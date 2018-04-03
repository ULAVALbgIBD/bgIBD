

%%

segment = 0;
for i = 1:length(plink_pvalue)
    if( plink_pvalue(i, 1) > segment )
        segment = plink_pvalue(i, 1);
        temp(segment) = plink_pvalue(i, 2);
    else
        if( temp(plink_pvalue(i, 1)) > plink_pvalue(i, 2) )
            temp(plink_pvalue(i, 1)) = plink_pvalue(i, 2);
        end
    end
end

%%
plot(temp, '.');
xlim([0, length(temp)]);


%%
% segment haplotype solution into segments, with new coding

count2 =stat(JL_all_chr3_tdt);
%count_ignore_normal = stat(JL_all_chr3_tdt);

%%

%generate p value (anova test) for each segment

for i = 5565:5565
    p(i) = fstat(count{i});
    p2(i) = fstat(count2{i});
    %p_in(i) = fstat(count_ignore_normal{i});
end

%%

plot(p2, '.');
xlim([0, 6510]);


%%

fid = fopen('G:\Paper\ZRHC_Missing\ZRHC\ZRHC\data\JL_Bfile\result\p_tdt_nbn.txt', 'w');
for i = 1:6510
    fprintf(fid, '%d\t%f\n', i, p2(i));
end
fclose(fid);

%%

%%

load 'G:\Paper\ZRHC_Missing\ZRHC\ZRHC\data\JL_Bfile\result\JL_all_fam.txt';
%load 'G:\Paper\ZRHC_Missing\ZRHC\ZRHC\data\JL_Bfile\result\JL_all_chr3_linkage.txt';
fid = fopen('G:\Paper\ZRHC_Missing\ZRHC\ZRHC\data\JL_Bfile\result\JL_all_chr3_linkage.txt', 'r');
linkage = textscan(fid, '%d %d %d %d %d %d %d %d %d %d %d %d %*[^\n]');
fclose(fid);


%%

JL_linkage = [];

for i = 1:length(linkage)
    JL_linkage = [JL_linkage, linkage{i}];
end

%%
% put recombination from different families together
% organized by chromosome

for i = 1:22
    JL_chr_recom{i} = [];
end

for j = 1:length(linkage{1})
    JL_chr_recom{linkage{2}(j)} = union(JL_chr_recom{linkage{2}(j)}, linkage{4}(j)); 
end

%%

c = 3;

f = 0;  %family id

count_f = 0;
count_id = 0;
count_f_4 = 0;
count_f_3 = 0;

%JL_fam_ext{c} = [];

%JL_fam_ext{c}(:, 1:6) = JL_all_fam(:, 1:6);

for f = 1:1000
    


        temp4 = find(JL_all_fam(:, 1) == f);
        
        if( isempty(temp4) )
            continue;
        else
            count_f = count_f + 1;
            count_id = count_id + length(temp4);
            if( length(temp4) >= 4 )
                count_f_4 = count_f_4 + 1;
                JL_ds(count_f_4, 1:2) = JL_all_fam(temp4(3:4), 6);
            end
            if( length(temp4) == 3 )
                count_f_3 = count_f_3 + 1;
            end  
        end
        
        continue;
   
        temp1 = find( JL_linkage(:, 2) == c );
        temp2 = find( JL_linkage(temp1, 1) == f );
        temp3 = temp1(temp2);


        for a = 1:length(JL_chr_recom{c})
            for b = 1:length(temp3)
                if( JL_linkage(temp3(b), 3) <= JL_chr_recom{c}(a) && JL_chr_recom{c}(a) <= JL_linkage(temp3(b), 4))
                    for count = 1:length(temp4)
                        JL_fam_ext{c}(temp4(count), 6 + 2 * a - 1) = JL_linkage(temp3(b), 4 + 2 * count - 1);
                        JL_fam_ext{c}(temp4(count), 6 + 2 * a) = JL_linkage(temp3(b), 4 + 2 * count);
                    end
                end                   
            end          
        end

end

%%

% rename individual ids, to avoid duplication

max_id = max(max(JL_all_fam(:,2:4)));
JL_id_map(1:max_id) = 0;
count_id = 0;

for i = 1:length(JL_all_fam(:,2))
    for j = 2:4
        if( JL_all_fam(i,j) == 0 )
            continue;
        end
        if( JL_id_map(JL_all_fam(i,j)) == 0 )
            count_id = count_id + 1;
            JL_id_map(JL_all_fam(i,j)) = count_id;
        end
        JL_fam_ext{3}(i,j) = JL_id_map(JL_all_fam(i,j));
    end
end


%%

fid = fopen('G:\Paper\ZRHC_Missing\ZRHC\ZRHC\data\JL_Bfile\result\JL_all_fam_data.txt', 'w');

fprintf(fid, 'FAMILY\tIND\tFATHER\tMOTHER\tSEX\tTRAIT\t');
count = 0;
for j = 3
    fprintf(fid, '%d\t', count + 1: count + (length(JL_fam_ext{j})-6)/2 );
    count = count + (length(JL_fam_ext{j})-6)/2;
    fprintf(fid, '\n');

    for k = 1:(length(JL_fam_ext{j}(:,1)))
        temp = find( JL_fam_ext{j}(:, 1) == JL_fam_ext{j}(k, 1) );
%         if( length(temp) <= 3 || JL_fam_ext{j}(k, 1) > 100 )
%             k = k + length(temp) - 1;
%             continue;
%         end
        fprintf(fid, '%d\t', JL_fam_ext{j}(k, 1:6));
        for p = 7:2:length(JL_fam_ext{j}(k,:))
            fprintf(fid, '%d/%d\t', JL_fam_ext{j}(k, p), JL_fam_ext{j}(k, p+1));
        end
        fprintf(fid, '\n');
    end
    
end
fclose(fid);

%%

fid = fopen('G:\Paper\ZRHC_Missing\ZRHC\ZRHC\data\JL_Bfile\result\JL_marker.txt', 'w');

count = 0;
for j = 3
    for k = count + 1: count + (length(JL_fam_ext{3}(1,:))-6)/2
        fprintf(fid, '%d\n1\n2\n3\n4\n;;\n', k);
    end
    count = count + (length(JL_fam_ext{3}(1,:))-6)/2;
end



fclose(fid);

%%

load 'G:\Paper\ZRHC_Missing\ZRHC\ZRHC\data\JL_Bfile\result\JL_chr3_p.txt'
load 'G:\Paper\ZRHC_Missing\ZRHC\ZRHC\data\JL_Bfile\result\JL_chr3_m.txt'
load 'G:\Paper\ZRHC_Missing\ZRHC\ZRHC\data\JL_Bfile\result\chr3_markerlist.txt'

%%
map = chr3_markerlist(:,2);



%%

% recombination rate

% 112 families, each has 4 miosis, 170Mbps

[n, xout] = hist(map([p;m]), 17*2);;
line(xout, n*100/(5*112*4));


%%


all_rec = sort([(chr1_p(:,1)+chr1_p(:,2))./2; (chr1_m(:,1)+chr1_m(:,2))./2]);
for i = 1:length(all_rec)-1
    all_rec_int(i) = all_rec(i+1) - all_rec(i);
end
cdfplot((all_rec_int));

%%

%plot((all_rec_int(1:100)), '.');
bar(all_rec_int(1:100));

%%



JL_chr6_p_bp = [map(JL_chr3_p(:, 1)), map(JL_chr3_p(:, 2))];
JL_chr6_m_bp = [map(JL_chr3_m(:, 1)), map(JL_chr3_m(:, 2))];


%%

hold on;
JL_chr3_resolution = [JL_chr3_p_bp(:,2)-JL_chr3_p_bp(:,1); JL_chr3_m_bp(:,2)-JL_chr3_m_bp(:,1)];
%cdfplot((JL_chr3_resolution));
temp = find(JL_chr3_resolution < 100000);
hist(JL_chr3_resolution(temp), 30);



%%

load 'G:\Paper\ZRHC_Missing\ZRHC\ZRHC\data\ucsd\whole_genome\Tao\result\stat\hotspots\new_hotspots.txt';


%%

count = 0;
for i = 1:length(output.paths.output_breakpoint)
    temp(count + 1:count + length(output.paths.output_breakpoint{i})) = parameters.sampled_markerlist(output.paths.output_breakpoint{i},2);
    count = count + length(output.paths.output_breakpoint{i});
end

temp1 = find(new_hotspots(:,1) == 22);

output.paths.hotspot_overlap = hotspot_overlap(temp, new_hotspots(temp1,2));

sum(new_hotspots(temp1,5))

plot_recombination(new_hotspots(temp1,2));

clear temp1 temp i count;

%%

%association with copy number variation

JL_chr6_bp = [JL_chr6_p_bp;JL_chr6_m_bp];



%%

%permutation to get a p value

perm_times = 1000;

nz = [];

for i = 1:perm_times

JL_chr6_cnv_overlap = hotspot_overlap_perm(JL_chr6_bp, chr6_cnv(:, 1:2), map);

nz(i) = length(find(JL_chr6_cnv_overlap == 1));

end

%%

JL_chr6_cnv_overlap = hotspot_overlap(JL_chr6_bp, chr6_cnv(:, 1:2));

a_nz = length(find(JL_chr6_cnv_overlap == 1));

p_value_cnv = length(find(nz >= a_nz))/perm_times;

%%

overlap_rec = [];

j = 0;
for i = 1:length(JL_chr6_cnv_overlap)
    if( JL_chr6_cnv_overlap(i) ~= 0 )
        j = j + 1;
        overlap_rec(j,1:2) = JL_chr6_bp(i,1:2);
    end
end

%%
JL_chr6_p_cnv_overlap = hotspot_overlap(JL_chr6_p_bp, chr6_cnv(:, 1:2));
JL_chr6_m_cnv_overlap = hotspot_overlap(JL_chr6_m_bp, chr6_cnv(:, 1:2));

%%

coverage = sum(chr6_cnv(:,2)-chr6_cnv(:,1))/(max(map)-min(map));

%%

%recombination interval p
%individual by individual

p = length(JL_map);
JL_interval = [];
inc = 0;
for i = 1:length(chr1_p_bp)
    temp = mean(chr1_p_bp(i, 1:2));
    if( temp > p )
        inc = inc + 1;
        JL_interval(inc) = temp - p;
    end
    p = temp;
end


%%

plot_interval(JL_interval, 0);

%%

load 'G:\Paper\ZRHC_Missing\ZRHC\ZRHC\data\JL_Bfile\result\JL_all_chr3_hf_s.txt';

%%
load 'G:\Paper\ZRHC_Missing\ZRHC\ZRHC\data\JL_Bfile\result\JL_all_chr3_hf_s_oldpool.txt';


%%

hn = JL_all_chr3_hf_s_oldpool(:,2);

%%

hn = nf_oldpool;

%%
hn = JL_all_chr3_hf_s(:,2);

%%
hn = nf; %number covering 0.9 of frequency
%%
count = 0;
step = 47;
haplotype = [];
for i = 1:step:length(hn)-step
    count = count + 1;
    haplotype(count, 1) = sum(hn(i:i+step-1))/(step);
    haplotype(count, 2) = map((i+step/2-1)*20+1);
end

%%

line(haplotype(:,2), haplotype(:,1), 'Color','r', 'LineStyle', '-');

%%
%distribution of haplotype diversity

plot(JL_all_chr3_hf_s(:,2));

%%



load 'G:\Paper\ZRHC_Missing\ZRHC\ZRHC\data\JL_Bfile\result\JL_all_chr3_hf_oldpool.txt';

%%

JL_all_chr3_hf_c_oldpool = cummulative(JL_all_chr3_hf_oldpool);

%%

for i = 1:1627
    nf_oldpool(i) = find( JL_all_chr3_hf_c_oldpool(i,:) > 0.8*max(JL_all_chr3_hf_c_oldpool(i,:)), 1 );
end

%%
load 'G:\Paper\ZRHC_Missing\ZRHC\ZRHC\data\JL_Bfile\result\JL_all_chr3_hf.txt';


%%
for i = 1:1627
    nf(i) = find( JL_all_chr3_hf_c(i,:) > 0.9, 1 );
end

%%

%distribution of haplotype frequency
temp = find( JL_all_chr3_hf(:, 1) == 22 );
cdfplot( JL_all_chr3_hf(temp, 2) );

%%

%cummulative haplotype frequency distribution for each segment i

JL_all_chr3_hf_c = cummulative(JL_all_chr3_hf);

%%

%for each haplotype top x percent haplotype coverage of total frequency

for i = 1:1627
    threshold = int32(JL_all_chr3_hf_s(i, 2) * 0.5);
    JL_all_chr3_hf_c_p(i, 1) = JL_all_chr3_hf_c(i, threshold);
end
plot( JL_all_chr3_hf_c_p(:, 1) );


%%

load 'G:\Paper\ZRHC_Missing\ZRHC\ZRHC\data\JL_Bfile\data\freedom\all_chr3_2.txt_IBS.txt'

%%
load 'G:\Paper\ZRHC_Missing\ZRHC\ZRHC\data\JL_Bfile\data\freedom\3all_chr3_2.txt_IBS_familyof3.txt'

%%

load 'G:\Paper\ZRHC_Missing\ZRHC\ZRHC\data\JL_Bfile\data\freedom\all_chr3_2.txt_IBS_hetero.txt'

%%

data = X3all_chr3_2(:,2:32551);
%data = all_chr3_2(:, 3:32552);
%deno = sum(all_chr3_2(:,2));
deno = 81;
%deno = 112;
freedom = sum(data);
count = 0;
step = 1000;
freedom_smooth = [];
for i = 1:step:length(freedom)-step
    count = count + 1;
    freedom_smooth(count, 1) = sum(freedom(i:i+step-1))/(step*deno);
    freedom_smooth(count, 2) = map(i+step/2);
end

%scatter(freedom_smooth(:,2), freedom_smooth(:,1), '-x');
line(freedom_smooth(:,2), freedom_smooth(:,1), 'Color','b');

%%

%marker density
[n, xout] = hist(map, 170);
line(xout, n);


%%









