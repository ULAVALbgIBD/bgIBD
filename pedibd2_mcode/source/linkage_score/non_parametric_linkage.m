

function [output, error] = non_parametric_linkage(assignment, family)

error = 0;
output = [];

[alleles_all, ~, error] = check_allsegments(assignment, family);
if( error ~= 0 )
    disp('error in global IBD');
    return;
end


disp(' ');
disp('computing nonparametric linkage score ...');

time = cputime;

[z_score, p_value, error] = non_parametric_additive(alleles_all, family.structure);

output.z_score = z_score;
output.p_value = p_value;

display(['nonparametric linkage computing time: ', num2str(cputime - time), ' seconds']);

end

