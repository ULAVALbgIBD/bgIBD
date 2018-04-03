function [result] = compare_haplotype(ref_genotype, output_genotype, input_genotype, i)

            % currently valid for only one individual
            
            bitIND = false(size(ref_genotype,2),1);
            bitIND(i) = true;
            
            input_NOmissing = input_genotype(:,:,1) ~= 0 ...
                & input_genotype(:,:,2) ~= 0;
            input_homo = input_NOmissing ...
                & input_genotype(:,:,1) == input_genotype(:,:,2);
            input_hetero = input_NOmissing ...
                & input_genotype(:,:,1) ~= input_genotype(:,:,2);
            
            output_NOmissing = output_genotype(:,:,1) ~= 0 ...
                & output_genotype(:,:,2) ~= 0;
            output_homo = output_NOmissing ...
                & output_genotype(:,:,1) == output_genotype(:,:,2);
            output_hetero = output_NOmissing ...
                & output_genotype(:,:,1) ~= output_genotype(:,:,2);
            output_correct1 = output_hetero ...
                & all(output_genotype == ref_genotype, 3);
            output_correct2 = output_hetero ...
                & all(output_genotype == ref_genotype(:,:,[2,1]), 3);
            output_correct = output_NOmissing ...
                & (all(output_genotype == ref_genotype, 3) ...
                | all(output_genotype == ref_genotype(:,:,[2,1]), 3));
            
            
            ref_NOmissing = ref_genotype(:,:,1) ~= 0 ...
                & ref_genotype(:,:,2) ~= 0;
            
            result(1) = nnz(input_hetero(:,bitIND));
            result(2) = nnz(input_homo(:,bitIND));
            result(3) = nnz(input_NOmissing(:,bitIND));       
            
            result(4) = nnz(output_NOmissing(:,bitIND));
            result(5) = nnz(output_hetero(:,bitIND)); % output hetero
            result(6) = max( nnz(output_correct1(:,bitIND)), ...
                nnz(output_correct2(:,bitIND)) ); % output INCORRECT hetero
            
            result(7) = nnz(output_correct(:,bitIND));
            
end





