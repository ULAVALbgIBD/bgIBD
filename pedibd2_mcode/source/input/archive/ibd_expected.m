function de_ibd = ibd_expected(ibs_data, epl)



% 0: ibs0
% 1: ibs1
% 2: ibs2
% 3: missing genotype

ibs0 = 1;
ibs1 = 2;
ibs2 = 3;
missing = 4;

ibd0 = 1;
ibd1 = 2;
ibd2 = 3;

Nmarkers = length(ibs_data);
ce(1:3,1:4) = 0; 
% cummulated emmission

for i = 1:Nmarkers
    ce = ce + epl{i};
end

ce = ce';

N(ibs0) = nnz(ibs_data == 0);
N(ibs1) = nnz(ibs_data == 1);
N(ibs2) = nnz(ibs_data == 2);
N(missing) = nnz(ibs_data == 3);

N = N';

% N = N(1:3);
% ce = ce(1:3,1:3);
% % de_ibd = (ce^(-1))*N;
% 
% temp = ((ce'*ce)^(-1))*(ce');
% count = 0;
% de_ibd(1:3,1:1000) = 0;
% for i = 1:1000
%     dt1 = 0.05 - rand()/10;
%     dt2 = 0.05 - rand()/10;
%     a = sum(N);
%     n(1) = N(1) * (1 + dt1);
%     n(2) = N(2) * (1 + dt2);
%     n(3) = a - N(1) - N(2);
%     de_ibd(:,i) = temp * n';
% end

de_ibd = ((ce'*ce)^(-1))*(ce'*N);

% temp = ce * [-1,-1;1,-2;0,1];
% de_ibd = ((temp'*temp)^(-1))* ( temp'*(N - ce*[1;0;0]) );

