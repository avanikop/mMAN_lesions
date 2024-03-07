%% two birds repeat dist
%% bird5
fileID = fopen('D:\analysis\data_for_elife_mMAN\Source_data_1\bird5_prelesion','r');
b5_pre = fscanf(fileID, '%s')
fclose(fileID);
fileID = fopen('D:\analysis\data_for_elife_mMAN\Source_data_1\bird5_postlesion','r');
b5_post = fscanf(fileID, '%s')
fclose(fileID);

%

[startIndex,endIndex]=regexp(b5_pre,'cf{1,}');
len_f_repeats=(endIndex-startIndex);
[startIndex2,endIndex2]=regexp(b5_post,'cf{1,}');
len_f_repeats2=(endIndex2-startIndex2);

%% bird7 

fileID = fopen('D:\analysis\data_for_elife_mMAN\Source_data_1\bird7_prelesion','r');
b7_pre = fscanf(fileID, '%s')
fclose(fileID);
fileID = fopen('D:\analysis\data_for_elife_mMAN\Source_data_1\bird7_postlesion','r');
b7_post = fscanf(fileID, '%s')
fclose(fileID);

% counting repeats
[startIndex,endIndex]=regexp(b7_pre,'glb+');
len_b_repeats=(endIndex-startIndex)-1;
[startIndex2,endIndex2]=regexp(b7_post,'glb+');
len_b_repeats2=(endIndex2-startIndex2)-1;


%%
figure
tiledlayout(2,1)

sm=3
% Top plot
nexttile
bins = 0:2:40;
bincenters = 1:2:39;
newf=histcounts(len_f_repeats,bins,'Normalization','probability')
newf2=histcounts(len_f_repeats2,bins,'Normalization','probability')
plot(bincenters,smooth(newf,sm))
hold on
plot(bincenters,smooth(newf2,sm))
title('bird5 Length of f-repeats (normalized)','FontSize',16)
legend('Prelesion','Postlesion')
ylim([0 0.4])
box off

% Bottom plot
nexttile
bins = 0:1:10;
bincenters = 0.5:1:9.5;

newf=histcounts(len_b_repeats,bins,'Normalization','probability')
newf2=histcounts(len_b_repeats2,bins,'Normalization','probability')
plot(bincenters,smooth(newf,sm))
hold on
plot(bincenters,smooth(newf2,sm))
title('bird7 Length of b-repeats (normalized)','FontSize',16)
legend('Prelesion','Postlesion')
ylim([0 0.4])
box off
