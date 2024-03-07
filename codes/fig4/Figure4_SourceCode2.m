%% code for making fig4D cv of repeat number
%% bird3
fileID = fopen('D:\analysis\data_for_elife_mMAN\Source_data_1\bird3_prelesion','r');
b3_pre = fscanf(fileID, '%s')
fclose(fileID);
fileID = fopen('D:\analysis\data_for_elife_mMAN\Source_data_1\bird3_postlesion','r');
b3_post = fscanf(fileID, '%s')
fclose(fileID);
%% bird4

fileID = fopen('D:\analysis\data_for_elife_mMAN\Source_data_1\bird4_prelesion','r');
b4_pre = fscanf(fileID, '%s')
fclose(fileID);
fileID = fopen('D:\analysis\data_for_elife_mMAN\Source_data_1\bird4_postlesion','r');
b4_post = fscanf(fileID, '%s')
fclose(fileID);

%% bird5
fileID = fopen('D:\analysis\data_for_elife_mMAN\Source_data_1\bird5_prelesion','r');
b5_pre = fscanf(fileID, '%s')
fclose(fileID);
fileID = fopen('D:\analysis\data_for_elife_mMAN\Source_data_1\bird5_postlesion','r');
b5_post = fscanf(fileID, '%s')
fclose(fileID);
%% bird6

fileID = fopen('D:\analysis\data_for_elife_mMAN\Source_data_1\bird6_prelesion','r');
b6_pre = fscanf(fileID, '%s')
fclose(fileID);
fileID = fopen('D:\analysis\data_for_elife_mMAN\Source_data_1\bird6_postlesion','r');
b6_post = fscanf(fileID, '%s')
fclose(fileID);
%% bird7 

fileID = fopen('D:\analysis\data_for_elife_mMAN\Source_data_1\bird7_prelesion','r');
b7_pre = fscanf(fileID, '%s')
fclose(fileID);
fileID = fopen('D:\analysis\data_for_elife_mMAN\Source_data_1\bird7_postlesion','r');
b7_post = fscanf(fileID, '%s')
fclose(fileID);

%% bird 3
% b repeats

[startIndex,endIndex]=regexp(b3_pre,'[b]+');
b3_r1_pre=(endIndex-startIndex)+1;
[startIndex2,endIndex2]=regexp(b3_post,'[b]+');
b3_r1_po=(endIndex2-startIndex2)+1;

CV_pre{1} = nanstd(b3_r1_pre)/nanmean(b3_r1_pre);
CV_po{1} = nanstd(b3_r1_po)/nanmean(b3_r1_po);

% c repeats

[startIndex,endIndex]=regexp(b3_pre,'[c]+');
b3_r2_pre=(endIndex-startIndex)+1;
[startIndex2,endIndex2]=regexp(b3_post,'[c]+');
b3_r2_po=(endIndex2-startIndex2)+1;

CV_pre{2} = nanstd(b3_r2_pre)/nanmean(b3_r2_pre);

CV_po{2} = nanstd(b3_r2_po)/nanmean(b3_r2_po);


%% bird4

[startIndex,endIndex]=regexp(b4_pre,'ddb+');
b4_r1_pre=(endIndex-startIndex)-1;
[startIndex2,endIndex2]=regexp(b4_post,'ddb+');
b4_r1_post=(endIndex2-startIndex2)-1;

CV_pre{3} = nanstd(b4_r1_pre)/nanmean(b4_r1_pre);

CV_po{3} = nanstd(b4_r1_post)/nanmean(b4_r1_post);

%% bird5

% length of f repeats
[startIndex,endIndex]=regexp(b5_pre,'cf{1,}');
b5_r1_pre=(endIndex-startIndex);
[startIndex2,endIndex2]=regexp(b5_post,'cf{1,}');
b5_r1_po=(endIndex2-startIndex2);

CV_pre{4} = nanstd(b5_r1_pre)/nanmean(b5_r1_pre);

CV_po{4} = nanstd(b5_r1_po)/nanmean(b5_r1_po);


%% bird6
% mcbbb
[startb endb]=regexp(b6_pre,'mc[ba]+');
b6_r1_pre=endb-startb+1;
[startb2 endb2]=regexp(b6_post,'mc[ba]+');
b6_r1_po=endb2-startb2+1;

CV_pre{5} = nanstd(b6_r1_pre)/nanmean(b6_r1_pre);

CV_po{5} = nanstd(b6_r1_po)/nanmean(b6_r1_po);

%% bird7

% counting b+ 

[startIndex,endIndex]=regexp(b7_pre,'glb+');
b7_r1_pre=(endIndex-startIndex)-1;
[startIndex2,endIndex2]=regexp(b7_post,'glb+');
b7_r1_post=(endIndex2-startIndex2)-1;


CV_pre{6} = nanstd(b7_r1_pre)/nanmean(b7_r1_pre);

CV_po{6} = nanstd(b7_r1_post)/nanmean(b7_r1_post);


%% coefficient of variation
% stddev/mean it shows the extent of variability of the distribution.
% Higher the Cv, greater the dispersion

CV_all=[CV_pre{:}; CV_po{:}];
CV_all=CV_all';

%%
figure
markers=['x','x','s','d','^','h'];

boxplot(CV_all,'Whisker',1,'Colors','br')
hold on
for ip=1:length(CV_all)
    plot(1,CV_all(ip,1),'Marker',markers(ip),'MarkerSize',12,'MarkerEdgeColor','b');

    plot(2,CV_all(ip,2),'Marker',markers(ip),'MarkerSize',12,'MarkerEdgeColor','r');
end

for m=1:length(CV_all)
    plot([1,2], [CV_all(m,1),CV_all(m,2)],'Color',[.7 .7 .7], 'LineWidth', 0.5)

end
ax=gca
ax.FontSize=16
xticklabels({'Prelesion','Postlesion'})
%ylim([0,1])
%yticklabels([0:0.1:1])
ylabel('Coefficient of variation','FontSize',16)
title('Coefficient of variation for repeat distributions (n=6)','FontSize',16)
box off
hold off
daspect([6 1 1])

% testing cv
[p,h,stats] = signrank(CV_all(:,1),CV_all(:,2),'alpha',0.05)
