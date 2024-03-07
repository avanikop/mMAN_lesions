% code for plotting fig4 C
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

% c repeats

[startIndex,endIndex]=regexp(b3_pre,'[c]+');
b3_r2_pre=(endIndex-startIndex)+1;
[startIndex2,endIndex2]=regexp(b3_post,'[c]+');
b3_r2_po=(endIndex2-startIndex2)+1;


%% bird4

[startIndex,endIndex]=regexp(b4_pre,'ddb+');
b4_r1_pre=(endIndex-startIndex)-1;
[startIndex2,endIndex2]=regexp(b4_post,'ddb+');
b4_r1_post=(endIndex2-startIndex2)-1;


%% bird5

% length of f repeats
[startIndex,endIndex]=regexp(b5_pre,'cf{1,}');
b5_r1_pre=(endIndex-startIndex);
[startIndex2,endIndex2]=regexp(b5_post,'cf{1,}');
b5_r1_po=(endIndex2-startIndex2);


%% bird6
% mcbbb
[startb endb]=regexp(b6_pre,'mc[ba]+');
b6_r1_pre=endb-startb+1;
[startb2 endb2]=regexp(b6_post,'mc[ba]+');
b6_r1_po=endb2-startb2+1;


%% bird7

% counting b+ 

[startIndex,endIndex]=regexp(b7_pre,'glb+');
b7_r1_pre=(endIndex-startIndex)-1;
[startIndex2,endIndex2]=regexp(b7_post,'glb+');
b7_r1_post=(endIndex2-startIndex2)-1;


%% figure for average length of repeats

meanlengthspre=[nanmean(b3_r1_pre);nanmean(b3_r2_pre);nanmean(b4_r1_pre); nanmean(b5_r1_pre);...
    nanmean(b6_r1_pre);nanmean(b7_r1_pre)];
meanlengthspo=[nanmean(b3_r1_po);nanmean(b3_r2_po);nanmean(b4_r1_post); nanmean(b5_r1_po);...
    nanmean(b6_r1_po);nanmean(b7_r1_post)];
len_all=[meanlengthspre meanlengthspo];
% stats on length of repeats

[p,h,stats]=signrank(meanlengthspre,meanlengthspo)
% comparing repeat length dists with chisq
% just doing manually
[h1,p1]=kstest2(b3_r1_pre,b3_r1_po)
[h2,p2]=kstest2(b3_r2_pre,b3_r2_po)
[h3,p3]=kstest2(b4_r1_pre,b4_r1_post)
[h4,p4]=kstest2(b5_r1_pre,b5_r1_po)
[h5,p5]=kstest2(b6_r1_pre,b6_r1_po)
[h6,p6]=kstest2(b7_r1_pre,b7_r1_post)

%%
figure
markers=['x','x','s','d','^','h'];
boxplot(len_all,'Whisker',1,'Colors','br')
hold on
for il=1:length(len_all)
plot(1,len_all(il,1),'Marker',markers(il),'MarkerSize',12,'MarkerEdgeColor','b');
hold on
plot(2,len_all(il,2),'Marker',markers(il),'MarkerSize',12,'MarkerEdgeColor','r');
end

for m=1:length(len_all)
    plot([1,2], [len_all(m,1),len_all(m,2)],'Color',[.7 .7 .7], 'LineWidth', 0.5)
    hold on
end
ax=gca
ax.FontSize=16
xticklabels({'Prelesion','Postlesion'})
%ylim([0,1])
%yticklabels([0:0.1:1])
ylabel('Average length of repeats','FontSize',16)
title('Average length of repeats (n=6)','FontSize',16)
box off
hold off
daspect([1 6 1])