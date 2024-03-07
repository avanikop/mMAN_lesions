% for analysis of chunks 3d
% I look for Prelesion chunks within postlesion data
%%
% importing textfiles and reading them:
% replaced repeats files:
%% bird1
fileID = fopen('D:\analysis\data_for_elife_mMAN\Source_data_1\bird1_prelesion_REPR','r');
b1_pre_R = fscanf(fileID, '%s')
fclose(fileID);
fileID = fopen('D:\analysis\data_for_elife_mMAN\Source_data_1\bird1_postlesion_REPR','r');
b1_post_R = fscanf(fileID, '%s')
fclose(fileID);
%% bird2

fileID = fopen('D:\analysis\data_for_elife_mMAN\Source_data_1\bird2_prelesion_REPR','r');
b2_pre_R = fscanf(fileID, '%s')
fclose(fileID);
fileID = fopen('D:\analysis\data_for_elife_mMAN\Source_data_1\bird2_postlesion_REPR','r');
b2_post_R = fscanf(fileID, '%s')
fclose(fileID);
%% bird3
fileID = fopen('D:\analysis\data_for_elife_mMAN\Source_data_1\bird3_prelesion_REPR','r');
b3_pre_R = fscanf(fileID, '%s')
fclose(fileID);
fileID = fopen('D:\analysis\data_for_elife_mMAN\Source_data_1\bird3_postlesion_REPR','r');
b3_post_R = fscanf(fileID, '%s')
fclose(fileID);
%% bird4

fileID = fopen('D:\analysis\data_for_elife_mMAN\Source_data_1\bird4_prelesion_REPR','r');
b4_pre_R = fscanf(fileID, '%s')
fclose(fileID);
fileID = fopen('D:\analysis\data_for_elife_mMAN\Source_data_1\bird4_postlesion_REPR','r');
b4_post_R = fscanf(fileID, '%s')
fclose(fileID);

%% bird5
fileID = fopen('D:\analysis\data_for_elife_mMAN\Source_data_1\bird5_prelesion_REPR','r');
b5_pre_R = fscanf(fileID, '%s')
fclose(fileID);
fileID = fopen('D:\analysis\data_for_elife_mMAN\Source_data_1\bird5_postlesion_REPR','r');
b5_post_R = fscanf(fileID, '%s')
fclose(fileID);
%% bird6

fileID = fopen('D:\analysis\data_for_elife_mMAN\Source_data_1\bird6_prelesion_REPR','r');
b6_pre_R = fscanf(fileID, '%s')
fclose(fileID);
fileID = fopen('D:\analysis\data_for_elife_mMAN\Source_data_1\bird6_postlesion_REPR','r');
b6_post_R = fscanf(fileID, '%s')
fclose(fileID);
%% bird7 

fileID = fopen('D:\analysis\data_for_elife_mMAN\Source_data_1\bird7_prelesion_REPR','r');
b7_pre_R = fscanf(fileID, '%s')
fclose(fileID);
fileID = fopen('D:\analysis\data_for_elife_mMAN\Source_data_1\bird7_postlesion_REPR','r');
b7_post_R = fscanf(fileID, '%s')
fclose(fileID);

%% bird1 
[chunks2,~,seqforchunks,chunks2replace,labelidx]=seq_chunkextractionfunc(b1_pre_R,0);
close all

bird1{1,1}=chunkconsistency(seqforchunks,'5.....'); %pDpacc
%xy
bird1{2,1}=chunkconsistency(b1_pre_R,'x.');
%Yir
unq=unique(seqforchunks);
allletters=[char(97:122),char(65:90),'0123456789']; %in case some extra chars are needed, i added char 60:64
avletters1=allletters(~ismember(allletters,unq));
Ysym=avletters1(1);
seqforchunks2=regexprep(seqforchunks,'Yi',Ysym);
bird1{3,1}=chunkconsistency(seqforchunks2,[Ysym,'.']);


[chunks2,~,seqforchunks,chunks2replace,labelidx2]=seq_chunkextractionfunc(b1_post_R,0);
close all

% (r)pdddpacc
bird1{1,2}=chunkconsistency(seqforchunks,'2.....');

%xy
bird1{2,2}=chunkconsistency(b1_post_R,'x.');

%Yir
unq=unique(seqforchunks);
allletters=[char(97:122),char(65:90),'0123456789']; %in case some extra chars are needed, i added char 60:64
avletters1=allletters(~ismember(allletters,unq));
Ysym=avletters1(1);
seqforchunks2=regexprep(seqforchunks,'Yi',Ysym);
bird1{3,2}=chunkconsistency(seqforchunks2,[Ysym,'.']);

%% bird2

[chunks2,~,seqforchunks,chunks2replace,labelidx]=seq_chunkextractionfunc(b2_pre_R,0);
close all

bird2{1,1}=chunkconsistency(seqforchunks,'8.....');%ccllfk
bird2{2,1}=chunkconsistency(seqforchunks,'j..'); %jgd

[chunks2,~,seqforchunks,chunks2replace,labelidx2]=seq_chunkextractionfunc(b2_post_R,0);
close all

bird2{1,2}=chunkconsistency(seqforchunks,'4.....',6); %ccllfk have to give length because of [12]
bird2{2,2}=chunkconsistency(seqforchunks,'j..'); %jgd

%% bird3

[chunks2,~,seqforchunks,chunks2replace,labelidx]=seq_chunkextractionfunc(b3_pre_R,0);
close all

bird3{1,1}=chunkconsistency(seqforchunks,'4..'); %ec+b+
bird3{2,1}=chunkconsistency(seqforchunks,'a..'); %ahe

unq=unique(seqforchunks);
allletters=[char(97:122),char(65:90),'0123456789']; %in case some extra chars are needed, i added char 60:64
avletters1=allletters(~ismember(allletters,unq));
Ysym=avletters1(1);
seqforchunks2=regexprep(seqforchunks,'Yi',Ysym);
bird3{3,1}=chunkconsistency(seqforchunks2,[Ysym,'.']); %Yie5 

[chunks2,~,seqforchunks,chunks2replace,labelidx2]=seq_chunkextractionfunc(b3_post_R,0);
close all

bird3{1,2}=chunkconsistency(seqforchunks,'2..',3); %%% ec+b+! had to write length because of [Xf]
bird3{2,2}=chunkconsistency(seqforchunks,'[876]..',3); %ahe

unq=unique(seqforchunks);
allletters=[char(97:122),char(65:90),'0123456789']; %in case some extra chars are needed, i added char 60:64
avletters1=allletters(~ismember(allletters,unq));
Ysym=avletters1(1);
seqforchunks2=regexprep(seqforchunks,'Yi',Ysym);
bird3{3,2}=chunkconsistency(seqforchunks2,[Ysym,'.']);%Yie5 although not present in postlesion%% wh87or32
%% bird4

[chunks2,~,seqforchunks,chunks2replace,labelidx]=seq_chunkextractionfunc(b4_pre_R,0);
close all

bird4{1,1}=chunkconsistency(seqforchunks,'F.'); % Fc
bird4{2,1}=chunkconsistency(seqforchunks,'e....'); %edhDB
bird4{3,1}=chunkconsistency(seqforchunks,'Z.'); %gg
unq=unique(seqforchunks);
allletters=[char(97:122),char(65:90),'0123456789']; %in case some extra chars are needed, i added char 60:64
avletters1=allletters(~ismember(allletters,unq));
Ysym=avletters1(1);
seqforchunks2=regexprep(seqforchunks,'Yi',Ysym);
bird4{4,1}=chunkconsistency(seqforchunks2,[Ysym,'.']);%Yig5


[chunks2,~,seqforchunks,chunks2replace,labelidx2]=seq_chunkextractionfunc(b4_post_R,0);
close all

bird4{1,2}=chunkconsistency(seqforchunks,'F.'); % Fc
bird4{2,2}=chunkconsistency(seqforchunks,'e....'); %edhDB
bird4{3,2}=chunkconsistency(seqforchunks,'W.'); %gg this chunk seems to disappear 
unq=unique(seqforchunks);
allletters=[char(97:122),char(65:90),'0123456789']; %in case some extra chars are needed, i added char 60:64
avletters1=allletters(~ismember(allletters,unq));
Ysym=avletters1(1);
seqforchunks2=regexprep(seqforchunks,'Yi',Ysym);
bird4{4,2}=chunkconsistency(seqforchunks2,[Ysym,'.']);
%% bird5

[chunks2,~,seqforchunks,chunks2replace,labelidx]=seq_chunkextractionfunc(b5_pre_R,0);
close all

bird5{1,1}=chunkconsistency(seqforchunks,'C.'); %CF
bird5{2,1}=chunkconsistency(seqforchunks,'l.'); %lb
bird5{3,1}=chunkconsistency(seqforchunks,'V..'); %gac

unq=unique(seqforchunks);
allletters=[char(97:122),char(65:90),'0123456789']; %in case some extra chars are needed, i added char 60:64
avletters1=allletters(~ismember(allletters,unq));
Ysym=avletters1(1);
seqforchunks2=regexprep(seqforchunks,'Yi',Ysym);
bird5{4,1}=chunkconsistency(seqforchunks2,[Ysym,'.']);%Yig6

[chunks2,~,seqforchunks,chunks2replace,labelidx2]=seq_chunkextractionfunc(b5_post_R,0);
close all

bird5{1,2}=chunkconsistency(seqforchunks,'C.'); %CF
bird5{2,2}=chunkconsistency(seqforchunks,'[lQ].',2); %lb this chunk disappears after lesion!, l is split into 2 based on previous seq
bird5{3,2}=chunkconsistency(seqforchunks,'W..'); %g2ac2

unq=unique(seqforchunks);
allletters=[char(97:122),char(65:90),'0123456789']; %in case some extra chars are needed, i added char 60:64
avletters1=allletters(~ismember(allletters,unq));
Ysym=avletters1(1);
seqforchunks2=regexprep(seqforchunks,'Yi',Ysym);
bird5{4,2}=chunkconsistency(seqforchunks2,[Ysym,'.']);%Yig6

%% bird6

[chunks2,~,seqforchunks,chunks2replace,labelidx]=seq_chunkextractionfunc(b6_pre_R,0);
close all

bird6{1,1}=chunkconsistency(seqforchunks,'k..'); %kcl
bird6{2,1}=chunkconsistency(seqforchunks,'4.'); %cB

unq=unique(seqforchunks);
allletters=[char(97:122),char(65:90),'0123456789']; %in case some extra chars are needed, i added char 60:64
avletters1=allletters(~ismember(allletters,unq));
Ysym=avletters1(1);
seqforchunks2=regexprep(seqforchunks,'Yi',Ysym);
bird6{3,1}=chunkconsistency(seqforchunks2,[Ysym,'.']);%Yim2

[chunks2,~,seqforchunks,chunks2replace,labelidx2]=seq_chunkextractionfunc(b6_post_R,0);
close all

bird6{1,2}=chunkconsistency(seqforchunks,'k..'); %kcl
bird6{2,2}=chunkconsistency(seqforchunks,'6.'); %cB

unq=unique(seqforchunks);
allletters=[char(97:122),char(65:90),'0123456789']; %in case some extra chars are needed, i added char 60:64
avletters1=allletters(~ismember(allletters,unq));
Ysym=avletters1(1);
seqforchunks2=regexprep(seqforchunks,'Yi',Ysym);
bird6{3,2}=chunkconsistency(seqforchunks2,[Ysym,'.']);%Yim2
%% bird7

[chunks2,~,seqforchunks,chunks2replace,labelidx]=seq_chunkextractionfunc(b7_pre_R,0);
close all

bird7{1,1}=chunkconsistency(seqforchunks,'F.'); %Fg1
bird7{2,1}=chunkconsistency(seqforchunks,'l..'); %lBd1
bird7{3,1}=chunkconsistency(seqforchunks,'ax..',4); %axyd

unq=unique(seqforchunks);
allletters=[char(97:122),char(65:90),'0123456789']; %in case some extra chars are needed, i added char 60:64
avletters1=allletters(~ismember(allletters,unq));
Ysym=avletters1(1);
seqforchunks2=regexprep(seqforchunks,'Yi',Ysym);
bird7{4,1}=chunkconsistency(seqforchunks2,[Ysym,'.']);%Yil
[chunks2,~,seqforchunks,chunks2replace,labelidx2]=seq_chunkextractionfunc(b7_post_R,0);
close all

bird7{1,2}=chunkconsistency(seqforchunks,'F.'); %Fg1
bird7{2,2}=chunkconsistency(seqforchunks,'l..'); %lBd1
bird7{3,2}=chunkconsistency(seqforchunks,'ax..',4); %axyd

unq=unique(seqforchunks);
allletters=[char(97:122),char(65:90),'0123456789']; %in case some extra chars are needed, i added char 60:64
avletters1=allletters(~ismember(allletters,unq));
Ysym=avletters1(1);
seqforchunks2=regexprep(seqforchunks,'Yi',Ysym);
bird7{4,2}=chunkconsistency(seqforchunks2,[Ysym,'.']);%Yil
%% all consistency
allcons=cat(1,bird1,bird2,bird3,bird4,bird5, bird6,bird7);

%% significance
allconsmat=cell2mat(allcons);
[p,h,stats] = signrank(allconsmat(:,1),allconsmat(:,2),'alpha',0.05);
%% consistency figure
f=figure
hold on
for m=1:size(bird1,1)
    plot([1,2], [bird1{m,1},bird1{m,2}],'Marker','o','MarkerSize',12,'MarkerEdgeColor','k','MarkerSize',8,'Color',[0.5 0.5 0.5],'LineStyle','-', 'LineWidth', 0.5)
    hold on
end
for m2=1:size(bird2,1)
    plot([1,2], [bird2{m2,1},bird2{m2,2}],'Marker','+','MarkerSize',12,'MarkerEdgeColor','k','Color','k','LineStyle','-', 'LineWidth', 0.5)
    hold on
end
for m3=1:size(bird3,1)
    plot([1,2], [bird3{m3,1},bird3{m3,2}],'Marker','x','MarkerSize',12,'MarkerEdgeColor','k','Color',[0.5 0.5 0.5],'LineStyle','-', 'LineWidth', 0.5)
    hold on
end
for m4=1:size(bird4,1)
    plot([1,2], [bird4{m4,1},bird4{m4,2}],'Marker','s','MarkerSize',12,'MarkerEdgeColor','k','Color',[0.5 0.5 0.5],'LineStyle','-', 'LineWidth', 0.5)
    hold on
end
for m5=1:size(bird5,1)
    plot([1,2], [bird5{m5,1},bird5{m5,2}],'Marker','d','MarkerSize',12,'MarkerEdgeColor','k','Color',[0.5 0.5 0.5],'LineStyle','-', 'LineWidth', 0.5)
    hold on
end
for m6=1:size(bird6,1)
    plot([1,2], [bird6{m6,1},bird6{m6,2}],'Marker','^','MarkerSize',12,'MarkerEdgeColor','k','Color',[0.5 0.5 0.5],'LineStyle','-', 'LineWidth', 0.5)
    hold on
end
for m7=1:size(bird7,1)
    plot([1,2], [bird7{m7,1},bird7{m7,2}],'Marker','s','MarkerSize',12,'MarkerEdgeColor','k','Color',[0.5 0.5 0.5],'LineStyle','-', 'LineWidth', 0.5)
    hold on
end
%plot([1,2],[median(allconsmat(:,1)),median(allconsmat(:,2))],'LineStyle','-', 'LineWidth', 3)
boxplot(cell2mat(allcons),'Whisker',1,'Colors','br','Symbol','.k')
xticklabels({'Prelesion','Postlesion'});
% ylim([0,2])
% yticklabels([0:0.2:2])
ylabel('Consistency','FontSize',14);
title('Chunks consistency','(n=23 chunks from 7 birds)','FontSize',14);
box off
birdchunks=[bird1;bird2;bird3;bird4;bird5;bird6;bird7];
numchunks=length(birdchunks);

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',14,'FontWeight','bold')
set(gca,'XTickLabelMode','auto')
xticklabels({'Prelesion','Postlesion'});

xlim([0.75 2.25])
daspect([6 1 1])