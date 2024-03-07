% TE vs baseline history dependence
% DATASET OF BIRDS MADE OF SAME SIZE
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
%% 
% to calculate transition entropy per branchpoint, we look at labels2 with
% chunks replaced. Then we look at the identity of the endpoints of these
% chunks and the identity of remaining syls based on contet (ie syllable
% that preceeds it). We find transition entropy at exactly these points and
% then find these points in postlesion to do the same.
%% BIRD1
[~,~,prelesionchunkseq,chunks2replacepre,labelidxpre,~,~,patterncellpre,labelspre]=seq_chunkextractionfunc(b1_pre_R,0);
close all

[~,~,postlesionchunkseq,chunks2replacepost,labelidxpost,~,~,~,labelspost]=seq_chunkextractionfunc(b1_post_R,0);
close all

bpprelesion={'cc','ir','wr','cr','xy','l','r','w'};
bppostlesion={'cc','ir','wr','cr','xy','l','r','w'};

[tebp_bird1,b1histdep]=transent_prevspost(b1_pre_R,bpprelesion,b1_post_R,bppostlesion); % n*2 matrix of transition entropy per branchpoint
[bird1_nai_ent,bird1_rep_ent]=calc_overallentropy(b1_pre_R,b1_post_R,prelesionchunkseq,postlesionchunkseq); % 2 1*2 vectors of 
% overall transition entropy before and after replacing states w/ chi sq analysis
%calculate duration of gap preceeding this point

%% BIRD2

[~,~,prelesionchunkseq,chunks2replacepre,labelidxpre,~,~,patterncellpre,labelspre]=seq_chunkextractionfunc(b2_pre_R,0);
close all

[~,~,postlesionchunkseq,chunks2replacepost,labelidxpost,~,~,~,labelspost]=seq_chunkextractionfunc(b2_post_R,0);
close all

bpprelesion={'gd','fk','i','[kh]h','dh'};
bppostlesion={'gd','fk','i','[kh]h','dh'};
[tebp_bird2,b2histdep]=transent_prevspost(b2_pre_R,bpprelesion,b2_post_R,bppostlesion,[2,2,1,2,2]); % n*2 matrix of transition entropy per branchpoint
[bird2_nai_ent,bird2_rep_ent]=calc_overallentropy(b2_pre_R,b2_post_R,prelesionchunkseq,postlesionchunkseq); % 2 1*2 vectors of 
% overall transition entropy before and after replacing states w/ chi sq analysis
%% bird3

[~,~,prelesionchunkseq,chunks2replacepre,labelidxpre,~,~,patterncellpre,labelspre]=seq_chunkextractionfunc(b3_pre_R,0);
close all

[~,~,postlesionchunkseq,chunks2replacepost,labelidxpost,~,~,~,labelspost]=seq_chunkextractionfunc(b3_post_R,0);
close all

bpprelesion={'B','he','ie','g','s','Be','Bf','ef'};
bppostlesion={'B','he','ie','g','s','Be','Bf','ef'};
[tebp_bird3,b3histdep]=transent_prevspost(b3_pre_R,bpprelesion,b3_post_R,bppostlesion); % n*2 matrix of transition entropy per branchpoint
[bird3_nai_ent,bird3_rep_ent]=calc_overallentropy(b3_pre_R,b3_post_R,prelesionchunkseq,postlesionchunkseq); % 2 1*2 vectors of 
% overall transition entropy before and after replacing states w/ chi sq analysis
%% bird4

[~,~,prelesionchunkseq,chunks2replacepre,labelidxpre,~,~,patterncellpre,labelspre]=seq_chunkextractionfunc(b4_pre_R,0);
close all

[~,~,postlesionchunkseq,chunks2replacepost,labelidxpost,~,~,~,labelspost]=seq_chunkextractionfunc(b4_post_R,0);
close all

bpprelesion={'gc','Fc','ig','gg','B','cg','a'};
bppostlesion={'gc','Fc','ig','gg','B','cg','a'};

[tebp_bird4,b4histdep]=transent_prevspost(b4_pre_R,bpprelesion,b4_post_R,bppostlesion); % n*2 matrix of transition entropy per branchpoint
[bird4_nai_ent,bird4_rep_ent]=calc_overallentropy(b4_pre_R,b4_post_R,prelesionchunkseq,postlesionchunkseq); % 2 1*2 vectors of 
% overall transition entropy before and after replacing states w/ chi sq analysis
%% bird5

[~,~,prelesionchunkseq,chunks2replacepre,labelidxpre,~,~,patterncellpre,labelspre]=seq_chunkextractionfunc(b5_pre_R,0);
close all


[~,~,postlesionchunkseq,chunks2replacepost,labelidxpost,~,~,~,labelspost]=seq_chunkextractionfunc(b5_post_R,0);
close all

bpprelesion={'F','ig','b','ac','gg','dg','cg','ff','df','d'};
bppostlesion={'F','ig','b','ac','gg','dg','cg','ff','df','d'};
[tebp_bird5,b5histdep]=transent_prevspost(b5_pre_R,bpprelesion,b5_post_R,bppostlesion); % n*2 matrix of transition entropy per branchpoint
[bird5_nai_ent,bird5_rep_ent]=calc_overallentropy(b5_pre_R,b5_post_R,prelesionchunkseq,postlesionchunkseq); % 2 1*2 vectors of 
% overall transition entropy before and after replacing states w/ chi sq analysis
%% bird6

[~,~,prelesionchunkseq,chunks2replacepre,labelidxpre,~,~,patterncellpre,labelspre]=seq_chunkextractionfunc(b6_pre_R,0);
close all

[~,~,postlesionchunkseq,chunks2replacepost,labelidxpost,~,~,~,labelspost]=seq_chunkextractionfunc(b6_post_R,0);
close all

bpprelesion={'Bj','cB','im','kcl','jm','gm','mj','lj','b','g','l','j'};
bppostlesion={'Bj','cB','im','kcl','jm','gm','mj','lj','b','g','l','j'};
[tebp_bird6,b6histdep]=transent_prevspost(b6_pre_R,bpprelesion,b6_post_R,bppostlesion); % n*2 matrix of transition entropy per branchpoint
[bird6_nai_ent,bird6_rep_ent]=calc_overallentropy(b6_pre_R,b6_post_R,prelesionchunkseq,postlesionchunkseq); % 2 1*2 vectors of 
% overall transition entropy before and after replacing states w/ chi sq analysis

%% bird7

[~,~,prelesionchunkseq,chunks2replacepre,labelidxpre,~,~,patterncellpre,labelspre]=seq_chunkextractionfunc(b7_pre_R,0);
close all

[~,~,postlesionchunkseq,chunks2replacepost,labelidxpost,~,~,~,labelspost]=seq_chunkextractionfunc(b7_post_R,0);
close all

bpprelesion={'Fg','C','ig','Bd','xyd','e'};
bppostlesion={'Fg','C','ig','Bd','xyd','e'};

[tebp_bird7,b7histdep]=transent_prevspost(b7_pre_R,bpprelesion,b7_post_R,bppostlesion,[2,1,2,2,3,1]); % n*2 matrix of transition entropy per branchpoint
[bird7_nai_ent,bird7_rep_ent]=calc_overallentropy(b7_pre_R,b7_post_R,prelesionchunkseq,postlesionchunkseq); % 2 1*2 vectors of 
% overall transition entropy before and after replacing states w/ chi sq analysis
%% plotting
%%  history dependence vs change in TEBP
tebp=[tebp_bird1;tebp_bird2;tebp_bird3;tebp_bird4; tebp_bird5; tebp_bird6; tebp_bird7];
tebpdiff_all=tebp(:,2)-tebp(:,1);

histdeppre_all=[b1histdep(:,1);b2histdep(:,1);b3histdep(:,1);b4histdep(:,1);...
    b5histdep(:,1);b6histdep(:,1);b7histdep(:,1)];
tebp_bird1_diff=-[tebp_bird1(:,1)-tebp_bird1(:,2)];
tebp_bird2_diff=-[tebp_bird2(:,1)-tebp_bird2(:,2)];
tebp_bird3_diff=-[tebp_bird3(:,1)-tebp_bird3(:,2)];
tebp_bird4_diff=-[tebp_bird4(:,1)-tebp_bird4(:,2)];
tebp_bird5_diff=-[tebp_bird5(:,1)-tebp_bird5(:,2)];
tebp_bird6_diff=-[tebp_bird6(:,1)-tebp_bird6(:,2)];
tebp_bird7_diff=-[tebp_bird7(:,1)-tebp_bird7(:,2)];

[correlationhistdep,pvalhistdep]=corr(tebpdiff_all,histdeppre_all);
% calculating and plotting regression line
x=histdeppre_all(:,1);
X=[ones(length(x),1),x];
y=tebpdiff_all;
b=X\y;
yCalc1=X*b;
Rsq=1-sum((y-yCalc1).^2)/sum((y - mean(y)).^2);

f5=figure('Name','SIGNS_History dependence vs change in TEBP')

scatter(b1histdep(:,1),tebp_bird1_diff,'o',...
              'LineWidth',1.5)
hold on
scatter(b2histdep(:,1),tebp_bird2_diff,'+',...
              'LineWidth',1.5)

scatter(b3histdep(:,1),tebp_bird3_diff,'x',...
              'LineWidth',1.5)

scatter(b4histdep(:,1),tebp_bird4_diff,'s',...
              'LineWidth',1.5)

scatter(b5histdep(:,1),tebp_bird5_diff,'d',...
              'LineWidth',1.5)

scatter(b6histdep(:,1),tebp_bird6_diff,'^',...
              'LineWidth',1.5)

scatter(b7histdep(:,1),tebp_bird7_diff,'h',...
              'LineWidth',1.5)
xlim([0 1.1])
h=refline(0,0)
h.Color='k'
h.LineStyle='--'
%legend('bird1','bird2','bird3','bird4','bird5','bird6','bird7')
title('Transition entropy change vs baseline history dependence','(n=56 branchpoints in n=7 birds)')
xlabel('History dependence')
ylabel('Change in entropy')
txt1={['Corr.coef.=',num2str(round(correlationhistdep,3))],['p=',num2str(round(pvalhistdep,3))]};
text(0.75,0.8,txt1);

hold off