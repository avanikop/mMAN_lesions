% MODIFIED 14.07.2023
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
%% overall transition entropy
transent_all=[bird1_rep_ent;bird2_rep_ent;bird3_rep_ent;bird4_rep_ent;bird5_rep_ent;bird6_rep_ent;bird7_rep_ent];
%transent_all_nai=[bird1_nai_ent;bird2_nai_ent;bird3_nai_ent;bird4_nai_ent;bird5_nai_ent;bird6_nai_ent;bird7_nai_ent];
f1=figure('Name','Total transition entropy (n=7)')
boxplot(transent_all,'Whisker',1,'Colors','br')
hold on
Markers={'o','+','x','s','d','^','h'};


for m=1:length(transent_all)
    plot(1,transent_all(m,1),Markers{m},'MarkerSize',12,'MarkerEdgeColor','b');
    plot(2,transent_all(m,2),Markers{m},'MarkerSize',12,'MarkerEdgeColor','r');
    plot([1,2], [transent_all(m,1),transent_all(m,2)],'-','Color',[0.5 0.5 0.5], 'LineWidth', 0.5)
end
xticklabels({'Prelesion','Postlesion'})

ylabel('Transition Entropy')
title('Total transition entropy (n=7)')
hold off
box off

daspect([3 1 1])
% transition entropy wilcoxon sign rank test
[p,h,stats] = signrank(transent_all(:,1),transent_all(:,2),'alpha',0.05)

