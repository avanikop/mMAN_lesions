function[chunks2,percY,seqforchunks,chunks2replace,labelidx,newseq,divprobtoplot2,patterncell2,labels2,numnewsyls]=seq_chunkextractionfunc(seqnew,~)
% This function tried to automate extraction of chunks based on chi sq
% analysis. First, we try to determine if a syllable's next transition
% depends on the syllable that comes before it using chi sq analysis. If it
% does, we relable the syllable and use it as a different 'state' of the
% same syllable. Then using my chunk exrtraction function, we create chunks
% with one-in-one-out branches, >80% transition prob as middle nodes.
% Bonferroni correcion added where alpha = 0.01/n where n = number of
% comparisons made 08.03.2023
% Detailed notes are in the description of the function. 
% This function gives you chunks and plots old transition diagram and new
% transition diagram automatically.
% Input= seqnew, a clean sequence string with intro notes replaced and
% start symbol 'Y' already present
% if you dont want a plot write 0 for plotting
% you need seqforchunks and chunks2replace for chunk consistency analysis
% and labelidx also
% newseq gives you the seq w/ chunks replaced 
% divprobtoplot2 = final transitionprob matrix
% and patterncell2 gives you the patterncell for the diagraph
% labels2 = what are the replaced labels
% numnewsyls = number of additional states
%%
unq=unique(seqnew);
% for all syllables
% probabilities of all syls to everything
for i=1:length(unq)
    for j=1:length(unq)
        transmat{i}{j}=[unq(i),unq(j)];
        countst{i}(j)=length(strfind(seqnew,[unq(i),unq(j)])); % counting WITH overlab
    end
    probst{i}=countst{i}/sum(countst{i});
end
%%
% prob of all syls to everything given each syl before
% make matrix of all by all
for xi=1:length(unq)
    for xj=1:length(unq)
        mat{xi,xj}=[unq(xi),unq(xj)];
    end
end
%%
for yi=1:length(unq)
    for yj=1:length(unq)
        for yk=1:length(unq)
            test{yi,yj}{yk}=[mat{yi,yj},unq(yk)]; %every column is x-yoursyl-y
            countstg{yi,yj}(yk)=length(strfind(seqnew,[mat{yi,yj},unq(yk)]));
        end
    end
end
%% chisq test
h = nan(length(unq));
chi = nan(length(unq));
p = nan(length(unq));
countst2=countst; %saving countst before I modify it to delete linear branches
% remove increadibly small branches which are likely mislabels
for ix=1:length(countst)
    num=countst{ix}./sum(countst{ix});
    pink=countst{ix};
    pink(num<0.01)=0;
    countst{ix}=pink;
end
nzeros=cellfun(@(x) length(nonzeros(x))>1,countst, 'UniformOutput',0);
nz2=cellfun(@(x) x==1,nzeros);
countst(~nz2)={0};
% remove all comparisons of column 1
colY=strfind(unq,'Y'); % this is the column of countst which corresponds to Y
countst{colY}=0; %because Y doesnt actually depend on anything
% filtering out branches of countstg which are below 0.1% of the total
% transitions
testcountstg=countstg; % 'saving' countstg
%testcountstg is the one that stays intact; countstg will be modified;
%go through all cells and delete whole cells if sum within the cell is less
%than 1% of the total times this syllable is observed
sumcols=sum(cellfun(@(x) sum(x),testcountstg),1);
for i=1:size(countstg,1)
    for j=1:size(countstg,2)
        if sum(testcountstg{i,j})/sumcols(j)>0.01 %filtering out cells of countstg which appear very rarely %10%
            newtestcnt{i,j}=testcountstg{i,j};
        else newtestcnt{i,j}=[];
        end
    end
end
countstg=newtestcnt;

% carrying out chi sq test:
for zi=1:length(unq)
    existingpostsyls = find(countst{zi}./sum(countst{zi})>0.01); %vector which gives indices of not-rare transitions; >0.01
    for zj=1:length(unq)
            %numstates should be number of possibilities of
            %preceeding syl
            numstates=sum(cellfun(@(x) ~isempty(x),countstg(:,zj)));
        if ~isempty(countstg{zj,zi})
            if any(countstg{zj,zi}(existingpostsyls)./sum(countst{zi}(existingpostsyls))>0.01)
                if length(nonzeros(countstg{zj,zi}(existingpostsyls)))>1 %because you need at least a vector of 2 for chi sq test
                    %if this conditional branch is not 1% f the total branches of
                    %that syllable 
                    testnew{zj,zi}=test{zj,zi}(existingpostsyls);
                    [h(zj,zi),chi(zj,zi),p(zj,zi)] = chisq_2dist(countst{zi}(existingpostsyls),countstg{zj,zi}(existingpostsyls),0.01/numstates);
                    %changed on 08.03.23, bonferroni correction
                    %if zi =1, first element of countst ie YY,Ya,Yb etc filtered to nonzero
                    %transitions vs all rows of countstg over first column (ie
                    %YYY,YYa,YYb etc) filtered to nonzero transitions
                elseif length(nonzeros(countstg{zj,zi}(existingpostsyls))) == 1
                    h(zj,zi) = 1; %automatically, it becomes a separate state in itself
                end
            end
        end
    end
end

%% keep only the ones with h=1
[a,b]=find(h==1);
cellremaining=cell(length(unq));
for i=1:length(a)
    cellremaining{a(i),b(i)}=mat{a(i),b(i)};
end
%% now do the relabelling thing
% find letters not present in unq
allletters=[char(97:122),char(65:90),'0123456789']; %in case some extra chars are needed, i added char 60:64
avletters1=allletters(~ismember(allletters,unq)); %available letters
avletters=avletters1(end:-1:1); %just to make them more noticable when replaced

%% collecting indices for replacing and replacing together
states=cellfun(@(x) ~isempty(x),cellremaining);
sumstates=sum(states,1);
torelabel=cellremaining(:,sumstates>=1);

tt=1; %idc for avletters
for ti=1:size(torelabel,2)
    stateidx=find(~(cellfun('isempty',torelabel(:,ti)))); %index of combos to relabel in the ti row
    for tj=1:length(stateidx)
        repseq=torelabel{stateidx(tj),ti};
        labelidx{tt,1}=strfind(seqnew,repseq)+1; %because you're looking to replace the 2nd syl
        labelidx{tt,2}=[avletters(tt)];
        labelidx{tt,3}={[repseq(2),num2str(tj)]};
        labelidx{tt,4}={repseq};
        tt=tt+1;
    end
end
%% see if the replaced syls actually are different WITHIN themselves
[a,b]=find(h==1);
remcountstg=cell(size(countstg));
for i=1:length(a)
    remcountstg{a(i),b(i)}=countstg{a(i),b(i)}; %remcountstg = counstg where h==1 ie the dists are different
end
states2=cellfun(@(x) ~isempty(x),remcountstg);
sumstates=sum(states2,1);
remcountstg=remcountstg(:,sumstates>1);
which_to_merge=cell(size(remcountstg));
names_to_merge=cell(size(remcountstg));
for ip=1:size(remcountstg,2)
    emptidx=~(cellfun('isempty',remcountstg(:,ip)));
    x=remcountstg(emptidx,ip);
    t=cellfun(@(y) y./sum(y),x,'UniformOutput',false); %actual percentages are stored here for a
    a=cellfun(@(y) find((y./sum(y))>0.05),x,'UniformOutput',false); %branches below 5% are ignored
    f=cellfun(@(y) find((y./sum(y))>0.05),remcountstg(:,ip),'UniformOutput',false); 
    a = a(:);
    [~,b,c] = unique(cellfun(@char,a,'un',0));
    lo = histc(c,1:max(c));
    loo = lo(:) > 1;
    out = [a(b(loo))];
    % condition for if length(out{ir}) = 1 has been added
    for ir=1:length(out)
        index=[];
        for k = 1:numel(f)
            try
                if f{k} == cell2mat(out(ir))
                    ind = (f{k} == cell2mat(out(ir)));
                    index=[index;k]; % found indices of remcountstg of cells with this particular output 
                end
            end
        end
        if length(out{ir})>1 %then do all the testing
            % test these against each other
            % what if there are more than 2 options? unlikely, but still..
            pairwise_combinations=nchoosek(index,2); %here you get combinations of indices
            numstates2=length(pairwise_combinations)/2; %modified on 31.03.2023
            for id=1:size(pairwise_combinations,1) %along rows
                data1=remcountstg{pairwise_combinations(id,1),ip};
                data2=remcountstg{pairwise_combinations(id,2),ip};
                h_testing_in_pairs=chisq_2dist(data1,data2,0.01/numstates2); %modified on 08.03.2023
                if h_testing_in_pairs == 0
                    which_to_merge{ir,ip}=[which_to_merge{ir,ip};pairwise_combinations(id,1),ip;pairwise_combinations(id,2),ip];
                end
            end
        elseif length(out{ir})==1 %merge them automatically if there is only one major branch
            for ig=1:length(index)
                which_to_merge{ir,ip}=[which_to_merge{ir,ip};index(ig),ip]; %make vector of all indices + col number
            end
        end
        which_to_merge{ir,ip}=unique(which_to_merge{ir,ip},'rows');
        try %in case which_to_merge_is empty
        rows=which_to_merge{ir,ip}(:,1);
        cols=which_to_merge{ir,ip}(:,2);
        names_to_merge{ir,ip}=[];
            for ix=1:length(rows)
                names_to_merge{ir,ip}=[string(names_to_merge{ir,ip});string(torelabel{rows(ix),cols(ix)})];
            end
        end
    end
end
% node for which_to_merge: ip= which COLUMN of
% torelabel to merge,go over ALL ir's; those are the
% different TYPES of combos available. Everything with the indices that sit
% inside which_to_merge{ir,ip} should be the SAME syllable 
%% now prune labelidx based on your which_to_merge matrix
% go over every element of names_to_merge
for num_merg=1:numel(names_to_merge)
    thiscell=names_to_merge{num_merg};
    repidx=[];
    for idx=1:length(thiscell) %because 2 or more syllables might be alike and in need of merge
        str=thiscell(idx);
        for collen=1:size(labelidx,1)
            if strfind(labelidx{collen,4}{:},str)
                repidx=[repidx;collen];
            end
        end
    end
    if ~isempty(repidx)
        torepwith=labelidx{repidx(1),3};
        torepwith2=labelidx{repidx(1),2};
        for replen=2:length(repidx)
            labelidx{repidx(replen),3}=torepwith;
            labelidx{repidx(replen),2}=torepwith2;
        end
    end
end
%% replacing stuff
for xi=1:size(labelidx,1) %along rows
    seqnew(labelidx{xi,1})=labelidx{xi,2};
end
% this is the seqnew I'd want for my chunks consistency analysis
seqforchunks=seqnew;
%%
[a,freq]=uniquestring(seqnew,0.9); % was 0.9 HERE changed from 0.05 to 1 for bubu but bk2bk10 needs <1
labels=cell(1,length(a));
%unq=unique(seqnew);
for unqidx=1:length(a)
    labels{unqidx}=a(unqidx);
end
%% create patterncell
patterncell=labels;
for ix=1:size(labelidx,1)
    idx=strcmp(labels,labelidx{ix,2});
    patterncell(idx)=labelidx{ix,3};
end
%%
[~,divprobtoplot,~]=calctransitionprob_fromsequence(labels,seqnew,1);
% delete tiny percentage branches from divprobtoplot before making the
% graph
%divprobtoplot(abs(divprobtoplot)<4)=0; % 5% and above branches only

g=seq_plot_digraph(divprobtoplot,patterncell,'TGraph','PRE',freq*500);
for im=1:length(patterncell)
    labelsnum{im}=num2str(im);
end
%g = digraph(divprobtoplot,labelsnum)
if ~exist('plotting','var')
f=seq_plot_digraph(divprobtoplot,labelsnum,'TGraph','PRE',freq*500);
%%uncomment this only if you want to troubleshoot the graph and want to
%%know the node numbers quickly

seq_plot_transitionmatrix(divprobtoplot,patterncell,1,'TMat','PRE'); 
end
%% chunk extraction
try
[paths]=chunkextraction(g,80); 
chunks = cellfun(@(x) [patterncell(x)],paths,'UniformOutput',false);
chunks2=cellfun(@(x) [x{:}],chunks,'UniformOutput',false); %chunk names to plot
[~,uidx] = unique(chunks2,'stable');
chunks2=chunks2(uidx);
chunks2plot=cellfun(@(x) [labels(x)],paths,'UniformOutput',false);
chunks2replace=cellfun(@(x) [x{:}],chunks2plot,'UniformOutput',false); %chunk names to replace
chunks2replace=chunks2replace(uidx);
%% create transition matrix
[newseq,pp,newchunknames] = ak_replacechunks(seqnew,chunks2replace);

%% create labelcell2
[a,freq]=uniquestring(newseq,0);
% adding a condition such that the Yi chunk is never removed
chunk_with_Y=cellfun(@(x) strfind(x,'Y'),pp,'UniformOutput',false);
num_Y=find(~cellfun(@isempty,chunk_with_Y));
Y_in_a=strfind(a,newchunknames(num_Y));
freqY=freq(Y_in_a);
percY=freqY*100; %in percentage

%here the threshold is changed such that Yi HAS to be in the graph.

[a,freq]=uniquestring(newseq,percY*0.70); %changed by Avani and Lena on 31.03.2023
labels2=cell(1,length(a));
for unqidx=1:length(a)
    labels2{unqidx}=a(unqidx);
end
[~,divprobtoplot2,~]=calctransitionprob_fromsequence(labels2,newseq,0);
%% create patterncell2 to plot in final transition graph
patterncell2=labels2;
for ix=1:size(labelidx,1)
    idx=strcmp(labels2,labelidx{ix,2});
    patterncell2(idx)=labelidx{ix,3};
end
% replace chunks to plot in final transition graph
for iy=1:length(newchunknames)
    %Here i use try in case one of the chunks is removed during freqency thresholding in line 250
    try
    idx2=strcmp(patterncell2,newchunknames(iy)); %for example, found pos of A
    patterncell2{idx2}=chunks2{iy};
    end
end
if ~exist('plotting','var')
% plotting
    seq_plot_digraph(divprobtoplot2,patterncell2,'TGraph','PRE',freq*100);
    seq_plot_transitionmatrix(divprobtoplot2,patterncell2,1,'TMat','PRE');
end
% i need an output seq which has all the 'new' states replaced
%% how many additional states
numnewsyls=length(unique(seqforchunks))-length(unq);
catch
    '!!!!SOME MASSIVE ERROR!!!!'
    chunks2replace=[];
    newseq=0;
    divprobtoplot2=0;
    patterncell2=[];
    labels2=[];
    numnewsyls=0;
    chunks2=0;
    percY=0;
end
end