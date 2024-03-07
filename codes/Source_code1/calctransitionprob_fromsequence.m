function[testcell,divprobtoplot,TPnew,freq]=calctransitionprob_fromsequence(labels,seq,eliminaterare)
% testcell= syl positions for divprob
% divprobtoplot = transition matrix
% TPnew = rawmtx/matrix of transition counts
% freq = frequency of the chunks
for i=1:length(labels)
    for j=1:length(labels)
        testcell{i,j}=[labels{i}, labels{j}];
        countcell{i,j}=length(regexp(seq,testcell{i,j}));
    end
end
TPnew=cell2mat(countcell);
% according to Lena's logic I should be only summing over rows
% tot=sum(nonzeros(TPnew));
tot=sum((TPnew),2);
sumcols=sum(TPnew,1);
grandtot=sum(tot);
freq=sumcols/grandtot;

if eliminaterare
%eliminate rare notes
    rarenodes = tot./sum(tot)<0.01; %0.025;
    rarenodes=rarenodes(rarenodes~=1); %so that the start is not removed
    tot(rarenodes)=0;
    TPnew(rarenodes,:)=0;
    TPnew(:,rarenodes)=0;
end

divprob=(TPnew)./repmat(tot,1,length(TPnew));
divprob = round(divprob*100);
divprobtoplot = divprob;
divprobtoplot(isnan(divprobtoplot))=0;