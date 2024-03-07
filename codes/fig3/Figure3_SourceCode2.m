%% chunk fig 3b
%% bird2

fileID = fopen('D:\analysis\data_for_elife_mMAN\Source_data_1\bird2_prelesion','r');
b2_pre = fscanf(fileID, '%s')
fclose(fileID);
fileID = fopen('D:\analysis\data_for_elife_mMAN\Source_data_1\bird2_postlesion','r');
b2_post = fscanf(fileID, '%s')
fclose(fileID);
%%
point=(strfind(b2_pre,'c')); %strfind > regexp
for i=1:length(point)
    if b2_pre(point(i)-1)=='c'
        point(i)=0; %so that only first c's remain
    end
end
point4=point(point~=0); %removing all second c's
point=point4;
% point=point-1; %alinging to the old regexp code where numbering starts at .c...
% point=(regexp(b2_pre,'c...'));
for i=1:5
    for j=1:size(point,2)
        point(i+1,j)=point(1,j)+i;
    end
end
listprelesion=b2_pre(point');

%%
point2=(strfind(b2_post,'c')); %strfind > regexp
for i=1:length(point2)
    if b2_post(point2(i)-1)=='c'
        point2(i)=0; %so that only first c's remain
    end
end
point3=point2(point2~=0); %removing all second c's
point2=point3;
%point2=point2-1; %alinging to the old regexp code where numbering starts at .c...
for i=1:5
    for j=1:size(point2,2)
        point2(i+1,j)=point2(1,j)+i;
    end
end
% since point2 is longer than b2_post
newpoint=point2';
newpoint= newpoint(1:end-1,:);
listpostlesion=b2_post(newpoint);
%% figure
% create number matrix with listprelesion
prelesionmatrix=zeros(size(listprelesion));
% if there's Y, eveyrthing after it has value NaN
% dccll has low values eg d=11, c=10,
% l=12,f=14,k=16,A=18,h=18,i=18,j=20,'g'=22,t=24
for i=1:size(listprelesion,1)
    for j= 1:size(listprelesion,2)
        n=listprelesion(i,j);
    switch n
        case 'b'
            postlesionmatrix(i,j)=11;
        case 'd'
            prelesionmatrix(i,j)=10;
        case 'c'
            prelesionmatrix(i,j)=9;
        case 'l'
            prelesionmatrix(i,j)=8;
        case 'f'
            prelesionmatrix(i,j)=7;
        case 'k'
            prelesionmatrix(i,j)=6;
        case 'A'
            prelesionmatrix(i,j)=5;
        case 'h'
            prelesionmatrix(i,j)=5;
        case 'i'
            prelesionmatrix(i,j)=5;
        case 'j'
            prelesionmatrix(i,j)=4;
        case 'g'
            prelesionmatrix(i,j)=3;
        case 't'
            prelesionmatrix(i,j)=2;
        case 'Y'
            prelesionmatrix(i,j)=1;
        otherwise
            prelesionmatrix(i,j)=2;
            n
    end
    end
end
% convert everything after y to 1
for i=1:size(prelesionmatrix,1)
    for j= 1:size(prelesionmatrix,2)
        if prelesionmatrix(i,j)==1
            prelesionmatrix(i,j:size(prelesionmatrix,2))=1;
        end
    end
end

%% figure
% create number matrix with listpostlesion
postlesionmatrix=zeros(size(listpostlesion));
% if there's Y, eveyrthing after it has value NaN
% dccll has low values eg d=12, c=11
% l=10,f=9,k=8,A=7,h=7,i=7,j=6,'g'=5,t=4 , Y=0
% highest is 12
for i=1:size(listpostlesion,1)
    for j= 1:size(listpostlesion,2)
        n=listpostlesion(i,j);
    switch n
        case 'b'
            postlesionmatrix(i,j)=11;
        case 'd'
            postlesionmatrix(i,j)=10;
        case 'c'
            postlesionmatrix(i,j)=9;
        case 'l'
            postlesionmatrix(i,j)=8;
        case 'f'
            postlesionmatrix(i,j)=7;
        case 'k'
            postlesionmatrix(i,j)=6;
        case 'A'
            postlesionmatrix(i,j)=5;
        case 'h'
            postlesionmatrix(i,j)=5;
        case 'i'
            postlesionmatrix(i,j)=5;
        case 'j'
            postlesionmatrix(i,j)=4;
        case 'g'
            postlesionmatrix(i,j)=3;
        case 't'
            postlesionmatrix(i,j)=2;
        case 'Y'
            postlesionmatrix(i,j)=1;
    end
    end
end
% convert everything after y to 1
for i=1:size(postlesionmatrix,1)
    for j= 1:size(postlesionmatrix,2)
        if postlesionmatrix(i,j)==1
            postlesionmatrix(i,j:size(postlesionmatrix,2))=1;
        end
    end
end
postlesionmatrix=postlesionmatrix(1:size(prelesionmatrix,1),:);
%% convert to prime
primeprelesion=zeros(size(prelesionmatrix));
for i=1:size(prelesionmatrix,1)
    primeprelesion(i,:)=nthprime(prelesionmatrix(i,:));
end
primepostlesion=zeros(size(postlesionmatrix));
for j=1:size(postlesionmatrix,1)
    primepostlesion(j,:)=nthprime(postlesionmatrix(j,:));
end
 %% figure
% figure('Name','Chunks Prelesion')
% imagesc(prelesionmatrix)
% colormap('jet')
% ax=gca
% xticks(0.5:1:7.5)
% ax.XGrid="on"
% %ax.XTickLabel={' ','1',' ','2',' ','3',' ','4',' ','5',' ','6',' ','7',' '}
% box off
% %%
% figure('Name','Chunks Postlesion')
% imagesc(postlesionmatrix)
% colormap('jet')
% ax=gca
% xticks(0.5:1:7.5)
% ax.XGrid="on"
% %% sankey diagram
% %prelesionsankey=prelesionmatrix/2-3;
% %prelesionsankey(find(prelesionsankey<0))=0;
% figure('Name','bk2bk10 Chunk figure')
% tiledlayout(1,2)
% nexttile
% CreateSankeyPlot(prelesionmatrix)
% title('Prelesion')
% nexttile
% CreateSankeyPlot(postlesionmatrix)
% title('Postlesion')
% %% bar plot
% 
% preunq=unique(primeprelesion,'rows','sorted');
% postunq=unique(primepostlesion,'rows','sorted');
% newprematrix=[];
% for ij=1:length(preunq)
%     num=find(ismember(primeprelesion,preunq(ij,:),'rows'));
%     toadd=primeprelesion(num,:);
%     newprematrix=[newprematrix;toadd];
% end
% newpostmatrix=[];
% for ij=1:length(postunq)
%     num=find(ismember(primepostlesion,postunq(ij,:),'rows'));
%     toadd=primepostlesion(num,:);
%     newpostmatrix=[newpostmatrix;toadd];
% end
% figure()
% imagesc(newprematrix)
% 
% figure()
% imagesc(newpostmatrix)
%% sort according to histogram
primeprelesion = prelesionmatrix;
primepostlesion = postlesionmatrix;

[ai,~,ci]=unique(primeprelesion,'rows','sorted');
chist = histc(ci,unique(ci));
chist = chist./(sum(chist));
outelts = ai(chist>0./100);
outfreq = chist(chist>0./100);
% arranging prelesionmatrix according to frequencies
onebigmatpre=[outfreq(ci),primeprelesion];
sortedbigmatpre=sortrows(onebigmatpre,'ascend');
newpresortedpre=sortedbigmatpre(:,2:end);

[aj,~,cj]=unique(primepostlesion,'rows','sorted');
chist = histc(cj,unique(cj));
chist = chist./(sum(chist));
outelts = aj(chist>0./100);
outfreq = chist(chist>0./100);
% arranging prelesionmatrix according to frequencies
onebigmatpo=[outfreq(cj),primepostlesion];
sortedbigmatpo=sortrows(onebigmatpo,'ascend');
newpresortedpo=sortedbigmatpo(:,2:end);

%% avani makes a colormap
avanicolormap=[255,255,255;141,167,123;128,128,255;35,35,9;170,170,255;246,18,97;222,205,135;255,128,229;90,156,254;74,181,99;255,204,0;111,111,145];
avanicolormap=avanicolormap./255;
%% plotting in a single plot
figure()
superbigmat=zeros(length(newpresortedpre),14);
superbigmat(1:length(newpresortedpre),1:6)=newpresortedpre;
superbigmat(1:length(newpresortedpo),size(superbigmat,2)-5:size(superbigmat,2))=newpresortedpo;
imagesc(superbigmat)
ax=gca;
ax.Colormap=avanicolormap;
box off
xspot=[1.5:size(superbigmat,2)-0.5];
xline(xspot,'color','w','LineWidth',2)
axis off
%% plot for legend
figure()
x=[1:11];
y=[1:11];
hold on
for i=1:length(x)
    scatter(x(i),y(i),'s','filled','MarkerEdgeColor',avanicolormap(i+1,:),'MarkerFaceColor',avanicolormap(i+1,:));
end
ax=gca;
ax.Colormap=avanicolormap(2:end,:);
legend({'Start/End','','g','j','h','k','f','l','c','d',''})