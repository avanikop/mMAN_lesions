function [h,chi2,p]=chisq_2dist(dist1,dist2,alpha)
%chisq gof
%source= https://www.youtube.com/watch?v=Cert35F-w4c
%dist1=observed values
%dist2=expectedvalues
num=find(dist1==0);
dist1(num)=[]; %removing zeros based on obs values
dist2(num)=[]; %removing zeros based on obs values
if isrow(dist1)
    dist1=dist1';
end
if isrow(dist2)
    dist2=dist2';
end
df=length(dist1)-1;
obs=[dist1,dist2];
rowtot=sum(obs,2); %along columns
coltot=sum(obs,1);
tot=sum(sum(obs));
exp=zeros(size(obs));
for i=1:size(obs,1)
    for j=1:size(obs,2)
        exp(i,j)=rowtot(i)*coltot(j)/(tot);
    end
end
diff=obs-exp;
diff2=diff.^2;
normdiff2=diff2./exp;
chi2=sum(sum(normdiff2));
p=1-chi2cdf(chi2,df);
if (~exist('alpha','var'))
    alpha=0.01;
end
if p<alpha
    h=1;
else
    h=0;
end
end