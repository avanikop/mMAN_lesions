%% fig3 sup1 a
gapmedian_bird1=csvread('D:\analysis\data_for_elife_mMAN\Figure3_figuresupplement1_sourcedata1\bird1_gap.csv')
gapmedian_bird2=csvread('D:\analysis\data_for_elife_mMAN\Figure3_figuresupplement1_sourcedata1\bird2_gap.csv')
gapmedian_bird3=csvread('D:\analysis\data_for_elife_mMAN\Figure3_figuresupplement1_sourcedata1\bird3_gap.csv')
gapmedian_bird4=csvread('D:\analysis\data_for_elife_mMAN\Figure3_figuresupplement1_sourcedata1\bird4_gap.csv')
gapmedian_bird5=csvread('D:\analysis\data_for_elife_mMAN\Figure3_figuresupplement1_sourcedata1\bird5_gap.csv')
gapmedian_bird6=csvread('D:\analysis\data_for_elife_mMAN\Figure3_figuresupplement1_sourcedata1\bird6_gap.csv')
gapmedian_bird7=csvread('D:\analysis\data_for_elife_mMAN\Figure3_figuresupplement1_sourcedata1\bird7_gap.csv')

delTE_bird1=csvread('D:\analysis\data_for_elife_mMAN\Figure3_figuresupplement1_sourcedata1\bird1_TE.csv')
delTE_bird2=csvread('D:\analysis\data_for_elife_mMAN\Figure3_figuresupplement1_sourcedata1\bird2_TE.csv')
delTE_bird3=csvread('D:\analysis\data_for_elife_mMAN\Figure3_figuresupplement1_sourcedata1\bird3_TE.csv')
delTE_bird4=csvread('D:\analysis\data_for_elife_mMAN\Figure3_figuresupplement1_sourcedata1\bird4_TE.csv')
delTE_bird5=csvread('D:\analysis\data_for_elife_mMAN\Figure3_figuresupplement1_sourcedata1\bird5_TE.csv')
delTE_bird6=csvread('D:\analysis\data_for_elife_mMAN\Figure3_figuresupplement1_sourcedata1\bird6_TE.csv')
delTE_bird7=csvread('D:\analysis\data_for_elife_mMAN\Figure3_figuresupplement1_sourcedata1\bird7_TE.csv')
%% statistics
gapmedian_pre=[gapmedian_bird1(:,1);gapmedian_bird2(:,1);gapmedian_bird3(:,1);...
    gapmedian_bird4(:,1);gapmedian_bird5(:,1);gapmedian_bird6(:,1);gapmedian_bird7(:,1)];
delTE=[delTE_bird1;delTE_bird2;delTE_bird3;delTE_bird4;...
    delTE_bird5;delTE_bird6;delTE_bird7];
[correlation,pval]=corr(gapmedian_pre,delTE);
% calculating and plotting regression line
x=gapmedian_pre;
X=[ones(length(x),1),x];
y=delTE;
b=X\y;
yCalc1=X*b;
Rsq=1-sum((y-yCalc1).^2)/sum((y - mean(y)).^2);
%% plotting
figure('Name','gap duration vs change in TE')
title('del Transition entropy vs gap duration')
hold on
scatter(gapmedian_bird1(:,1),delTE_bird1,'o')
scatter(gapmedian_bird2(:,1),delTE_bird2,'+')
scatter(gapmedian_bird3(:,1),delTE_bird3,'x')
scatter(gapmedian_bird4(:,1),delTE_bird4,'s')
scatter(gapmedian_bird5(:,1),delTE_bird5,'d')
scatter(gapmedian_bird6(:,1),delTE_bird6,'^')
scatter(gapmedian_bird7(:,1),delTE_bird7,'p')
yline(0,'k--')
xlabel('Pre lesion gap duration (ms)')
ylabel('Change in transition entropy within chunk')
plot(x,yCalc1,'--k')
txt1={['Corr.coef.=',num2str(correlation)],['p=',num2str(round(pval,3))],['Rsq=',num2str(round(Rsq,3))]}
text(90,1,txt1);
%text(90,1.1,txt2)
hold off
