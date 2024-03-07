function seq_plot_transitionmatrix(mtx,labels,cm,title1,preorpost)
% cm = 1: one direction for colormap; cm=2 two directions out from 0 

% modified from Lucas in get all trans prob function

f=figure; 
f.Name=[preorpost '_' 'Tmatrix'];
hold on;

syl_num=size(mtx,2);

plotmtx = mtx;
plotmtx(isnan(mtx))=0;
imagesc(plotmtx);

switch cm
    case 1
        colormap(flipud(gray))
    case 2
        % imagesc(mtx,'AlphaData',~isnan(mtx))
        
        colo = [[ones(15,1); linspace(1,0.1,15)'] [linspace(0,1,15)'; linspace(1,0,15)'] [linspace(0.3,1,15)'; ones(15,1)]];
    colormap(colo)

    otherwise
        colormap(gray)
end
textStrings = num2str(mtx(:),'%.0f');
textStrings(isnan(mtx),:)='-';
textStrings = strtrim(cellstr(textStrings));
[x,y]=meshgrid(1:syl_num);
hStrings = text(x(:),y(:),textStrings(:), 'HorizontalAlignment','center','Fontname','arial','fontsize',18,'fontweight','bold');
midValue = mean(get(gca,'CLim'));
textColors = repmat(mtx(:) > midValue,1,3);
set(hStrings,{'Color'},num2cell(textColors,2));
set(gca,'XTick',1:syl_num,'XTickLabel',labels, 'YTick',1:syl_num,...
    'YTickLabel',labels, 'TickLength',[0 0],'fontname','arial','fontweight','bold','fontsize',18);
ylabel('transition from:')
xlabel('transition to:')
title(title1)


