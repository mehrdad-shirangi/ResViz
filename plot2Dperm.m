% Plot 2D maps 
% June 7 2012
function plot2Dperm(m,nx, ny,well,welltype,fontsize,maxColor,minColor)
%Example: plot2Dperm(Ensem001,40,40,[5,5],[1])
% Example 2:
% uiopen('C:\Users\Mehrdad\Dropbox\Codes & Results\input files - well placement\4\4 Results-June 7- until 840 Days\Updated_Ensembles\Ensem002.DAT',1)
% plot2Dperm(Ensem002,40,40,[3,11;13,5;30,7],[1,-1,1]); colorbar('Fontsize',18)
% saveas2('4real.fig'); saveas2('4real.pdf');saveas2('4real.emf')
% Example 3:
% plot2Dperm(Ensem001,40,40,[3,11;13,5;30,7],[1,-1,1]); colorbar('Fontsize',18)
% saveas2('4map.fig'); saveas2('4map.pdf');saveas2('4map.emf')
% Example 4:
% plot2Dperm(Ensem001,40,40,[3,11;13,5;30,7;4,26],[1,-1,1,-1]); colorbar('Fontsize',18)
% saveas2('5map.fig'); saveas2('5map.pdf');saveas2('5map.emf')
% plot2Dperm(Ensem002,40,40,[3,11;13,5;30,7;4,26],[1,-1,1,-1]); colorbar('Fontsize',18)
% saveas2('5real.fig'); saveas2('5real.pdf');saveas2('5real.emf')
% 
% N = 5;
% well = [3,11;13,5;30,7;4,26];
% welltype = [1,-1,0,1,-1];
% figure(1); figure(2);
% for i = 1 : N
%   s =sprintf('%dEnsem002.dat',i)
%   m = load(s);
%   figure(1); subplot(2,3,i+1)
%   plot2Dperm(m,40,40,well(1:i,:),welltype(1:i),colorbar('Fontsize',18)
%   s = sprintf('%dreal',i);
%   s1 = [s '.fig'];   s2 = [s '.emf'];    s3 = [s '.pdf']
%   saveas2(s1); saveas2(s2); saveas2(s3);
%   s =sprintf('%dEnsem001.dat',i)
%   m = load(s);
%   figure(2); subplot(2,3,i+1)
%   plot2Dperm(m,40,40,well(1:i,:),welltype(1:i),colorbar('Fontsize',18)
%   s = sprintf('%dmap',i);
%   s1 = [s '.fig'];   s2 = [s '.emf'];    s3 = [s '.pdf']
%   saveas2(s1); saveas2(s2); saveas2(s3);
% end
y = reshape(m,nx,ny)';
imagesc(y)
[nw , s] = size(well);
caxis([minColor maxColor]);
% set(gca,'CLim',[-2 7.5]);
inj = 1;
prd = 1;
hard = 1;
hold on
for i = 1 : nw
    if welltype(i) ==-1 % injector
        s = sprintf(' I-%d',inj);
%         text(well(i,1) - 0.5,well(i,2),['\otimes' s],'FontSize',fontsize,'FontWeight','Bold'); % ,'color','white'
%         plot(well(i,1),well(i,2),'kO','MarkerSize',14,'MarkerFaceColor','k')
        plot(well(i,1) - 0.5,well(i,2),'vk','MarkerSize',fontsize,'MarkerFaceColor','k')
        text(well(i,1) + 1,well(i,2),[s],'FontSize',fontsize,'FontWeight','Bold'); % ,'color','white'
        inj = inj+1;
    elseif  welltype(i) == 1 % producer
        s = sprintf(' P-%d',prd);
%         text(well(i,1) - 0.5,well(i,2),['\oslash' s],'FontSize',fontsize,'FontWeight','Bold');        
%         scatter(well(i,1) - 0.5,well(i,2),'circle','filled','color','black')
        plot(well(i,1) - 0.5,well(i,2),'ok','MarkerSize',fontsize,'MarkerFaceColor','k')
        text(well(i,1) + 1,well(i,2),[s],'FontSize',fontsize,'FontWeight','Bold'); % ,'color','white'
        prd = prd + 1;
    elseif  welltype(i) == 2 % HardData
        s = sprintf(' H-%d',hard);
        text(well(i,1) - 0.5,well(i,2),['\o' s],'FontSize',fontsize,'FontWeight','Bold');        
        hard = hard + 1;        
    end
end

