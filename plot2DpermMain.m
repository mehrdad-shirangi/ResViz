function plot2DpermMain(sol,jens,Nx,Ny,maxColor,minColor,tStep,iStepVec,maxnwells,NwDrill)
% plot2DpermMain(sol,jens,Nx,Ny)
% plots the evolution of permeability field for realization # j 
% 

close all

set(0,'defaultaxesfontsize',14);
set(0,'defaulttextfontsize',14);

FontSize = 14;
% ================================================
Nz = 1; 
% maxnwells = 8;
binWells = sol(1:maxnwells);
locs     = sol(maxnwells+1:3*(maxnwells));
locs     = reshape(locs,2,length(locs)/2)';
locs = [locs zeros(length(binWells),1)];    % Added column for drill order
injLocs = locs((binWells==-1),:);
prodLocs = locs((binWells==1),:);
ninj = size(injLocs,1);
nprod = size(prodLocs,1);
injLocs = [injLocs, zeros(ninj,1)];
prodLocs = [prodLocs, zeros(nprod,1)];

% Add drill order
drillNumber = 0;
injIndx    = 0;
prodIndx   = 0;
for i = 1:length(binWells)
    if binWells(i) < 0
        drillNumber = drillNumber + 1;
        injIndx     = injIndx     + 1;
        injLocs(injIndx,3) = drillNumber;
        locs(i,end) =drillNumber; 
        injLocs(injIndx,4) = i;                
    elseif binWells(i) > 0
        drillNumber = drillNumber + 1;
        prodIndx    = prodIndx    + 1;
        prodLocs(prodIndx,3) = drillNumber;
        locs(i,end) =drillNumber;      
        prodLocs(prodIndx,4) = i;                
    end
end
% nctrlpers = length(ctrlPeriodVec);
% ctrls = sol(3*(maxnwells)+1:end);
% ctrls = reshape(ctrls,maxnwells,nctrlpers);
% 
% injCtrls  = ctrls((binWells==-1),:);
% prodCtrls = ctrls((binWells==1),:);

% ================================================
        
well = locs; 
welltype =binWells; 

Hard = [];
Hardtype = [];  %[2,2,2];
figure(1);

%     set(gcf,'units','normalized','outerposition',[0 0 1 1]);

M = 2; N = 3; % dimension of the subplot

% if N == 4    
%     set(gcf,'units','normalized','outerposition',[0 0 0.75 1]);    
% else
%     set(gcf,'units','normalized','outerposition',[0 0 0.55 1]);    
    set(gcf,'units','normalized','outerposition',[0 0 0.55 0.7]);    
% end

  s = 'mtrue.in';
  m = load(s); 
  figure(1); subplot(M,N,1)
  plot2Dperm(m,Nx,Ny,well,welltype,FontSize,maxColor,minColor); %colorbar('Fontsize',18)
  title('Truth','FontSize',22)
      


  iStep = iStepVec(1);
  s =['results.' num2str(iStep) '/mprior' num2str(iStep) '_' num2str(jens) '.in'];  
  m = load(s);
  figure(1); subplot(M,N,2)
  plot2Dperm(m,Nx,Ny,[],[],FontSize,maxColor,minColor); %colorbar('Fontsize',18)
  title('0 days','FontSize',22)

  for iStep = iStepVec
      s =['results.' num2str(iStep) '/mcond' num2str(iStep) '_' num2str(jens) '.out'];
      m = load(s);
      figure(1); subplot(M,N,1+iStep)
%       plot2Dperm(m,Nx,Ny,well(1:iStep-1,:),welltype(1:iStep-1),FontSize); %colorbar('Fontsize',18)
      plot2Dperm(m,Nx,Ny,well(1:min(iStep*NwDrill,size(well,1)),:),welltype(1:min(iStep*NwDrill,length(welltype))),FontSize,maxColor,minColor); %colorbar('Fontsize',18)
%       t = sum(binWells(1:iStep-1)~=0) * tStep;
      t = (iStep-1) * tStep;
      title(sprintf('%d Days',t),'FontSize',22)
    %   s = sprintf('%dreal',i);
    %   s1 = [s '.fig'];   s2 = [s '.emf'];    s3 = [s '.pdf']
    %   saveas2(s1); saveas2(s2); saveas2(s3);
  end

%   figure(2)
%   s ='MAP--';
%   s1 = [s '.fig'];   s2 = [s '.emf'];    s3 = [s '.pdf']
%   saveas2(s1); saveas2(s2); saveas2(s3);
% 
  figure(1)
  s =['evol_real'  num2str(jens)];
  s1 = [s '.fig'];   s2 = [s '.emf'];    s3 = [s '.pdf'];
  saveas2(s1); saveas(gcf,s,'emf'); saveas2(s3);

