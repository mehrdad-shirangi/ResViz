function plotAllperms(sol,iStep,Ne, jMAP, Nx,Ny,maxColor,minColor)
% This function plots all the permeability fields at a given iStep
% There are Ne RML realizations and MAP estimate (total Ne + 1
% realizations)
% Written by Mehrdad Gharib Shirangi, May 
tStep = 210;

close all

set(0,'defaultaxesfontsize',14);
set(0,'defaulttextfontsize',14);

FontSize = 14;
% ================================================
% Nx = 50; Ny = 50; 
Nz = 1; maxnwells = 8;
binWells = sol(1:maxnwells);
locs     = sol(maxnwells+1:3*(maxnwells));
locs     = reshape(locs,2,length(locs)/2)';
locs = [locs zeros(length(binWells),1)];    % Added column for drill order
injLocs = locs((binWells==-1),:);
prodLocs = locs((binWells==1),:);
ninj = size(injLocs,1);
nprod = size(prodLocs,1);
% Add drill order
drillNumber = 0;
injIndx    = 0;
prodIndx   = 0;
for i = 1:length(binWells)
    if binWells(i) < 0
        drillNumber = drillNumber + 1;
        injIndx     = injIndx     + 1;
        injLocs(injIndx,end) = drillNumber;
        locs(i,end) =drillNumber; 
    elseif binWells(i) > 0
        drillNumber = drillNumber + 1;
        prodIndx    = prodIndx    + 1;
        prodLocs(prodIndx,end) = drillNumber;
        locs(i,end) =drillNumber;         
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
if ( Ne < 6)
    M = 2; N = 3; % dimension of the subplot    
    set(gcf,'units','normalized','outerposition',[0 0 0.55 0.72]);
elseif ( Ne < 7)
    M = 3; N = 3; % dimension of the subplot    
    set(gcf,'units','normalized','outerposition',[0 0 0.55 1]);    
else
    set(gcf,'units','normalized','outerposition',[0 0 0.75 1]);        
    M = 3; N = 4; % dimension of the subplot
end

%   s= 'permx.in';
  s = 'mtrue.in';
  m = load(s); %m = log(m);
  figure(1); 
 if ( Ne < 7)  
    subplot(M,N,1)
 else
    subplot(M,N,2)     
 end
  plot2Dperm(m,Nx,Ny,well,welltype,FontSize,maxColor,minColor); %colorbar('Fontsize',18)
  title('Truth','FontSize',22)
      

failCounter = 0;  
if (iStep <2)   % plot prior realizations 
    iStep = 2;
    i = 0;
    for jens = [1:Ne]
%         disp('yes')
          i = i + 1 ;
          s =['results.' num2str(iStep) '/mprior' num2str(iStep) '_' num2str(jens) '.in'];  
% -----------------------------          
          fid = fopen(s, 'r');
          if fid == -1
              failCounter = failCounter + 1;  
              s =['results.' num2str(iStep) '/mprior' num2str(iStep) '_' num2str(Ne+failCounter) '.in'];                
%               continue
          end
%           else
            fclose(fid);
%           end
% -----------------------------
          m = load(s);
          figure(1); 
          if ( Ne < 6)
             subplot(M,N,1+i)     
          else
             subplot(M,N,2+i)
          end
          plot2Dperm(m,Nx,Ny,[],[],FontSize,maxColor,minColor); %colorbar('Fontsize',18)      
          if jens == jMAP
              title(['Prior Mean'],'FontSize',22)
          else
              title(['Prior # ' num2str(jens)],'FontSize',22)
          end
    end 
    iStep = 1;     
else    
    i = 0;    
    for jens = [1:Ne]
          i = i + 1 ;        
          s =['results.' num2str(iStep) '/mcond' num2str(iStep) '_' num2str(jens) '.out'];
% -----------------------------          
          fid = fopen(s, 'r');
          if fid == -1
              failCounter = failCounter + 1;  
              disp(['File ' s ' not found \n'])
              s =['results.' num2str(iStep) '/mcond' num2str(iStep) '_' num2str(Ne+failCounter) '.out'];                
%               continue
          else
              fclose(fid);             
          end
% -----------------------------  % Added April 7 2014
          m = load(s);
          figure(1); 
          if ( Ne < 6)
             subplot(M,N,1+i)     
          else
             subplot(M,N,2+i)
          end          
%           plot2Dperm(m,Nx,Ny,well(1:iStep-1,:),welltype(1:iStep-1),FontSize); %colorbar('Fontsize',18)
          plot2Dperm(m,Nx,Ny,well(1:min(iStep,size(well,1)),:),welltype(1:min(iStep,length(welltype))),FontSize,maxColor,minColor); %colorbar('Fontsize',18)
          if jens == jMAP
              title(['MAP'],'FontSize',22)
          else
              title(['RML # ' num2str(jens)],'FontSize',22)
          end
    end
end 
t = (iStep-1) * tStep;
%   figure(2)
%   s ='MAP--';
%   s1 = [s '.fig'];   s2 = [s '.emf'];    s3 = [s '.pdf']
%   saveas2(s1); saveas2(s2); saveas2(s3);
% 
  figure(1)
  s =['allreal'  num2str(iStep)];
  s1 = [s '.fig'];   s2 = [s '.emf'];    s3 = [s '.pdf'];
  saveas2(s1); saveas(gcf,s,'emf'); saveas2(s3);


%     %% Plotting each one separately
