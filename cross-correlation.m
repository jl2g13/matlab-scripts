% Timeseries correlation for PP and Chl (surface vs depth integrated)

% Zero-lag correlation coefficient
% [c,lags] = xcorr(randn(1000,1),rand(1000,1),0,'coeff');
% 	where c = scalar
% 	lags = 0.
% 
% then just plot the scalar

% Description: Pan-Arctic zero-lag cross correlation between surface (3m) and
% depth-integrated (whole column) PP 1-year timeseries (2005). The same for Chlorohpyll
% (3m) cross-correlated with depth-integrated to 95 m.

% This script deals with the months in a clumsy fashion (before I read
% month in as an extra dimension of model field and written at a time when
% I wanted to be using 5-day model output in general) (jl2g13, 7/5/14)

% Chl figure white patch (NaNs) may appear because water depth < depth of model
% input read and then to produce loaded fields - summed across them which for depths > water depths, sum = NaN (should
% NaNsum). If true why don't see in figure in other shallow Arctic shelves?
% (18/5/2014)


% (1) For PP

load /noc/users/jl2g13/MATLAB/model_fields_vars/2005/modppC
load /noc/users/jl2g13/MATLAB/model_fields_vars/2005/modsurfppC

    pos = 1;

    for j = 1:size(modppC,2)
        for i = 1:size(modppC,1) % for every 1/4 longditude
            
            % Jan-Mar
            a = modppC(i,j,1:18);   % Jan-Mar = 1-18, Apr-Jun = 19-36, Jul-Sep = 37-54, Oct-Dec 55-73
            b = modsurfppC(i,j,1:18);
            [c,lags] = xcorr(a,b,0,'coeff'); % zero lag coeff, lags = 0 
            corJM(:,pos) = c; % storing xcorr value for pairs (a,b) at every (1,lon_n)
            
            % Apr-Jun
            a = modppC(i,j,19:36);
            b = modsurfppC(i,j,19:36);
            [c,lags] = xcorr(a,b,0,'coeff'); % zero lag coeff, lags = 0 
            corAJ(:,pos) = c; % storing xcorr value for pairs (a,b) at every (1,lon_n)
            
            % Jul-Sep
            a = modppC(i,j,37:54);
            b = modsurfppC(i,j,37:54);
            [c,lags] = xcorr(a,b,0,'coeff'); % zero lag coeff, lags = 0 
            corJS(:,pos) = c; % storing xcorr value for pairs (a,b) at every (1,lon_n)
            
            % Oct-Dec
            a = modppC(i,j,55:73);
            b = modsurfppC(i,j,55:73);
            [c,lags] = xcorr(a,b,0,'coeff'); % zero lag coeff, lags = 0 
            corOD(:,pos) = c; % storing xcorr value for pairs (a,b) at every (1,lon_n)
            
            % Annual
            a = modppC(i,j,:);
            b = modsurfppC(i,j,:);
            [c,lags] = xcorr(a,b,0,'coeff'); % zero lag coeff, lags = 0 
            cor(:,pos) = c; % storing xcorr value for pairs (a,b) at every (1,lon_n)

            pos = pos + 1;
        end 
        fprintf('Reading %d\n', j)
    end    
        % stores cors in a single vector 1442 * 1021 long (size(cor) = 1 x 1472282) 

     
    % Reshaping cor vectors for plotting

        corJM = reshape(corJM,1442,1021); % taking columnwise vector and reshaping into array with dim = lon x lat = 1442 x 1021
        corJM (corJM == 0) = NaN;
        
        corAJ = reshape(corAJ,1442,1021); % taking columnwise vector and reshaping into array with dim = lon x lat = 1442 x 1021
        corAJ (corAJ == 0) = NaN;
        
        corJS = reshape(corJS,1442,1021); % taking columnwise vector and reshaping into array with dim = lon x lat = 1442 x 1021
        corJS (corJS == 0) = NaN;
        
        corOD = reshape(corOD,1442,1021); % taking columnwise vector and reshaping into array with dim = lon x lat = 1442 x 1021
        corOD (corOD == 0) = NaN;
        
        cor = reshape(cor,1442,1021); % annual
        cor (cor == 0) = NaN;
        
        
    % Preparing quarterly ice overlays

load /noc/users/jl2g13/MATLAB/model_fields_vars/2005/modice
    
        t1 = modice(:,:,1:18); % setting up ice array to add as contour (vars used as above in modice pcolor plots) - note contour = ice min
        iceJM = nansum(t1,3) ./ 18;
        
        t1 = modice(:,:,19:36); % setting up ice array to add as contour (vars used as above in modice pcolor plots) - note contour = ice min
        iceAJ = nansum(t1,3) ./ 18;
        
        t1 = modice(:,:,37:54); % setting up ice array to add as contour (vars used as above in modice pcolor plots) - note contour = ice min
        iceJS = nansum(t1,3) ./ 18;
        
        t1 = modice(:,:,55:73); % setting up ice array to add as contour (vars used as above in modice pcolor plots) - note contour = ice min
        iceOD = nansum(t1,3) ./ 19;
        
        iceendS = modice(:,:,54); % setting up ice array to add as contour (vars used as above in modice pcolor plots) - note contour = ice min


    % M-Map plotting 

        lattotal = double(1021);
        lontotal = double(1442);  
        ft='/noc/altix2/scratch/omfman/ORCA025-N201/means/1990/ORCA025-N201_1990m09P.nc';
        xx=ncread(ft,'nav_lon', [1 1], [lontotal lattotal]);
        yy=ncread(ft,'nav_lat', [1 1], [lontotal lattotal]);

        % removing spurious spikes in N Atl when plotting (beware true zeros, quick fix)
        land = xx(1,1);
        xx(xx == land) = NaN; xx(xx == 0) = NaN;
        land = yy(1,1);
        yy(yy == land) = NaN; yy(yy == 0) = NaN;
        
        
        % corAJ
        figure(1); odvpal(10);
        subplot(2,2,1)
        suptitle('Primary production surface vs depth-integrated cross-correlations');
            m_proj('stereographic','lat',90,'long',30,'radius',30)
            m_elev('contour',[-500 -500],'edgecolor','r'); hold on
            m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
            m_pcolor(xx,yy,corAJ); shading flat; hold on;
            [cs,h]=m_contour(xx,yy,iceAJ,[0.1],'edgecolor','k','LineWidth',1); hold on % 10%  ice cover contour 
            m_coast('patch',[.7 .7 .7],'edgecolor','k');
            caxis([0 1]); cb = colorbar; ylabel(cb,'Normalised cross-correlant')
            title ('Apr-Jun','Fontweight','bold')
           
        % corJS    
        subplot(2,2,2)    
%             figure; odvpal(10);
            m_proj('stereographic','lat',90,'long',30,'radius',30)
            m_elev('contour',[-500 -500],'edgecolor','r'); hold on
            m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
            m_pcolor(xx,yy,corJS); shading flat; hold on;
            [cs,h]=m_contour(xx,yy,iceJS,[0.1],'edgecolor','k','LineWidth',1); hold on % 10%  ice cover contour 
            m_coast('patch',[.7 .7 .7],'edgecolor','k');
            caxis([0 1]); cb = colorbar; ylabel(cb,'Normalised cross-correlant')
            title ('Jul-Sep','Fontweight','bold')
        
        % corOD    
        subplot(2,2,3)    
%             figure; odvpal(10);
            m_proj('stereographic','lat',90,'long',30,'radius',30)
            m_elev('contour',[-500 -500],'edgecolor','r'); hold on
            m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
            m_pcolor(xx,yy,corOD); shading flat; hold on;
            [cs,h]=m_contour(xx,yy,iceOD,[0.1],'edgecolor','k','LineWidth',1); hold on % 10%  ice cover contour 
            m_coast('patch',[.7 .7 .7],'edgecolor','k');
            caxis([0 1]); cb = colorbar; ylabel(cb,'Normalised cross-correlant')
            title ('Oct-Dec','Fontweight','bold')
        
        % annual cor    
        subplot(2,2,4)    
%             figure; odvpal(10);
            m_proj('stereographic','lat',90,'long',30,'radius',30)
            m_elev('contour',[-500 -500],'edgecolor','r'); hold on
            m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
            m_pcolor(xx,yy,cor); shading flat; hold on;
            [cs,h]=m_contour(xx,yy,iceendS,[0.1],'edgecolor','k','LineWidth',1); hold on % 10%  ice cover contour (end September) 
            m_coast('patch',[.7 .7 .7],'edgecolor','k');
            caxis([0 1]); cb = colorbar; ylabel(cb,'Normalised cross-correlant')
            title ('Annual','Fontweight','bold')
            
            % print -depsc -painters PPXcorr.eps
            
            clearvars corJM corAJ corJS corOD
            
% (2) For Chlorophyll 

load /noc/users/jl2g13/MATLAB/model_fields_vars/2005/modchl
load /noc/users/jl2g13/MATLAB/model_fields_vars/2005/modsurfchl

    pos = 1;

    for j = 1:size(modchl,2)
        for i = 1:size(modchl,1) % for every 1/4 longditude
            
%             % Jan-Mar
%             a = modchl(i,j,1:18);  % Jan-Mar = 1-18, Apr-Jun = 19-36, Jul-Sep = 37-54, Oct-Dec 55-73
%             b = modsurfchl(i,j,1:18);
%             [c,lags] = xcorr(a,b,0,'coeff'); % zero lag coeff, lags = 0 
%             corJM(:,pos) = c; % storing xcorr value for pairs (a,b) at every (1,lon_n) 
%             
%             % Apr-Jun
%             a = modchl(i,j,19:36);
%             b = modsurfchl(i,j,19:36);
%             [c,lags] = xcorr(a,b,0,'coeff'); % zero lag coeff, lags = 0 
%             corAJ(:,pos) = c; % storing xcorr value for pairs (a,b) at every (1,lon_n) 
%             
%             % Jul-Sep
%             a = modchl(i,j,37:54);
%             b = modsurfchl(i,j,37:54);
%             [c,lags] = xcorr(a,b,0,'coeff'); % zero lag coeff, lags = 0 
%             corJS(:,pos) = c; % storing xcorr value for pairs (a,b) at every (1,lon_n) 
%             
%             % Oct-Dec
%             a = modchl(i,j,55:73);
%             b = modsurfchl(i,j,55:73);
%             [c,lags] = xcorr(a,b,0,'coeff'); % zero lag coeff, lags = 0 
%             corOD(:,pos) = c; % storing xcorr value for pairs (a,b) at every (1,lon_n) 
%                         
            % Annual
            a = modchl(i,j,:);
            b = modsurfchl(i,j,:);
            [c,lags] = xcorr(a,b,0,'coeff'); % zero lag coeff, lags = 0 
            cor2(:,pos) = c; % storing xcorr value for pairs (a,b) at every (1,lon_n) 

            pos = pos + 1;
        end 
        fprintf('Reading %d\n', j)
    end    
    
    

    % stores cors in a single vector 1442 * 1021 long (size(cor) = 1 x 1472282) 
    
        corJM = reshape(corJM,1442,1021);
        corJM (corJM == 0) = NaN;
        
        corAJ = reshape(corAJ,1442,1021);
        corAJ (corAJ == 0) = NaN;
        
        corJS = reshape(corJS,1442,1021);
        corJS (corJS == 0) = NaN;
        
        corOD = reshape(corOD,1442,1021);
        corOD (corOD == 0) = NaN;
        
        cor2 = reshape(cor2,1442,1021); % taking columnwise vector and reshaping into array with dim = lon x lat = 1442 x 1021
        cor2 (cor2 == 0) = NaN;
        

    % M-Map plotting 
          
            % corJM
%             m_proj('stereographic','lat',90,'long',30,'radius',30)
%             m_elev('contour',[-500 -500],'edgecolor','r'); hold on
%             m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
%             m_pcolor(xx,yy,corJM); shading flat; hold on;
%             [cs,h]=m_contour(xx,yy,iceJM,[0.1],'edgecolor','k','LineWidth',1); hold on % 10%  ice cover contour 
%             m_coast('patch',[.7 .7 .7],'edgecolor','k');
%             caxis([0.7 1]); cb = colorbar; ylabel(cb,'cross-correlant')
%             title ('Annual Chl Xcorr','Fontweight','bold')
    
            
            % corAJ
            figure(2); odvpal(10);
%             figure; odvpal(10);
            subplot(2,2,1)
            suptitle('Chlorophyll surface vs depth-integrated cross-correlations');
            m_proj('stereographic','lat',90,'long',30,'radius',30)
            m_elev('contour',[-500 -500],'edgecolor','r'); hold on
            m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
            m_pcolor(xx,yy,corAJ); shading flat; hold on;
            [cs,h]=m_contour(xx,yy,iceAJ,[0.1],'edgecolor','k','LineWidth',1); hold on % 10%  ice cover contour 
            m_coast('patch',[.7 .7 .7],'edgecolor','k');
            caxis([0.7 1]); cb = colorbar; ylabel(cb,'cross-correlant')
            title ('Apr-Jun','Fontweight','bold')
            
            
            % corJS
%             figure; odvpal(10);
            subplot(2,2,2)
            m_proj('stereographic','lat',90,'long',30,'radius',30)
            m_elev('contour',[-500 -500],'edgecolor','r'); hold on
            m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
            m_pcolor(xx,yy,corJS); shading flat; hold on;
            [cs,h]=m_contour(xx,yy,iceJS,[0.1],'edgecolor','k','LineWidth',1); hold on % 10%  ice cover contour 
            m_coast('patch',[.7 .7 .7],'edgecolor','k');
            caxis([0.7 1]); cb = colorbar; ylabel(cb,'cross-correlant')
            title ('Jul-Sep','Fontweight','bold')
            
            % corOD
%             figure; odvpal(10);
            subplot(2,2,3)
            m_proj('stereographic','lat',90,'long',30,'radius',30)
            m_elev('contour',[-500 -500],'edgecolor','r'); hold on
            m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
            m_pcolor(xx,yy,corOD); shading flat; hold on;
            [cs,h]=m_contour(xx,yy,iceOD,[0.1],'edgecolor','k','LineWidth',1); hold on % 10%  ice cover contour 
            m_coast('patch',[.7 .7 .7],'edgecolor','k');
            caxis([0.7 1]); cb = colorbar; ylabel(cb,'cross-correlant')
            title ('Oct-Dec','Fontweight','bold')
    
    
            % Annual
%             figure; odvpal(10);
            subplot(2,2,4)
            m_proj('stereographic','lat',90,'long',30,'radius',30)
            m_elev('contour',[-500 -500],'edgecolor','r'); hold on
            m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
            m_pcolor(xx,yy,cor2); shading flat; hold on;
            [cs,h]=m_contour(xx,yy,iceendS,[0.1],'edgecolor','k','LineWidth',1); hold on % 10%  ice cover contour 
            m_coast('patch',[.7 .7 .7],'edgecolor','k');
            caxis([0.7 1]); cb = colorbar; ylabel(cb,'cross-correlant')
            title ('Annual','Fontweight','bold')
            
           % print -depsc -painters ChlXcorr.eps
            
            
% (3) Annual PP and Chl on same figure       


   % annual   
    
    figure; odvpal(10);
    
    % PP
    subplot(2,1,1)
    suptitle('Annual surface vs depth-integrated cross-correlations')
    m_proj('stereographic','lat',90,'long',30,'radius',30)
    m_elev('contour',[-500 -500],'edgecolor','r'); hold on
    m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
    m_pcolor(xx,yy,cor); shading flat; hold on;
    [cs,h]=m_contour(xx,yy,iceendS,[0.1],'edgecolor','k','LineWidth',1); hold on % 10%  ice cover contour (end September) 
    m_coast('patch',[.7 .7 .7],'edgecolor','k');
    caxis([0 1]); cb = colorbar; ylabel(cb,'Normalised cross-correlant')
    title ('Primary production','Fontweight','bold')
    
    % chl
    subplot(2,1,2)    
    m_proj('stereographic','lat',90,'long',30,'radius',30)
    m_elev('contour',[-500 -500],'edgecolor','r'); hold on
    m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
    m_pcolor(xx,yy,cor2); shading flat; hold on;
    [cs,h]=m_contour(xx,yy,iceendS,[0.1],'edgecolor','k','LineWidth',1); hold on % 10%  ice cover contour (end September) 
    m_coast('patch',[.7 .7 .7],'edgecolor','k');
    caxis([0 1]); cb = colorbar; ylabel(cb,'Normalised cross-correlant')
    title ('Chlorophyll','Fontweight','bold')
    
   % print -depsc -painters AnnXcorr.eps