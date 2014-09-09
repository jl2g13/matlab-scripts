% Subsurface primary production - Plot 1: surface and deep primary production

% created from 'satellite pp monthly' script on 8/8/14 (jl2g13)


% (1) Reading in 5-day model variable - PP   (need 5-day files for XCORR)

% Only considering 2005
yr = 2005; 
lattotal = double(1021);
latstart = double(1);
lontotal = double(1442);
depth = double(14); % read in top 100 m - need for XCORR
pos = 1;

   
dnomD = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/*d05D.nc', yr); % get a list of D (Diagnostic) files for this year
fnamesD = dir (dnomD);
[num, x] = size(fnamesD);


% Preallocating arrays for model output to be read into
modndpp = NaN(lontotal,lattotal,depth,num);
moddpp = NaN(lontotal,lattotal,depth,num);

    
    % Initiating extraction loop
    
for f = 1:1:num % looping from file 1 to num (73) - for each year


    % (2) Primary Production

    fnameD = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/%s', yr, fnamesD(f).name);
    fprintf('Reading %s\n', fnameD);

    % (i) For non-diatom production

    t1 = ncread(fnameD,'PRN3',[1 latstart depth 1], [lontotal lattotal depth 1]); % PRN is depth-integrated so a 2D model field
    t1 = squeeze(t1);
    t1 (t1 == 0) = NaN; % removing zero values (land)

    modndpp(:,:,:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)

    % (ii) For diatom production
    t1 = ncread(fnameD,'PRD3',[1 latstart depth 1], [lontotal lattotal depth 1]); % PRD is depth-integrated so a 2D model field
    t1 = squeeze(t1);
    t1 (t1 == 0) = NaN; % removing zero values (land)

    moddpp(:,:,:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)
    
    pos = pos + 1;

end

    
modppC = (modndpp + moddpp) * 12.011 * 6.625; % in mgC

clearvars modndpp moddpp t1 % clear up memory

modppC_int = squeeze(nansum(modppC,3)); % integrating over depth
modsurfppC = squeeze(modppC(:,:,1,:)); % taking top box = 3m

% ======
% PLOT 1
% ======

ANmodPPC = squeeze(nansum(modppC,4)); % Annual PP (integrating over all months in year)

clearvars modppC

plot1a = nansum(ANmodPPC(:,:,1:4),3); % surface PP = integrate top 4 boxs ~ top 20m
plot1b = nansum(ANmodPPC(:,:,5:9),3); % deep PP = integrate boxes 5--9 ~ 20-60 m depth 

corAn = NaN(1442*1021); % preallocating

pos = 1;
    
 for j = 1:size(modppC_int,2) % for every 1/4 lattitude
        for i = 1:size(modppC_int,1) % for every 1/4 longditude
            
            % Annual cross correlation
            a = modppC_int(i,j,:);
            b = modsurfppC(i,j,:);
            [c,lags] = xcorr(a,b,0,'coeff'); % zero lag coeff, lags = 0 
            corAn(:,pos) = c;
            
            pos = pos + 1;
            
        end
        fprintf('Reading lat j = %d of 1021\n', j)
 end

 corAn = reshape(corAn,1442,1021);

% General plotting preamble
    % Setting up M_Map model
    ft='/noc/altix2/scratch/omfman/ORCA025-N201/means/1990/ORCA025-N201_1990m09P.nc';
    xx=ncread(ft,'nav_lon', [1 1], [lontotal lattotal]);
    yy=ncread(ft,'nav_lat', [1 1], [lontotal lattotal]);

    % removing spurious spikes in N Atl when plotting (beware true zeros, quick fix)
    land = xx(1,1);
    xx(xx == land) = NaN; xx(xx == 0) = NaN;
    land = yy(1,1);
    yy(yy == land) = NaN; yy(yy == 0) = NaN;
  
    
    figure(1);
  
    % annual model PP 0-20 m
    
    subplot(1,3,1);
    m_proj('stereographic','lat',90,'long',30,'radius',30)
    m_elev('contour',[-500 -500],'edgecolor','r'); hold on
    m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
    m_pcolor(xx,yy,plot1a); shading flat; hold on; % top 4 boxs ~ top 20m
    m_coast('patch',[.7 .7 .7],'edgecolor','k');
    caxis([0 500]); c = colorbar; ylabel(c,'mgC/m2/yr'); % should be ~ 0-200 (1000 for whole column)
    title('Annual PP 0-20 m','FontWeight','bold')
    
    % annual model PP 20-60 m

  
    subplot(1,3,2);
    m_proj('stereographic','lat',90,'long',30,'radius',30)
    m_elev('contour',[-500 -500],'edgecolor','r'); hold on
    m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
    m_pcolor(xx,yy,plot1b); shading flat; hold on; % boxes 5--9 ~ 20-60 m depth 
    m_coast('patch',[.7 .7 .7],'edgecolor','k'); 
    caxis([0 500]); c = colorbar; ylabel(c,'mgC/m2/yr'); % should be ~ 0-200
    title('Annual PP 20-60 m','FontWeight','bold')
    
    % annual XCORR between surface (3m) and depth-int PP
    subplot(1,3,3);
    m_proj('stereographic','lat',90,'long',30,'radius',30)
    m_elev('contour',[-500 -500],'edgecolor','r'); hold on
    m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
    m_pcolor(xx,yy,corAn); shading flat; hold on; % XCORR 
    m_coast('patch',[.7 .7 .7],'edgecolor','k'); 
    caxis([0.95 1]); c = colorbar; ylabel(c,'dimensionless');
    title('Annual XCORR (surface vs depth integrated','FontWeight','bold')
    
    
    