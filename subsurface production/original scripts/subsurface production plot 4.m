% subsurface production
%        Plot 4: cross correlations of mod chl (20-60m) and mod pp (20-60m)
%        Plot: mod chl (20-60m) and mod pp (20-60, as fig 1) NOT REQUIRED



% =======
% PLOT 4A
% =======

% (1) Reading in monthly model variables - PP, PAR, SST

% Only considering 2005
yr = 2005; 
    
lattotal = double(1021);
latstart = double(1);
lontotal = double(1442);
depth = double(9); % top 20 m - only need surface Chl so save memory
pos = 1;


dnomD = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/*d05D.nc', yr); % get a list of D (Diagnostic) files for this year
dnomP = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/*d05P.nc', yr); % get a list of P (passive tracer) files for this year

fnamesD = dir (dnomD);
fnamesP = dir(dnomP);
[num, x] = size(fnamesD);


% Preallocating arrays for model output to be read into
modndpp = NaN(lontotal,lattotal,depth,num);
moddpp = NaN(lontotal,lattotal,depth,num);


% Split PP and chl read in to save memory

% (A) Primary production

    
for f = 1:1:num % loop

    % (2) Primary Production

    fnameD = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/%s', yr, fnamesD(f).name);
    fprintf('Reading %s\n', fnameD);

    % (i) For non-diatom production

    t1 = ncread(fnameD,'PRN3',[1 latstart 1 1], [lontotal lattotal depth 1]); % PRN is depth-integrated so a 2D model field
    t1 = squeeze(t1);
    t1 (t1 == 0) = NaN; % removing zero values (land)

    modndpp(:,:,:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)

    % (ii) For diatom production
    t1 = ncread(fnameD,'PRD3',[1 latstart 1 1], [lontotal lattotal depth 1]); % PRD is depth-integrated so a 2D model field
    t1 = squeeze(t1);
    t1 (t1 == 0) = NaN; % removing zero values (land)

    moddpp(:,:,:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)
    
    pos = pos + 1;

end


clearvars t1

modppC_int = squeeze(nansum(((modndpp(:,:,5:9,:) + moddpp(:,:,5:9,:)) * 12.011 * 6.625),3)); % in mgC. 5:9 = 20-60 m

clearvars modndpp moddpp





% (B) Chlorophyll

pos = 1;


% Preallocating
chn = NaN(lontotal,lattotal,depth,num);
chd = NaN(lontotal,lattotal,depth,num);

for f = 1:1:num % loop
    
    % (1) Chlorophyll

    fnameP = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/%s', yr, fnamesP(f).name);
    fprintf('Reading %s\n', fnameP);

    % (i) For non-diatom chlorophyll
    t1 = ncread(fnameP,'CHN',[1 latstart 1 1], [lontotal lattotal depth 1]); % PRN is depth-integrated so a 2D model field                
    t1 (t1 == 0) = NaN; % removing zero values (land)
    chn(:,:,:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)

    % (ii) For diatom chlorophyll
    t1 = ncread(fnameP,'CHD',[1 latstart 1 1], [lontotal lattotal depth 1]); % PRD is depth-integrated so a 2D model field
    t1 (t1 == 0) = NaN; % removing zero values (land)
    chd(:,:,:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)   


   pos = pos + 1;

end

clearvars t1

modchl_int = squeeze(nansum((chn(:,:,5:9,:) + chd(:,:,5:9,:)),3)); % integrating over depth 20-60m.

clearvars chn chd 





% (2) Calculating cross correlants (over lat-lon grid)


pos = 1;

 for j = 1:size(modppC_int,2) % for every 1/4 lattitude
        for i = 1:size(modppC_int,1) % for every 1/4 longditude
            
            % Annual cross correlation
            a = modppC_int(i,j,:);
            b = modchl_int(i,j,:);
            [c,lags] = xcorr(a,b,0,'coeff'); % zero lag coeff, lags = 0 
            corAn(:,pos) = c;
            
            pos = pos + 1;
            
        end
        fprintf('Reading lat j = %d of 1021\n', j)
 end

 plot4a = reshape(plot4a,1442,1021);
 
 
% (3) Plotting
 
 
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
  
    
    
  
    % Xcorr annual model Chl (20-60m) vs PP (20-60m)
    figure(4);
    
    odvpal(30);
    m_proj('stereographic','lat',90,'long',30,'radius',30)
    m_elev('contour',[-500 -500],'edgecolor','k'); hold on
    m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
    m_pcolor(xx,yy,plot4a); shading flat; hold on; % top 4 boxs ~ top 20m
    m_coast('patch',[.7 .7 .7],'edgecolor','k');
    caxis([0.5 1]); c = colorbar; ylabel(c,'dimensionless'); % should be ~ 0-200 (1000 for whole column)
    title('Cross correlants for 20-60m model Chl and PP','FontWeight','bold')
  
    print -depsc -painters xcorr_modchlpp_20_60m.eps

    
    % complementary: plot mod chl (-20 m) and mod pp (20-60 m; repeat of
    % fig 1b) - won't show much, just that where chl is high, pp is high and vice versa (same as Fig 4) - don't expect a simple linear relationship between chl and pp  
    
    plot5a = nansum(modchl_int,3); % integrate over year
    plot5b = nansum(modppC_int,3); % integrate over year
    
    figure(5);
    odvpal(50);
    subplot(1,2,1);
    m_proj('stereographic','lat',90,'long',30,'radius',30)
    m_elev('contour',[-500 -500],'edgecolor','r'); hold on
    m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
    m_pcolor(xx,yy,plot5a); shading flat; hold on; % top 4 boxs ~ top 20m
    m_coast('patch',[.7 .7 .7],'edgecolor','k');
    caxis([0 150]); c = colorbar; ylabel(c,'dimensionless'); % should be ~ 0-200 (1000 for whole column)
    title('Chl 20-60 m','FontWeight','bold')
    
    subplot(1,2,2);
    m_proj('stereographic','lat',90,'long',30,'radius',30)
    m_elev('contour',[-500 -500],'edgecolor','r'); hold on
    m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
    m_pcolor(xx,yy,plot5b); shading flat; hold on; % top 4 boxs ~ top 20m
    m_coast('patch',[.7 .7 .7],'edgecolor','k');
    caxis([0 1500]); c = colorbar; ylabel(c,'dimensionless'); % should be ~ 0-200 (1000 for whole column)
    title('PP 20-60 m','FontWeight','bold')
    
    print -depsc -painters modchlpp_20_60m.eps
    
    