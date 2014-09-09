% Subsurface primary production - Plot 1: surface and deep primary
% production and their cross-correlation
    % Fig 1: 0-20 and 20-60 PP annual and for post bloom (Aug-Sep)
    % Fig 2: Subsurface PP proportion by season
    % Fig 3: Subsurface PP proportion by month (May-Aug)
    
    
% Needs to be run on linux machine: iapetus - otherwise runs out of memory    

% created from 'satellite pp monthly' script on 8/8/14 (jl2g13)


% (1) Reading in 5-day model variable - PP   (need 5-day files for XCORR)

% Only considering 2005
yr = 2005;

lattotal = double(1021);
latstart = double(1);
lontotal = double(1442);
depth = double(10); % 14 = read in top 100 m - need for XCORR (could read in depth-int PP var and depth = 1)
pos = 1;

   
dnomD = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/*d05D.nc', yr); % get a list of D (Diagnostic) files for this year
fnamesD = dir (dnomD);
[num, x] = size(fnamesD);


% Preallocating arrays for model output to be read into
modndpp = NaN(lontotal,lattotal,depth,num);
moddpp = NaN(lontotal,lattotal,depth,num);

    
    % Initiating extraction loop
    
for f = 1:1:num % looping from file 1 to num (73) - for each year


    % (1) Primary Production

    fnameD = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/%s', yr, fnamesD(f).name);
    fprintf('Reading %s\n', fnameD);

    % (i) For non-diatom production

    t1 = ncread(fnameD,'PRN3',[1 latstart depth 1], [lontotal lattotal depth 1]); % PRN is depth-integrated so a 2D model field
    t1 = squeeze(t1);
     t1 (t1 == 0) = NaN; % removing zero values (land). some actual zeros
%     in central basin

    modndpp(:,:,:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)

    % (ii) For diatom production
    t1 = ncread(fnameD,'PRD3',[1 latstart depth 1], [lontotal lattotal depth 1]); % PRD is depth-integrated so a 2D model field
    t1 = squeeze(t1);
     t1 (t1 == 0) = NaN; % removing zero values (land). some actual zeros
%     in central basin

    moddpp(:,:,:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)
    
   
    % check values reasonable and most PP in 60 m: read in depth-integrated PP
    
    t1 = ncread(fnameD,'PRD',[1 latstart 1], [lontotal lattotal 1]);
    t1 = squeeze(t1);
    t1 (t1 == 0) = NaN; 
    
    prd(:,:,pos) = t1;
    
    t1 = ncread(fnameD,'PRN',[1 latstart 1], [lontotal lattotal 1]);
    t1 = squeeze(t1);
    t1 (t1 == 0) = NaN; 
    
    prn(:,:,pos) = t1;
    
    pos = pos + 1;
    
end

modppC = (modndpp + moddpp) * 12.011 * 6.625; % in mgC

% remove anywhere water depth < 60 m -> NaN: take all depths down to filter depth in column; sum over
% them, if any of them are NaN (rock) then t1->0 at this point -> NaN.
t1 = modppC(:,:,:,1) * 0; % only want filter at top 60 m -> not over whole grid. nb. n=10, is z~60m. 
t2 = sum(isfinite(t1),3); % chl = finite in water, = NaN below seafloor
kbathy = t2; % model water depth in terms of k model levels
kbathy (kbathy == 0) = NaN; % from here just times whatever model field you want to filter by kbathy
kbathy (kbathy >= 1) = 1; % when sums over finite values each x-y water column gets value 1-10 if > 60m deep (would be 10 if didn't set == 0 -> NaN at read in above)


pp = (prn + prd) * 12.011 * 6.625; % mgC


clearvars modndpp moddpp t1 prn prd % clear up memory


% % Quick check fields read in correctly
% t1 = nansum(modppC_int,3);
% t2 = nansum(modsurfppC,3);
% 
% 
% % t1/t2
%  m_proj('stereographic','lat',90,'long',30,'radius',30)
%  m_elev('contour',[-500 -500],'edgecolor','r'); hold on
%  m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
%  m_pcolor(xx,yy,t2); shading flat; hold on;
%  m_coast('patch',[.7 .7 .7],'edgecolor','k'); 
%  caxis([0 0.5]); c = colorbar; ylabel(c,'');
%  title('','FontWeight','bold')
%  
%  % -> look like they're *1000 too low (i.e. in g not mg)



% ======
% PLOT 1
% ======
pp_an = (nansum(pp,3)) .* (365/73);

modppC_an = (nansum(modppC,4)) .* (365/73); % Annual PP (integrating over all months in year). days in yr/no files. mgC/m3/yr
    
    
% clearvars modppC

plot1a_t = squeeze(nansum(modppC_an(:,:,1:4),3) * (20/4)); % surface PP = integrate top 4 boxs ~ top 20m. depth-interval/noboxes. mgC/m2/yr
plot1b_t = squeeze(nansum(modppC_an(:,:,5:10),3) * (40/6)); % deep PP = integrate boxes 5-10 ~ 20-60 m depth. mgC/m2/yr

prop_t = plot1b_t ./ (plot1a_t + plot1b_t);

plot1a = plot1a_t .* kbathy; % kbathy time indepedent: = 1 where water > 20 deep, = 0 where water < 20m deep --- removing places with water < 20m deep
plot1b = plot1b_t .* kbathy; % as above but relevant depth = 60 m
prop = prop_t .* kbathy;

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
    
    
    % check plot: pp_an - depth-int annual-int pp field
    
    figure;
    odvpal(50);
    
    m_proj('stereographic','lat',90,'long',30,'radius',30)
    m_elev('contour',[-500 -500],'edgecolor','r'); hold on
    m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
    m_pcolor(xx,yy,pp_an); shading flat; hold on;
    m_coast('patch',[.7 .7 .7],'edgecolor','k');
    caxis([0 100000]); cb = colorbar; ylabel(cb,'mgC/m2/yr'); % 0-100 gC/m2/yr
    

    % annual model PP 0-20 m
    figure(1);
    odvpal(50);
    
     subplot(3,1,1);
    m_proj('stereographic','lat',90,'long',30,'radius',30)
    m_elev('contour',[-500 -500],'edgecolor','r'); hold on
    m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
    m_pcolor(xx,yy,plot1a); shading flat; hold on; % top 4 boxs ~ top 20m
    m_coast('patch',[.7 .7 .7],'edgecolor','k');
    caxis([0 1000]); cb = colorbar; ylabel(cb,'mgC/m2/yr');
%     title('Annual PP 0-20 m','FontWeight','bold')
    title('a.','position',[-0.5 0.6],'FontWeight','bold')
    
  
    
    % annual model PP 20-60 m

    subplot(3,1,2);
    m_proj('stereographic','lat',90,'long',30,'radius',30)
    m_elev('contour',[-500 -500],'edgecolor','r'); hold on
    m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
    m_pcolor(xx,yy,plot1b); shading flat; hold on; % boxes 5--9 ~ 20-60 m depth 
    m_coast('patch',[.7 .7 .7],'edgecolor','k'); 
    caxis([0 100]); cb = colorbar; ylabel(cb,'mgC/m2/yr');
%     title('Annual PP 20-60 m','FontWeight','bold')
    title('b.','position',[-0.5 0.6],'FontWeight','bold')  % 'position': manually set title position using [x y] vector 
    
    
    % proportion of model PP that is subsurface (20-60 m)
    
    subplot(3,1,3);
    m_proj('stereographic','lat',90,'long',30,'radius',30)
    m_elev('contour',[-500 -500],'edgecolor','r'); hold on
    m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
    m_pcolor(xx,yy,prop); shading flat; hold on; % boxes 5--9 ~ 20-60 m depth 
    m_coast('patch',[.7 .7 .7],'edgecolor','k'); 
    caxis([0 0.15]); cb = colorbar; ylabel(cb,'dimensionless');
%     title('Annual proportion 20-60 m','FontWeight','bold')
    title('c.','position',[-0.5 0.6],'FontWeight','bold')
    
    
    print -depsc -painters annualPP_20_60_prop.eps
    
    
    
    % manually set figure layout
    
    
    figure;
    odvpal(50);
    
    positionVector1 = [0.15, 0.6, 0.3, 0.4]; % first two [left bottom width height] (values 0-1) . left-bottom positions bottom left corner of figure, width-height to set size
    subplot('Position',positionVector1)

    m_proj('stereographic','lat',90,'long',30,'radius',30)
    m_elev('contour',[-500 -500],'edgecolor','r'); hold on
    m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
    m_pcolor(xx,yy,plot1a); shading flat; hold on; % top 4 boxs ~ top 20m
    m_coast('patch',[.7 .7 .7],'edgecolor','k');
    caxis([0 1000]); cb = colorbar; ylabel(cb,'mgC/m2/yr');
%     title('Annual PP 0-20 m','FontWeight','bold')
    title('a.','position',[-0.5 0.6],'FontWeight','bold')

    
    
    positionVector2 = [0.55, 0.6, 0.3, 0.4];
    subplot('Position',positionVector2)
    
    m_proj('stereographic','lat',90,'long',30,'radius',30)
    m_elev('contour',[-500 -500],'edgecolor','r'); hold on
    m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
    m_pcolor(xx,yy,plot1b); shading flat; hold on; % boxes 5--9 ~ 20-60 m depth 
    m_coast('patch',[.7 .7 .7],'edgecolor','k'); 
    caxis([0 100]); cb = colorbar; ylabel(cb,'mgC/m2/yr');
%     title('Annual PP 20-60 m','FontWeight','bold')
    title('b.','position',[-0.5 0.6],'FontWeight','bold')     
    
    
    positionVector3 = [0.3, 0.1, 0.4, 0.5];
    subplot('Position',positionVector3)
    
    m_proj('stereographic','lat',90,'long',30,'radius',30)
    m_elev('contour',[-500 -500],'edgecolor','r'); hold on
    m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
    m_pcolor(xx,yy,prop); shading flat; hold on; % boxes 5--9 ~ 20-60 m depth 
    m_coast('patch',[.7 .7 .7],'edgecolor','k'); 
    caxis([0 0.15]); cb = colorbar; ylabel(cb,'dimensionless');
%     title('Annual proportion 20-60 m','FontWeight','bold')
    title('c.','position',[-0.5 0.6],'FontWeight','bold')
    
    
    print -depsc -painters fig1.eps
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
% ========
% PLOT 2    
% ========
modppC_int = squeeze(nansum(modppC,3)); % integrating over depth. mgC/m2/d
modsurfppC = squeeze(modppC(:,:,1,:) * 3); % taking top box = 3m. mgC/m2/d

pos = 1;

 for j = 1:size(modppC,2) % for every 1/4 lattitude
        for i = 1:size(modppC,1) % for every 1/4 longditude
            
            % Annual cross correlation
            a = modppC_int(i,j,:);
            b = modsurfppC(i,j,:);
            [c,lags] = xcorr(a,b,0,'coeff'); % zero lag coeff, lags = 0 
            corAn(:,pos) = c;
            
            pos = pos + 1;
            
        end
        fprintf('Reading lat = %d of 1021\n', j)
 end

plot1c = reshape(corAn,1442,1021);    
    
    % annual XCORR between surface (3m) and depth-int PP
    figure(2);
    
    % Issue - plotting nans in water depth < 100 m (depth = 14) and in
    % central basin (perhaps PP = 0 at <= 100 m depth -> NaN so corAn -> NaN)
    
    m_proj('stereographic','lat',90,'long',30,'radius',30)
    m_elev('contour',[-100 -100],'edgecolor','r'); hold on
    m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
    m_pcolor(xx,yy,plot1c); shading flat; hold on; % XCORR 
    m_coast('patch',[.7 .7 .7],'edgecolor','k'); 
    caxis([0 1]); cb = colorbar; ylabel(cb,'units');
    title('Annual XCORR (surface vs depth integrated','FontWeight','bold')
    
    print -depsc -painters xcorr_modPP_surf_vs_col.eps
    
    
    % SUPPLEMENTARY FIGURE - Annual pecentage of production that happens 20-60 m. PLOT1A/PLOT1B
%     
%     figure(3);
%     odvpal(50);
%     
%     supplot = plot1b ./ (plot1b + plot1b);
%     
%     % Issue - plotting nans in water depth < 100 m (depth = 14) and in
%     % central basin (perhaps PP = 0 at <= 100 m depth -> NaN so corAn -> NaN)
%     
%     m_proj('stereographic','lat',90,'long',30,'radius',30)
%     m_elev('contour',[-100 -100],'edgecolor','r'); hold on
%     m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
%     m_pcolor(xx,yy,supplot); shading flat; hold on;
%     m_coast('patch',[.7 .7 .7],'edgecolor','k'); 
%     caxis([0 20]); cb = colorbar; ylabel(cb,'dimensionless');
%     title('Plot1a/plot1b','FontWeight','bold')
%     
%     print -depsc -painters plot1a_div_plot1b.eps
%     
%     
%     


    % SUPPLEMENTARY FIGURE - Autumn (Jul-Sep) PP: surface, subsurface and percentage subsurface
    
    
    
    sur_depth = 20; % surface depth range (0-20 m)
    sur_files = 4; % no of boxes in surface depth range
    sub_depth = 40; % subsurface depth range (20-60 m)
    sub_files = 5; % number of boxes in subsurface depth range
    
    JS_files = 18; % number of files in Jul-Sep
    JS = 92; % no of days in Jul-Sep
    
    
    subsurf_JS = (nansum(nansum(modppC(:,:,5:9,37:54),4),3)) .* (sub_depth * JS / sub_files * JS_files); % mgC/m2/season
    surf_JS = (nansum(nansum(modppC(:,:,1:4,37:54),4),3)) .* (sur_depth * JS / sur_files * JS_files); % mgC/m2/season
    
    subsurf_JS (subsurf_JS == 0) = NaN;
    
    prop_JS = subsurf_JS ./ (surf_JS + subsurf_JS); % proportion of Jul-Sep production that happens subsurface
    
    figure(6);
    odvpal(50);
    suptitle('Post-bloom (July-September)');
    
    subplot(3,1,1);
    m_proj('stereographic','lat',90,'long',30,'radius',30)
    m_elev('contour',[-500 -500],'edgecolor','r'); hold on
    m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
    m_pcolor(xx,yy,subsurf_JS); shading flat; hold on; % top 4 boxs ~ top 20m
    m_coast('patch',[.7 .7 .7],'edgecolor','k');
    caxis([0 500]); cb = colorbar; ylabel(cb,'mgC/m2/season');
    title('Subsurface PP (20-60m)','FontWeight','bold')
    
    
    subplot(3,1,2);
    m_proj('stereographic','lat',90,'long',30,'radius',30)
    m_elev('contour',[-500 -500],'edgecolor','r'); hold on
    m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
    m_pcolor(xx,yy,surf_JS); shading flat; hold on; % top 4 boxs ~ top 20m
    m_coast('patch',[.7 .7 .7],'edgecolor','k');
    caxis([0 5000]); cb = colorbar; ylabel(cb,'mgC/m2/season');
    title('Surface PP (0-20m)','FontWeight','bold')
    
    subplot(3,1,3);
    m_proj('stereographic','lat',90,'long',30,'radius',30)
    m_elev('contour',[-500 -500],'edgecolor','r'); hold on
    m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
    m_pcolor(xx,yy,prop_JS); shading flat; hold on; % top 4 boxs ~ top 20m
    m_coast('patch',[.7 .7 .7],'edgecolor','k');
    caxis([0 0.1]); cb = colorbar; ylabel(cb,'dimensionless');
    title('Proportion of PP subsurface','FontWeight','bold')
    
    
    print -depsc -painters JS_PP_20_60_prop.eps
    
    
    
    % SUPPLEMENTARY FIGURE - Percentage of production that happens 20-60 m
    % for AJ, JS, OD
    
    
    % April to June
    AJ_files = 18; % number of files in Apr-Jun
    AJ = 91; % number of days in Apr-Jun
    
    subsurf_AJ = (nansum(nansum(modppC(:,:,5:9,19:36),4),3)) .* (sub_depth * AJ / sub_files * AJ_files); % mgC/m2/season
    surf_AJ = (nansum(nansum(modppC(:,:,1:4,19:36),4),3)) .* (sur_depth * AJ / sur_files * AJ_files); % mgC/m2/season
    prop_AJ = subsurf_AJ ./ (surf_AJ + subsurf_AJ);
    
    % Jul to Sep - part 3 of sup fig immediately above
    
    % Oct - Dec     
    OD_files = 19; % number of files in Jul-Sep
    OD = 92;
    
    subsurf_OD = (nansum(nansum(modppC(:,:,5:9,55:73),4),3)) .* (sub_depth * OD / sub_files * OD_files); % mgC/m2/season
    surf_OD = (nansum(nansum(modppC(:,:,1:4,55:73),4),3)) .* (sur_depth * OD / sur_files * OD_files); % mgC/m2/season
    prop_OD = subsurf_OD ./ (surf_OD + subsurf_OD);
    
    
    % Plotting: Subsurface PP proportion by season
    
    
    figure(7);
    odvpal(50);
    suptitle('Subsurface PP proportion by season');
    
    subplot(3,1,1);
    m_proj('stereographic','lat',90,'long',30,'radius',30)
    m_elev('contour',[-500 -500],'edgecolor','r'); hold on
    m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
    m_pcolor(xx,yy,prop_AJ); shading flat; hold on; % top 4 boxs ~ top 20m
    m_coast('patch',[.7 .7 .7],'edgecolor','k');
    caxis([0 0.15]); cb = colorbar; ylabel(cb,'dimless');
    title('Apr-Jun','FontWeight','bold')
    
    
    subplot(3,1,2);
    m_proj('stereographic','lat',90,'long',30,'radius',30)
    m_elev('contour',[-500 -500],'edgecolor','r'); hold on
    m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
    m_pcolor(xx,yy,prop_JS); shading flat; hold on; % top 4 boxs ~ top 20m
    m_coast('patch',[.7 .7 .7],'edgecolor','k');
    caxis([0 0.15]); cb = colorbar; ylabel(cb,'dimless');
    title('Jul-Sep','FontWeight','bold')
    
    subplot(3,1,3);
    m_proj('stereographic','lat',90,'long',30,'radius',30)
    m_elev('contour',[-500 -500],'edgecolor','r'); hold on
    m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
    m_pcolor(xx,yy,prop_OD); shading flat; hold on; % top 4 boxs ~ top 20m
    m_coast('patch',[.7 .7 .7],'edgecolor','k');
    caxis([0 0.15]); cb = colorbar; ylabel(cb,'dimensionless');
    title('Oct-Dec','FontWeight','bold')
    
    
    print -depsc -painters seasonal_subPP_prop.eps
    
    
    % SUPPLEMENTARY FIGURE - Percentage of production that happens 20-60 m
    % for May-Aug by month
    
%     lenM = [31 28 31 30 31 30 31 31 30 31 30 31];
%     filesM = [6 5 7 6 6 6 6 6 6 6 6 7];
    
    % May
    M_files = 6; % number of files in month
    M = 31; % number of days in month
    subsurf_M = (nansum(nansum(modppC(:,:,5:9,25:30),4),3)) .* (sub_depth * M / sub_files * M_files); % mgC/m2/month
    surf_M = (nansum(nansum(modppC(:,:,1:4,25:30),4),3)) .* (sur_depth * M / sur_files * M_files); % mgC/m2/month
    prop_M = subsurf_M ./ (surf_M + subsurf_M);
    
    % Jun
    Jn_files = 6; % number of files in month
    Jn = 30; % number of days in month
    subsurf_Jn = (nansum(nansum(modppC(:,:,5:9,31:36),4),3)) .* (sub_depth * Jn / sub_files * Jn_files); % mgC/m2/month
    surf_Jn = (nansum(nansum(modppC(:,:,1:4,31:36),4),3)) .* (sur_depth * Jn / sur_files * Jn_files); % mgC/m2/month
    prop_Jn = subsurf_Jn ./ (surf_Jn + subsurf_Jn);
    
    % Jul
    Jl_files = 6; % number of files in month
    Jl = 31; % number of days in month
    subsurf_Jl = (nansum(nansum(modppC(:,:,5:9,37:42),4),3)) .* (sub_depth * Jl / sub_files * Jl_files); % mgC/m2/month
    surf_Jl = (nansum(nansum(modppC(:,:,1:4,37:42),4),3)) .* (sur_depth * Jl / sur_files * Jl_files); % mgC/m2/month
    prop_Jl = subsurf_Jl ./ (surf_Jl + subsurf_Jl);
    
    % Aug
    Au_files = 6; % number of files in month
    Au = 31; % number of days in month
    subsurf_Au = (nansum(nansum(modppC(:,:,5:9,43:48),4),3)) .* (sub_depth * Au / sub_files * Au_files); % mgC/m2/month
    surf_Au = (nansum(nansum(modppC(:,:,1:4,43:48),4),3)) .* (sur_depth * Au / sur_files * Au_files); % mgC/m2/month
    prop_Au = subsurf_Au ./ (surf_Au + subsurf_Au);
    
    
    % correct for places where surf_Month = 0; places where surf_Month = 0 -> NaN then prop_Month -> NaN (not 1) where surf_Month = 0 (I think)
    for i = 1:size(surf_M,1)
        for j = 1:size(surf_M,2)
            
              if (surf_M(i,j) == 0) && (subsurf_M(i,j) == 0), prop_M(i,j) = NaN; % if surface and subsurface PP = 0, prop = NaN
              elseif (surf_M(i,j) == 0) && (subsurf_M(i,j) ~= 0), prop_M(i,j) = 1; % if surface = 0 and subsurface not 0, prop = 1
              end
              
              if (surf_Jn(i,j) == 0) && (subsurf_Jn(i,j) == 0), prop_Jn(i,j) = NaN; % if surface and subsurface PP = 0, prop = NaN
              elseif (surf_Jn(i,j) == 0) && (subsurf_Jn(i,j) ~= 0), prop_Jn(i,j) = 1; % if surface = 0 and subsurface not 0, prop = 1
              end
              
              if (surf_Jl(i,j) == 0) && (subsurf_Jl(i,j) == 0), prop_Jl(i,j) = NaN; % if surface and subsurface PP = 0, prop = NaN
              elseif (surf_Jl(i,j) == 0) && (subsurf_Jl(i,j) ~= 0), prop_Jl(i,j) = 1; % if surface = 0 and subsurface not 0, prop = 1
              end
              
              if (surf_Au(i,j) == 0) && (subsurf_Au(i,j) == 0), prop_Au(i,j) = NaN; % if surface and subsurface PP = 0, prop = NaN
              elseif (surf_Au(i,j) == 0) && (subsurf_Au(i,j) ~= 0), prop_Au(i,j) = 1; % if surface = 0 and subsurface not 0, prop = 1
              end
             
        end
    end
    
    
    % Plotting: Subsurface PP proportion by month
    
    
    figure(8);
    odvpal(50);
    suptitle('Subsurface PP proportion by month');
    
%     subplot(2,2,1);
    m_proj('stereographic','lat',90,'long',30,'radius',30)
    m_elev('contour',[-100 -100],'edgecolor','r'); hold on
    m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
    m_pcolor(xx,yy,prop_M); shading flat; hold on; % top 4 boxs ~ top 20m
    m_coast('patch',[.7 .7 .7],'edgecolor','k');
    caxis([0 1]); cb = colorbar; ylabel(cb,'dimless');
    title('May','FontWeight','bold')
    
    
    subplot(2,2,2);
    m_proj('stereographic','lat',90,'long',30,'radius',30)
    m_elev('contour',[-500 -500],'edgecolor','r'); hold on
    m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
    m_pcolor(xx,yy,prop_Jn); shading flat; hold on; % top 4 boxs ~ top 20m
    m_coast('patch',[.7 .7 .7],'edgecolor','k');
    caxis([0 0.2]); cb = colorbar; ylabel(cb,'dimless');
    title('Jun','FontWeight','bold')
    
    subplot(2,2,3);
    m_proj('stereographic','lat',90,'long',30,'radius',30)
    m_elev('contour',[-500 -500],'edgecolor','r'); hold on
    m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
    m_pcolor(xx,yy,prop_Jl); shading flat; hold on; % top 4 boxs ~ top 20m
    m_coast('patch',[.7 .7 .7],'edgecolor','k');
    caxis([0 0.2]); cb = colorbar; ylabel(cb,'dimless');
    title('Jul','FontWeight','bold')
    
    subplot(2,2,4);
    m_proj('stereographic','lat',90,'long',30,'radius',30)
    m_elev('contour',[-500 -500],'edgecolor','r'); hold on
    m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
    m_pcolor(xx,yy,prop_Au); shading flat; hold on; % top 4 boxs ~ top 20m
    m_coast('patch',[.7 .7 .7],'edgecolor','k');
    caxis([0 0.2]); cb = colorbar; ylabel(cb,'dimless');
    title('Aug','FontWeight','bold')
    
    
    print -depsc -painters monthly_subPP_prop.eps