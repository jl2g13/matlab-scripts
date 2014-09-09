% Final figure 1

% Created by jl2g13 (04/09/13) from oringinal scripts 'subsurface
% production plot 1' and 'subsurface production plot 2'.

% Plot has 4 panels as follows
% a. PP 0-20m.
% b. PP 20-60m
% c. PP proportion 20-60m
% d. bar plot of DIN/PP



% Read in production fields
% -------------------------


% (1) Reading in 5-day model variable - PP   (need 5-day files for XCORR)

% Only considering 2005
yr = 2005;

lattotal = double(1021);
latstart = double(1);
lontotal = double(1442);
depth = double(10); % 14 = read in top 60 m
pos = 1;
months = 12;


% Preallocating arrays for model output to be read into
modndpp = NaN(lontotal,lattotal,depth,months);
moddpp = NaN(lontotal,lattotal,depth,months);
modN = NaN(lontotal,lattotal,depth,months);

    
  for m = 1:12
    if m < 10
        dnomP = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/ORCA025-N201_%dm0%dP.nc', yr, yr, m);
        dnomD = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/ORCA025-N201_%dm0%dD.nc', yr, yr, m);
        dnomT = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/ORCA025-N201_%dm0%dT.nc', yr, yr, m);
    else
        dnomP = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/ORCA025-N201_%dm%dP.nc', yr, yr, m);
        dnomD = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/ORCA025-N201_%dm%dD.nc', yr, yr, m);
        dnomT = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/ORCA025-N201_%dm%dT.nc', yr, yr, m);
    end


    fnamesD = dir(dnomD);
    fnamesP = dir(dnomP);
    fnamesT = dir(dnomT);
    
    fprintf('Reading month %d\n', m);    
    

    
    % (1) Primary Production

    fnameD = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/%s', yr, fnamesD.name);
    fprintf('Reading %s\n', fnameD);

    % (i) For non-diatom production

    t1 = ncread(fnameD,'PRN3',[1 latstart 1 1], [lontotal lattotal depth 1]);
    t1 (t1 == 0) = NaN; % removing zero values (land)

    modndpp(:,:,:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)

    % (ii) For diatom production
    t1 = ncread(fnameD,'PRD3',[1 latstart 1 1], [lontotal lattotal depth 1]);
    t1 (t1 == 0) = NaN; % removing zero values (land)

    moddpp(:,:,:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)
    
    
    
    % (2) Dissolved Inorganic Nitrogen   
    
    fnameP = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/%s', yr, fnamesP.name);
    fprintf('Reading %s\n', fnameP);

    % (a) Depth-resolved

    % Reading in nitrogen model field - DIN
    t1 = ncread(fnameP,'DIN',[1 latstart 1 1], [lontotal lattotal depth 1]);
    t1 (t1 == 0) = NaN; % removing zero values (land) - DIN might actually be zero
    modN(:,:,:,pos) = t1;
    
    
    % (3) For Mixed Layer Depth
    
    fnameT = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/%s', yr, fnamesT.name);
    fprintf('Reading %s\n', fnameT);

    % Reading in mixed layer depth model field
    t1 = ncread(fnameT,'somxl010',[1 latstart 1], [lontotal lattotal 1]); % mixing layer depth -> Tfile,somixhgt (not so meaningful for monthly averages)
    t1 = squeeze(t1);
    t1 (t1 == 0) = NaN; % removing zero values (land)

    modmld(:,:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)

    
    pos = pos + 1;

end

modndpp (modndpp == 9.969209968386869e+36) = NaN; % what on earth is this flag? - netcdf fill value, probably occurs when read in problem
moddpp (moddpp == 9.969209968386869e+36) = NaN;
modN (modN == 9.969209968386869e+36) = NaN; 
modmld (modmld == 9.969209968386869e+36) = NaN; 

modppC = (modndpp + moddpp) * 12.011 * 6.625; % in mgC

clearvars modndpp moddpp t1

% remove anywhere water depth < 60 m -> NaN: take all depths down to filter depth in column; sum over
% them, if any of them are NaN (rock) then t1->0 at this point -> NaN.
depthfilter = 4; % nemo_dep(4) ~ 20m

t1 = modN(:,:,1:depthfilter,1) * 0; % only want filter at top 20 m -> not over whole grid. nb. n=4, is z~60m. 
t2 = sum(isfinite(t1),3); % chl = finite in water, = NaN below seafloor
kbathy = t2; % model water depth in terms of k model levels
kbathy (kbathy < depthfilter) = NaN; % from here just times whatever model field you want to filter by kbathy
kbathy (kbathy == depthfilter) = 1; % when sums over finite values each x-y water column gets value 1-10 if > 60m deep (would be 10 if didn't set == 0 -> NaN at read in above)
% NOTE IF A DIN VALUE IN THE WATER COLUMN -> 0 THEN THIS WATER COLUMN WILL
% GET FILTERED OUT - DOESN'T HAPPEN 1ST-5TH JAN.

% -----------------------------
% check kbathy filter - should be NaN in water depth < 60, 1 everywhere
% else

ft='/noc/altix2/scratch/omfman/ORCA025-N201/means/1990/ORCA025-N201_1990m09P.nc';
xx=ncread(ft,'nav_lon', [1 1], [lontotal lattotal]);
yy=ncread(ft,'nav_lat', [1 1], [lontotal lattotal]);

figure;
odvpal(50);

m_proj('stereographic','lat',90,'long',30,'radius',30)
m_elev('contour',[-500 -500],'edgecolor','r'); hold on
m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
m_pcolor(xx,yy,kbathy); shading flat; hold on;
m_coast('patch',[.7 .7 .7],'edgecolor','k');
caxis([0 2]); cb = colorbar; ylabel(cb,'');
% ----------------------------


clearvars t1 t2 % prn prd % clear up memory




% Creating bar plot variables
% ---------------------------


% A. Creating MLD filter
% --------------------
% Filtering out areas where mixed layer depth > 100 m (i.e. Atlantic inflow)
mldfilt = max(modmld,[],3); % select maximum value along dim = 3 (month)

mldfilt (mldfilt > 100) = NaN;
mldfilter = ~isnan(mldfilt);



% % --------------
% % check mldfilter filter
% 
% ft='/noc/altix2/scratch/omfman/ORCA025-N201/means/1990/ORCA025-N201_1990m09P.nc';
% xx=ncread(ft,'nav_lon', [1 1], [lontotal lattotal]);
% yy=ncread(ft,'nav_lat', [1 1], [lontotal lattotal]);
% 
% figure;
% odvpal(50);
% 
% m_proj('stereographic','lat',90,'long',30,'radius',30)
% m_elev('contour',[-500 -500],'edgecolor','r'); hold on
% m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
% m_pcolor(xx,yy,mldfilter); shading flat; hold on;
% m_coast('patch',[.7 .7 .7],'edgecolor','k');
% caxis([0 2]); cb = colorbar; ylabel(cb,'');
% % --------------

% filter dim must match modppC and modN to multiply
filtrep = repmat(mldfilter,[1 1 10 12]);

filtPP = modppC .* filtrep;
filtN = modN .* filtrep;




% B: whole Arctic
% ---------------

% integrate over depth
modppC_20 = squeeze(nansum(modppC(:,:,1:4,:),3) ./ 4); % integrating over 0-20 m depth
modN_20 = nansum(modN(:,:,1:4,:),3) ./ 4;

modppC_60 = nansum(modppC(:,:,5:9,:),3) ./ 5; % integrating over 20-60 m depth
modN_60 = nansum(modN(:,:,5:9,:),3) ./ 5;



% intergrate over area (x,y)
modppC_20 = squeeze(nansum(nansum(modppC_20(:,831:end,:),2),1)) ./ ((1021-830) * 1442); % integrating over lat and lon of Arctic (lat > 831 ~ 66° N = Arctic circ ('~' because grid tripolar)) (mgC/m3/d)
modppC_60 = squeeze(nansum(nansum(modppC_60(:,831:end,:),2),1)) ./ ((1021-830) * 1442); 

modN_20 = squeeze(nansum(nansum(modN_20(:,831:end,:),2),1)) ./ ((1021-830) * 1442); % integrating over lat and lon of Arctic (lat > 831 ~ 66° N = Arctic circ ('~' because grid tripolar)) (mmolN/m3)
modN_60 = squeeze(nansum(nansum(modN_60(:,831:end,:),2),1)) ./ ((1021-830) * 1442);





% C: Arctic with MLD filter applied
% ---------------

% integrate over depth
filtppC_20 = squeeze(nansum(filtPP(:,:,1:4,:),3) ./ 4); % integrating over 0-20 m depth
filtN_20 = nansum(filtN(:,:,1:4,:),3) ./ 4;

filtppC_60 = nansum(filtPP(:,:,5:9,:),3) ./ 5; % integrating over 20-60 m depth
filtN_60 = nansum(filtN(:,:,5:9,:),3) ./ 5;


% intergrate over area (x,y)
filtppC_20 = squeeze(nansum(nansum(filtppC_20(:,831:end,:),2),1)) ./ ((1021-830) * 1442); % integrating over lat and lon of Arctic (lat > 831 ~ 66° N = Arctic circ ('~' because grid tripolar)) (mgC/m3/d)
filtppC_60 = squeeze(nansum(nansum(filtppC_60(:,831:end,:),2),1)) ./ ((1021-830) * 1442); 

filtN_20 = squeeze(nansum(nansum(filtN_20(:,831:end,:),2),1)) ./ ((1021-830) * 1442); % integrating over lat and lon of Arctic (lat > 831 ~ 66° N = Arctic circ ('~' because grid tripolar)) (mmolN/m3)
filtN_60 = squeeze(nansum(nansum(filtN_60(:,831:end,:),2),1)) ./ ((1021-830) * 1442);








% Creating PP plot fields
% -----------------------

modppC_an = (nansum(modppC,4)) .* (365/months); % Annual PP (integrating over all months in year). days in yr/no files. mgC/m3/yr
    
    
% clearvars modppC

plot1a = squeeze(nansum(modppC_an(:,:,1:4),3) * (20/4) / 1000); % surface PP = integrate top 4 boxs ~ top 20m. depth-interval/noboxes. gC/m2/yr
plot1b_t = squeeze(nansum(modppC_an(:,:,5:10),3) * (40/6) / 1000); % deep PP = integrate boxes 5-10 ~ 20-60 m depth. gC/m2/yr

prop_t = plot1b_t ./ (plot1a + plot1b_t);

plot1b = plot1b_t .* kbathy; % as above but relevant depth = 60 m
prop = prop_t .* kbathy;


% ======
% PLOT 1
% ======
% M_map plotting preamble

% Setting up M_Map model
ft='/noc/altix2/scratch/omfman/ORCA025-N201/means/1990/ORCA025-N201_1990m09P.nc';
xx=ncread(ft,'nav_lon', [1 1], [lontotal lattotal]);
yy=ncread(ft,'nav_lat', [1 1], [lontotal lattotal]);

% removing spurious spikes in N Atl when plotting (beware true zeros, quick fix)
land = xx(1,1);
xx(xx == land) = NaN; xx(xx == 0) = NaN;
land = yy(1,1);
yy(yy == land) = NaN; yy(yy == 0) = NaN;



% PP and DIN bar plot

% plotting preamble
NumStacksPerGroup = 2; % N and PP
NumGroupsPerAxis = 12; % months
NumStackElements = 2; % depth windows
% groupLabels = { 'Jan'; 'Feb'; 'Mar'; 'Apr'; 'May'; 'Jun'; 'Jul'; 'Aug'; 'Sep'; 'Oct'; 'Nov'; 'Dec' }; % Months
groupLabels = { 1:12 }; % Months


% manually set figure layout    
figure;
odvpal(50);

% PP 0-20m
positionVector1 = [0.15, 0.6, 0.3, 0.4]; % first two [left bottom width height] (values 0-1) . left-bottom positions bottom left corner of figure, width-height to set size
subplot('Position',positionVector1)

m_proj('stereographic','lat',90,'long',30,'radius',30)
m_elev('contour',[-500 -500],'edgecolor','r'); hold on
m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
m_pcolor(xx,yy,plot1a); shading flat; hold on; % top 4 boxs ~ top 20m
m_coast('patch',[.7 .7 .7],'edgecolor','k');
caxis([0 100]); cb = colorbar; ylabel(cb,'gC/m2/yr'); % 0-100 gC/m2/yr
%     title('Annual PP 0-20 m','FontWeight','bold')
title('a.','position',[-0.5 0.6],'FontWeight','bold')


% PP 20-60m
positionVector2 = [0.55, 0.6, 0.3, 0.4];
subplot('Position',positionVector2)

m_proj('stereographic','lat',90,'long',30,'radius',30)
m_elev('contour',[-500 -500],'edgecolor','r'); hold on
m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
m_pcolor(xx,yy,plot1b); shading flat; hold on; % boxes 5--9 ~ 20-60 m depth 
m_coast('patch',[.7 .7 .7],'edgecolor','k'); 
caxis([0 30]); cb = colorbar; ylabel(cb,'gC/m2/yr');
%     title('Annual PP 20-60 m','FontWeight','bold')
title('b.','position',[-0.5 0.6],'FontWeight','bold')     


% PP prop 20-60m
positionVector3 = [0.15, 0.1, 0.35, 0.5];
subplot('Position',positionVector3)

m_proj('stereographic','lat',90,'long',30,'radius',30)
m_elev('contour',[-500 -500],'edgecolor','r'); hold on
m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
m_pcolor(xx,yy,prop); shading flat; hold on; % boxes 5--9 ~ 20-60 m depth 
m_coast('patch',[.7 .7 .7],'edgecolor','k'); 
caxis([0 0.6]); cb = colorbar; ylabel(cb,'dimensionless');
%     title('Annual proportion 20-60 m','FontWeight','bold')
title('c.','position',[-0.5 0.6],'FontWeight','bold')


% bar plot din/pp
stackData = reshape([modN_60 modppC_60 modN_20 modppC_20], 12,2,2); % dim ordering needs to be (NumGroupsPerAxis,NumStacksPerGroup,NumStackElements) = [v1_s1 v2_s1 v1_s2 v2_s2 v1_s3 v2_s3 etc]

positionVector4 = [0.6, 0.2, 0.3, 0.3];
subplot('Position',positionVector4)

plotBarStackGroups(stackData, groupLabels);
set(gca,'FontSize',9) 
set(gcf,'Position',[100 100 720 650]) 
grid off 
set(gca,'Layer','top') % put grid lines on top of stacks
% title('Surface and subsurface primary production','FontWeight','bold'); % Add title and axis labels
xlabel('Month','fontsize',6);
ylabel('PP (mgC/m3/d) or DIN (mmolN/m3)','fontsize',6);
% ylabel({'Primary production (mgC/m3/d)','DIN (mmolN/m3)'},'fontsize',10); % cell array to split title over two lines
title('d.','position',[-2.5 5.5],'FontWeight','bold')
figd = (gcf); % inset figure function required handles of two figures to plot together







% Aux figure - bar plot with filtered bar plot as inset figure



% inset bar plot of filtered din/pp
% stackData = reshape([modN_60 modppC_60 modN_20 modppC_20], 12,2,2); % dim ordering needs to be (NumGroupsPerAxis,NumStacksPerGroup,NumStackElements) = [v1_s1 v2_s1 v1_s2 v2_s2 v1_s3 v2_s3 etc]
% 
% figure;
% plotBarStackGroups(stackData, groupLabels);
% set(gca,'FontSize',9) 
% set(gcf,'Position',[100 100 720 650]) 
% grid off 
% set(gca,'Layer','top') % put grid lines on top of stacks
% % title('Surface and subsurface primary production','FontWeight','bold'); % Add title and axis labels
% xlabel('Month','fontsize',6);
% ylabel('PP (mgC/m3/d) or DIN (mmolN/m3)','fontsize',6);
% % ylabel({'Primary production (mgC/m3/d)','DIN (mmolN/m3)'},'fontsize',10); % cell array to split title over two lines
% title('d.','position',[-2.5 5.5],'FontWeight','bold')
% figd = (gcf); % inset figure function required handles of two figures to plot together
% 
% 
% 
% % preamble
% stackData = reshape([filtN_60 filtppC_60 filtN_20 filtppC_20], 12,2,2); % other preamble as initialised above
% 
% figure;
% plotBarStackGroups(stackData, groupLabels);
% set(gca,'FontSize',9) 
% set(gcf,'Position',[100 100 720 650]) 
% grid off
% set(gca,'Layer','top') % put grid lines on top of stacks
% % title('Surface and subsurface primary production','FontWeight','bold'); % Add title and axis labels
% % ylabel({'Primary production (mgC/m3/d)','DIN (mmolN/m3)'},'fontsize',10); % cell array to split title over two lines
% title('e.','position',[-2.5 5.5],'FontWeight','bold')
% fige = (gcf); % inset figure function required handles of two figures to plot together
% 
% 
% [h_m h_i]=inset(figd,fige);
% 
% 
% 
% print -depsc -painters monthly_pp_DIN.eps_mld_filtered.eps






