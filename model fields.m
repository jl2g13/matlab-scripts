% Plotting model fields

% Reading in model variables

pos = 1; % Initialsing index for where to put data in array


for y = 1:1:1 % Looping over one year -> all years (full possible interval is 1990-2006)
    yr = 2004 + y;
    
    %     % Prescribing lattitude extent of interest
    lattotal = double(1021); % grid axis 1021 long
    %     latstart = double(825);  % corresponds to 65°N (coord2index function)
    latstart = double(1);
    %     latrange = double(lattotal - latstart + 1);  % +1 because range must include both bounds of range
    
    % Prescribing longditude extent of interest
    lontotal = double(1442); % Pan-Arctic so all longditudes taken
    
    % Prescribing depth to read model variable in over for 3D fields.
    depth = 1;  % 1~3 m, 14~100 m.
    depth_res = 14;
    
    
    % Getting list of all files for this year
    
    % (1) get a list of P (Passive Tracers) files for this year
    dnomP = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/*d05P.nc', yr);
    fnamesP = dir (dnomP);
    
    % (2) get a list of D (Diagnostic) files for this year
    dnomD = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/*d05D.nc', yr);
    fnamesD = dir (dnomD);
    
    % (3) get a list of T (SST, seaice etc.) files for this year
    dnomT = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/*d05T.nc', yr);
    fnamesT = dir (dnomT);
    
    % (4) get a list of I (Ice) files for this year
    dnomI = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/*d05I.nc', yr);
    fnamesI = dir (dnomI);
    
    
    
    % Finding number of  files to loop over in each yearly folder (73 for 5 day means; n.b. needs to be same for all file types)
    [num, x] = size(fnamesP);
    
    
    % Preallocating arrays for model output to be read into
    
      modsst = zeros(lontotal,lattotal,num);
      modsurfpar = zeros(lontotal,lattotal,num);
      modndpp = zeros(lontotal,lattotal,num);
      moddpp = zeros(lontotal,lattotal,num);
      modice = zeros(lontotal,lattotal,num);
%       modsurfchn = zeros(lontotal,lattotal,depth,num);
%       modsurfchd = zeros(lontotal,lattotal,depth,num);
      modchn = zeros(lontotal,lattotal,depth_res,num);
      modchd = zeros(lontotal,lattotal,depth_res,num);
      modsurfndpp = zeros(lontotal,lattotal,depth,num);
      modsurfdpp = zeros(lontotal,lattotal,depth,num);
      modmld = zeros(lontotal,lattotal,num);
      modicethick = zeros(lontotal,lattotal,num);
      modsurfN = zeros(lontotal,lattotal,depth,num);
      modsurfFe = zeros(lontotal,lattotal,depth,num);
      modsurfSi = zeros(lontotal,lattotal,depth,num);

    
    % Initiating extraction loop
    
    for f = 1:1:num % looping from file 1 to num (73) - for each year

        
        
        % (1) For Chlorophyll - saved modchl variable is sum of top 13 boxes (~95 m)
        
        % Extracting model surface chlorophyll
        fnameP = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/%s', yr, fnamesP(f).name);
        fprintf('Reading %s\n', fnameP);
        
%         % (a) Surface model
%         
%         % (i) For non-diatom chlorophyll
%         
%         t1 = ncread(fnameP,'CHN',[1 latstart 1 1], [lontotal lattotal depth 1]); % can set this to read in just N of Arctic Circle
%         t1 = squeeze(t1);
%         t1 (t1 == 0) = NaN; % removing zero values (land)
%         
%         modsurfchn(:,:,:,pos) = t1; % lat, lon, no of 5-day means (73/yr)
%         
%         % (ii) For diatom chlorophyll
%         
%         t1 = ncread(fnameP,'CHD',[1 latstart 1 1], [lontotal lattotal depth 1]);
%         t1 = squeeze(t1);
%         t1 (t1 == 0) = NaN; % removing zero values (land)
%         
%         modsurfchd(:,:,:,pos) = t1;
%         
%         % (b) Depth integrated
%         
%         % (i) For non-diatom chlorophyll
%         
%         t1 = ncread(fnameP,'CHN',[1 latstart 1 1], [lontotal lattotal depth_res 1]); % can set this to read in just N of Arctic Circle
%         t1 = squeeze(t1);
%         t1 (t1 == 0) = NaN; % removing zero values (land)
%         
%         modchn(:,:,:,pos) = t1; % lat, lon, no of 5-day means (73/yr)
%         
%         % (ii) For diatom chlorophyll
%         
%         t1 = ncread(fnameP,'CHD',[1 latstart 1 1], [lontotal lattotal depth_res 1]);
%         t1 = squeeze(t1);
%         t1 (t1 == 0) = NaN; % removing zero values (land)
%         
%         modchd(:,:,:,pos) = t1;
%         
%         
%         
%         
%       % (2) For Photosynthetically Active Radiation
%         
        % Extracting model surface PAR
        fnameD = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/%s', yr, fnamesD(f).name);
        fprintf('Reading %s\n', fnameD);
%         
%         % Reading in surface PAR model field
%         t1 = ncread(fnameD,'MED_QSR',[1 latstart 1], [lontotal lattotal 1]); % MED_QSR ('top level radiation) is a 2D model field, has no depth dimension
%         t1 = squeeze(t1);
%         t1 (t1 == 0) = NaN; % removing zero values (land)
%         
%         totalshortwave(:,:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)
%         
%         
%         
%         % (3) Primary Production 
% 
%         % (a) Model depth-integrated
%                 
%         % (i) For non-diatom production
%         
%         t1 = ncread(fnameD,'PRN',[1 latstart 1], [lontotal lattotal 1]); % PRN is depth-integrated so a 2D model field
%         t1 = squeeze(t1);
%         t1 (t1 == 0) = NaN; % removing zero values (land)
%         
%         modndpp(:,:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)
%         
%         % (ii) For diatom production
%         t1 = ncread(fnameD,'PRD',[1 latstart 1], [lontotal lattotal 1]); % PRD is depth-integrated so a 2D model field
%         t1 = squeeze(t1);
%         t1 (t1 == 0) = NaN; % removing zero values (land)
%         
%         moddpp(:,:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)
% 

        % (b) Depth-resolved - modintppC is PP summed over top 13 boxes

        % (i) For non-diatom production

        t1 = ncread(fnameD,'PRN3',[1 latstart 1 1], [lontotal lattotal depth 1]);
        t1 = squeeze(t1);
        t1 (t1 == 0) = NaN; % removing zero values (land)
         
        modsurfndpp(:,:,:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)

        % (ii) For diatom production
        t1 = ncread(fnameD,'PRD3',[1 latstart 1 1], [lontotal lattotal depth 1]); % PRD is depth-integrated so a 2D model field
        t1 = squeeze(t1);
        t1 (t1 == 0) = NaN; % removing zero values (land)
    
        modsurfdpp(:,:,:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)
%         
%         
%         
%         
%         % (4) For Sea Surface Temperature
%         
%         % Extracting model sea surface temperature
%         fnameT = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/%s', yr, fnamesT(f).name);
%         fprintf('Reading %s\n', fnameT);
%         
%         % Reading in sea surface temperature model field
%         t1 = ncread(fnameT,'sosstsst',[1 latstart 1], [lontotal lattotal 1]); % sosstsst is a 2D model field, has no depth dimension
%         t1 = squeeze(t1);
%         t1 (t1 == 0) = NaN; % removing zero values (land)
%         
%         modsst(:,:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)
%         
%         
%         
%         % (5) For Sea Ice Cover
%         
%         % Reading in sea ice cover model field
%         t1 = ncread(fnameT,'soicecov',[1 latstart 1], [lontotal lattotal 1]); % sea ice cover is a 2D model field, has no depth dimension
%         t1 = squeeze(t1);
%         % t1 (t1 == 0) = NaN; % removing zero values (land)
%         
%         modice(:,:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)
%         
%         
%         
%         % (6) For Sea Ice Thickness
%         
%         % Extracting model sea ice thickness
%         fnameI = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/%s', yr, fnamesI(f).name);
%         fprintf('Reading %s\n', fnameI);
%         
%         % Reading in sea ice thickness model field
%         t1 = ncread(fnameI,'iicethic',[1 latstart 1], [lontotal lattotal 1]); % sea ice thickness is a 2D model field, has no depth dimension
%         t1 = squeeze(t1);
%         % t1 (t1 == 0) = NaN; % removing zero values (land)
%         
%         modicethick(:,:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)
%         
% 
%         
%         % (7) For Mixed Layer Depth
% 
%         % Reading in mixed layer depth model field
%         t1 = ncread(fnameT,'somxl010',[1 latstart 1], [lontotal lattotal 1]); % mixing layer depth -> Tfile,somixhgt (not so meaningful for monthly averages)
%         t1 = squeeze(t1);
%         t1 (t1 == 0) = NaN; % removing zero values (land)
%         
%         modmld(:,:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)
%         
%         
%         
%         % (8) For Nitrogen        
%         
%         % (a) Depth-resolved - saved modN variable is sum of top 13 boxes (~95 m), modsurfN is for top box (3m)
%         
%         % Reading in surface nitrogen model field - DIN
%         t1 = ncread(fnameP,'DIN',[1 latstart 1 1], [lontotal lattotal depth 1]);
%         t1 = squeeze(t1);
%         % t1 (t1 == 0) = NaN; % removing zero values (land) - DIN might actually be zero
%         
%         modsurfN(:,:,:,pos) = t1;        
%         
%         
%         
%         % (9) For Iron
%         
%         % Reading in surface dissolved iron model field - just the top 3 m
%         t1 = ncread(fnameP,'FER',[1 latstart 1 1], [lontotal lattotal depth 1]);
%         t1 = squeeze(t1);
%         % t1 (t1 == 0) = NaN; % removing zero values (land) - FER might actually be zero
%         
%         modsurfFe(:,:,:,pos) = t1;   
%         
%         
%         
%          % (10) For Silicate
%          
%         % Reading in surface dissolved iron model field - just the top 3 m
%         t1 = ncread(fnameP,'SIL',[1 latstart 1 1], [lontotal lattotal depth 1]);
%         t1 = squeeze(t1);
%         % t1 (t1 == 0) = NaN; % removing zero values (land) - FER might actually be zero
%         
%         modsurfSi(:,:,:,pos) = t1;
%         
                 
        pos = pos + 1; % resetting position counter for next iteration
         
    end
end

modsurfpar = totalshortwave .* 0.43; % Fudge conversion from total shortwave radiation to PAR (used previously by Andrew)
modchl = nansum((modchn + modchd),3); % integrating over depth
modsurfchl = squeeze(modchl(:,:,1,:));
% modsurfchl = modsurfchd + modsurfchn; % adding diatom and non-diatom chl components
modpp = modndpp + moddpp; % adding diatom and non-diatom primary production components
modppC = modpp * 12.011 * 6.625; % conversion from mmolN (model) to mgC (satellites) is * 14 (g/mol) * 6.625 (assuming Redfieldian)

modsurfpp = modsurfndpp + modsurfdpp;
modsurfppC = modsurfpp * 12.011 * 6.625; % surface field of PP

% Emptying memory
clearvars totalshortwave modsurfchd modsurfchn modndpp moddpp modpp modsurfndpp modsurfdpp modsurfpp modchn modchd


% Plotting preamble

ft='/noc/altix2/scratch/omfman/ORCA025-N201/means/1990/ORCA025-N201_1990m09P.nc';
xx=ncread(ft,'nav_lon', [1 1], [lontotal lattotal]);
yy=ncread(ft,'nav_lat', [1 1], [lontotal lattotal]);


% removing spurious spikes in N Atl when plotting (beware true zeros, quick fix)
land = xx(1,1);
xx(xx == land) = NaN; xx(xx == 0) = NaN;
land = yy(1,1);
yy(yy == land) = NaN; yy(yy == 0) = NaN;





% Plot 1 UML and surface nitrogen


    t1 = modmld(:,:,37:42); % Jul t = 37-42, Sep t = 49-54  
    t2 = nansum(t1,3) ./ 6; % division is by number of files summed (37-42) because files are in units mgChl/m2/day - adding 6 files of daily rate together won't give you the daily rate (averages are sum(individuals)/number(individuals)!)

    t1 = modsurfN(:,:,37:42);
    t3 = nansum(t1,3) ./ 6;


    subplot(2,1,1)
        odvpal(10);
        m_proj('stereographic','lat',90,'long',30,'radius',30)
        m_elev('contour',[-500 -500],'edgecolor','r'); hold on
        m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
        m_pcolor(xx,yy,t2); shading flat; hold on;
        m_coast('patch',[.7 .7 .7],'edgecolor','k');
        caxis([0 50]); cb = colorbar; ylabel(cb, 'MLD (m)');    % Sep plot - caxis = 0 20
        title ('Model Mixed Layer Depth, Jul 2005','FontWeight','bold');
        

    subplot(2,1,2)
        odvpal(10);
        m_proj('stereographic','lat',90,'long',30,'radius',30)
        m_elev('contour',[-500 -500],'edgecolor','r'); hold on
        m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
        m_pcolor(xx,yy,t3); shading flat; hold on;
        m_coast('patch',[.7 .7 .7],'edgecolor','k');
        caxis([0 5]); cb = colorbar; ylabel(cb, 'DIN (mmol-N/m3)');   % is it /m3 - taken a 3m box of it with this average value, not summed over the 3m N inventory
        title ('Model Surface DIN, Jul 2005','FontWeight','bold');



% Plot 2 Ice cover and thickness

    
    t1 = modice(:,:,37:42);
    t4 = nansum(t1,3) ./ 6;
    
    t1 = modicethick(:,:,37:42);
    t5 = nansum(t1,3) ./ 6;
    
    
    subplot(2,1,1)    
        odvpal(10);
        m_proj('stereographic','lat',90,'long',30,'radius',30)
        m_elev('contour',[-500 -500],'edgecolor','r'); hold on
        m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
        m_pcolor(xx,yy,t4); shading flat; hold on;
        m_coast('patch',[.7 .7 .7],'edgecolor','k');
        caxis([0 1]); cb = colorbar; ylabel(cb,'Fractional ice cover')
        title ('Model ice cover, Jul 2005','Fontweight','bold')
        
        
    subplot(2,1,2)    
        odvpal(10);
        m_proj('stereographic','lat',90,'long',30,'radius',30)
        m_elev('contour',[-500 -500],'edgecolor','r'); hold on
        m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
        m_pcolor(xx,yy,t5); shading flat; hold on;
        m_coast('patch',[.7 .7 .7],'edgecolor','k');
        caxis([0 5]); cb = colorbar; ylabel(cb,'Ice thickness (m)')
        title ('Model ice thickness, Jul 2005','Fontweight','bold')  
        
        
        
% Plot 3 Depth-integrated (100 m) and surface chlorophyll 


    t1 = modchl(:,:,37:42);
    t6 = nansum(t1,3) ./ 6;
    
    t1 = modsurfchl(:,:,37:42);
    t7 = nansum(t1,3) ./ 6;
    
    t1 = modice(:,:,37:42); % setting up ice array to add as contour (vars used as above in modice pcolor plots)
    t4 = nansum(t1,3) ./ 6; 
    
    
    subplot(2,1,1)    
        odvpal(10);
        m_proj('stereographic','lat',90,'long',30,'radius',30)
        m_elev('contour',[-500 -500],'edgecolor','r'); hold on
        m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
        m_pcolor(xx,yy,t6); shading flat; hold on;
        [cs,h]=m_contour(xx,yy,t4,[0.1],'edgecolor','k','LineWidth',1); hold on % 10%  ice cover contour 
        m_coast('patch',[.7 .7 .7],'edgecolor','k');
        caxis([0 10]); cb = colorbar; ylabel(cb,'Chl (mgChl/m2/day)')   % chl is in /m2 because you integrated over depth to produce variable modchl
        title ('Model depth-integrated chlorophyll (95m), Jul 2005','Fontweight','bold')
        
        
    subplot(2,1,2)    
        odvpal(10);
        m_proj('stereographic','lat',90,'long',30,'radius',30)
        m_elev('contour',[-500 -500],'edgecolor','r'); hold on
        m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
        m_pcolor(xx,yy,t7); shading flat; hold on;
        [cs,h]=m_contour(xx,yy,t4,[0.1],'edgecolor','k','LineWidth',1); hold on % 10%  ice cover contour 
        m_coast('patch',[.7 .7 .7],'edgecolor','k');
        caxis([0 2]); cb = colorbar; ylabel(cb,'Chl (mgChl/m3/day)')
        title ('Model surface chlorophyll (3m), Jul 2005','Fontweight','bold')
             
        
        
% Plot 4 Depth-integrated (100 m) and surface primary production
        

    t1 = modppC(:,:,37:42);
    t8 = nansum(t1,3) ./ 6;
    
    t1 = modsurfppC(:,:,37:42);
    t9 = nansum(t1,3) ./ 6;
    
    t1 = modice(:,:,37:42); % setting up ice array to add as contour (vars used as above in modice pcolor plots)
    t4 = nansum(t1,3) ./ 6; 
    
    
    subplot(2,1,1)    
        odvpal(10);
        m_proj('stereographic','lat',90,'long',30,'radius',30)
        m_elev('contour',[-500 -500],'edgecolor','r'); hold on
        m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
        m_pcolor(xx,yy,t8); shading flat; hold on;
        [cs,h]=m_contour(xx,yy,t4,[0.1],'edgecolor','k','LineWidth',1); hold on % 10%  ice cover contour 
        m_coast('patch',[.7 .7 .7],'edgecolor','k');
        caxis([0 1000]); cb = colorbar; ylabel(cb,'PP (mgC/m2/day)')   % chl is in /m2 because you integrated over depth to produce variable modchl
        title ('Model depth-integrated primary production (95m), Jul 2005','Fontweight','bold')
        
        
    subplot(2,1,2)    
        odvpal(10);
        m_proj('stereographic','lat',90,'long',30,'radius',30)
        m_elev('contour',[-500 -500],'edgecolor','r'); hold on
        m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
        m_pcolor(xx,yy,t9); shading flat; hold on;
        [cs,h]=m_contour(xx,yy,t4,[0.1],'edgecolor','k','LineWidth',1); hold on % 10%  ice cover contour 
        m_coast('patch',[.7 .7 .7],'edgecolor','k');
        caxis([0 100]); cb = colorbar; ylabel(cb,'PP (mgC/m3/day)')
        title ('Model surface primary production (3m), Jul 2005','Fontweight','bold')
                
        
% Plot 5 Checks: what fraction of PP occurs at > 95 m and how much DIN is there in upper 95 m

    t1 = modppC(:,:,37:42);
    t2 = nansum(t1,3) ./ 6;
    
    t1 = modintppC(:,:,37:42);
    t3 = nansum(t1,3) ./ 6;
    
    t4 = (t2 - t3) ./ t2;    % fraction of PP at > 95 m
    
    % ----
    
    t5 = modN(:,:,37:42);
    t6 = nansum(t1,3) ./ 6;
    
    
    subplot(2,1,1) 
        odvpal(10);
        m_proj('stereographic','lat',90,'long',30,'radius',30)
        m_elev('contour',[-500 -500],'edgecolor','r'); hold on 
        m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
        m_pcolor(xx,yy,t4); shading flat; hold on;
        m_coast('patch',[.7 .7 .7],'edgecolor','k');
        caxis([0 1]); cb = colorbar; ylabel(cb,'fraction of total PP')
        title('Fraction of Model PP below 95 m, Jul 2005','FontWeight','bold')
    
    % --> in July basically all of the PP is at >100 m
    

    subplot(2,1,2)    
        odvpal(10);
        m_proj('stereographic','lat',90,'long',30,'radius',30)
        m_elev('contour',[-500 -500],'edgecolor','r'); hold on
        m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
        m_pcolor(xx,yy,t6); shading flat; hold on;
        m_coast('patch',[.7 .7 .7],'edgecolor','k');
        caxis([0 100]); cb = colorbar; ylabel(cb,'DIN (mmol-N/m2)')
        title ('Model DIN over top 95m, Jul 2005','Fontweight','bold')
        
        
% Plot 6 Checks: what are Fe and Si fields like
       
    t1 = modsurfFe(:,:,49:54);
    t2 = nansum(t1,3) ./ 6;
    
    t1 = modsurfSi(:,:,49:54);
    t3 = nansum(t1,3) ./ 6;
    
    subplot(2,1,1) 
        odvpal(10);
        m_proj('stereographic','lat',90,'long',30,'radius',30)
        m_elev('contour',[-500 -500],'edgecolor','r'); hold on 
        m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
        m_pcolor(xx,yy,t2); shading flat; hold on;
        m_coast('patch',[.7 .7 .7],'edgecolor','k');
        caxis([0 0.002]); cb = colorbar; ylabel(cb,'Dissolved Fe (mmolFe/m3)') % July Fe structure: use caxis([0.001 0.0013]).
        title('Dissolved Fe in surface waters (3m), Sep 2005','FontWeight','bold')
        % limitation at Fe < 0.33 x10^-3 (i.e. purple in current plot) -
        % Yool 2010 table 1 - No Fe limitation (and Fe pretty homogenous, slighhtly lower on Atl side) in July.
        % to see Fe
        
        
    subplot(2,1,2) 
        odvpal(10);
        m_proj('stereographic','lat',90,'long',30,'radius',30)
        m_elev('contour',[-500 -500],'edgecolor','r'); hold on 
        m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
        m_pcolor(xx,yy,t3); shading flat; hold on;
        m_coast('patch',[.7 .7 .7],'edgecolor','k');
        caxis([0 20]); cb = colorbar; ylabel(cb,'Silicate (mmolSi/m3)')
        title('Silicate in surface waters (3m), Sep 2005','FontWeight','bold')
        % limitation at Si < 0.75 (i.e. purple in current plot) -
        % Yool 2010 table - some Si limitation in Barents Sea (Atl inflow)
        
        % -> According to model - Arctic is N limited (no Fe or Si limitaion)
        % in growing season, in good agreement with obs.
        
        % Sept fields look very similar to those in July
        


        
% ============
% Notes


% Additional plots


% (1) Logarithmic surface chlorophyll
        
        

%     t1 = log(modsurfchl(:,:,37:42));
%     t2 = nansum(t1,3) ./6; % division is by number of files summed (37-42), to keep units as /day
% 
%     
%     
%     yt = [0.01 0.1 1 10]; % setting colorbar spacing
%     lyt = log10(yt);
%     
%     figure; odvpal(10);
%     m_proj('stereographic','lat',90,'long',30,'radius',30)
%     m_elev('contour',[-500 -500],'edgecolor','r'); hold on
%     m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
%     m_pcolor(xx,yy,t2); shading flat; hold on;
%     m_coast('patch',[.7 .7 .7],'edgecolor','k');
%     caxis([lyt(1) lyt(end)]); cb = colorbar; set(cb, 'YTick', lyt, 'YTickLabel', yt); ylabel(cb, 'log Chl (mg Chl/m2/day)');
%     title ('Model Surface (3m) Chlorophyll, July 2005','FontWeight','bold');
%     
    

        
% (2) Model Chl with ice-cover contour overlain

        % Clear model grid smears over N Atl

        % Set up chl array
        t1 = modsurfchl(:,:,37:42);
        t2 = nansum(t1,3) ./6; % division is by number of files summed (37-42) because files are in units mgChl/m2/day - adding 6 files of daily rate together won't give you the daily rate (averages are sum(individuals)/number(individuals)!)

        % Set up ice cover array
        t1 = modice(:,:,37:42);
        t3 = nansum(t1,3) ./ 6;  

        figure; odvpal(10);
        m_proj('stereographic','lat',90,'long',30,'radius',30)
        m_elev('contour',[-500 -500],'edgecolor','r'); hold on
        m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
        [cs,h]=m_contour(xx,yy,t3,[0.1],'edgecolor','k','LineWidth',2); hold on % 0.1 specifies which contour to draw (i.e. currently draws ice margin at 10%)
        % clabel(cs,h,v,'fontsize',6);            % clabelm(cs,hv) rotates contour labels and inserts them in the contour lines, can spcify vector v which specifies which contours to plot (instead of [0.1] in line above, these are labelled; after this are name-value pairings to specify visual properties
        m_pcolor(xx,yy,t2); shading flat; hold on;
        m_coast('patch',[.7 .7 .7],'edgecolor','k');
        caxis([0 2]); cb = colorbar; ylabel(cb, 'Chl (mg Chl/m2/day)');
        title ('Model Chlorophyll,Jul 2005','FontWeight','bold');