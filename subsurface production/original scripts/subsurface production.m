% Subsurface production plots

% created from 'satellite pp monthly' script on 8/8/14 (jl2g13)


% (1) Reading in monthly model variables - CHL, PP, PAR, SST

% Only considering 2005
yr = 2005; 
    
lattotal = double(1021);
latstart = double(1);
lontotal = double(1442);
depth = double(21); % top 200 m
pos = 1;
    
    % combined
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

        fnamesP = dir(dnomP);
        fnamesD = dir(dnomD);
        fnamesT = dir(dnomT);
        [num,x] = size(fnamesP);
        
%         data(:,m) = size(fnamesP,1); % for plotting need to know number of files with data (non-padded) to maintain /day units


        fprintf('Reading month %d\n', m);

        
        % (1) Chlorophyll (surface)

        fnameP = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/%s', yr, fnamesP.name);
        fprintf('Reading %s\n', fnameP);

        % (i) For non-diatom chlorophyll
        t1 = ncread(fnameP,'CHN',[1 latstart 1 1], [lontotal lattotal depth 1]); % PRN is depth-integrated so a 2D model field                
        t1 (t1 == 0) = NaN; % removing zero values (land)
        chn(:,:,:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)

        % (ii) For diatom chlorophyll
        t1 = ncread(fnameP,'CHD',[1 latstart 1 1], [lontotal (lattotal) depth 1]); % PRD is depth-integrated so a 2D model field
        t1 (t1 == 0) = NaN; % removing zero values (land)
        chd(:,:,:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)       


        % (2) Primary Production

        fnameD = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/%s', yr, fnamesD.name);
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


        % (3) Sea Surface Temperature

        % Extracting model sea surface temperature
        fnameT = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/%s', yr, fnamesT.name);
        fprintf('Reading %s\n', fnameT);

        % Reading in sea surface temperature model field
        t1 = ncread(fnameT,'sosstsst',[1 latstart 1], [lontotal lattotal 1]); % sosstsst is a 2D model field, has no depth dimension
        t1 = squeeze(t1);
        t1 (t1 == 0) = NaN; % removing zero values (land)

        sst(:,:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)


        % (4) Photosynthetically Active Radiation

        % Reading in surface PAR model field
        t1 = ncread(fnameD,'MED_QSR',[1 latstart 1], [lontotal lattotal 1]); % MED_QSR ('top level radiation) is a 2D model field, has no depth dimension
        t1 = squeeze(t1);
        t1 (t1 == 0) = NaN; % removing zero values (land)

        totalshortwave(:,:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)


        pos = pos + 1;

    end
    
modppC = (modndpp + moddpp) * 12.011 * 6.625; % in gC
modchl = chn + chd;
par = totalshortwave .* 0.43;




% ======
% PLOT 1
% ======


ANmodPPC = nansum(modppC,4); % Annual PP (integrating over all months in year)

modppC_DI = nansum(modppC,3); % Whole column PP (integrating over all depths)
modsurfPP = squeeze(modppC(:,:,1:4,:)); % surface PP (3m)

    
 for j = 1:size(modppC,2) % for every 1/4 lattitude
        for i = 1:size(modppC,1) % for every 1/4 longditude
            
            % Annual cross correlation
            a = modppC_DI(i,j,:);   % Jan-Mar = 1-18, Apr-Jun = 19-36, Jul-Sep = 37-54, Oct-Dec 55-73
            b = modsurfppC(i,j,:);
            [c,lags] = xcorr(a,b,0,'coeff'); % zero lag coeff, lags = 0 
            corAN(:,pos) = c;
            
            pos = pos + 1;
            
        end
        fprintf('Reading %d\n', j)
 end


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
    plot1a = nansum(ANmodppC(:,:,1:4),3); % integrate top 4 boxs ~ top 20m
    
    subplot(1,3,1);
    m_proj('stereographic','lat',90,'long',30,'radius',30)
    m_elev('contour',[-500 -500],'edgecolor','r'); hold on
    m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
    m_pcolor(xx,yy,plot1a); shading flat; hold on; % top 4 boxs ~ top 20m
    m_coast('patch',[.7 .7 .7],'edgecolor','k');
    ylabel(c,'Annual PP 0-20 m');
    title(tsrt,'FontWeight','bold')
    
    % annual model PP 20-60 m
    plot1b = nansum(ANmodppC(:,:,5:9),3); % integrate boxes 5--9 ~ 20-60 m depth 
  
    subplot(1,3,2);
    m_proj('stereographic','lat',90,'long',30,'radius',30)
    m_elev('contour',[-500 -500],'edgecolor','r'); hold on
    m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
    m_pcolor(xx,yy,ANmodppC(:,:,5:9); shading flat; hold on; % boxes 5--9 ~ 20-60 m depth 
    m_coast('patch',[.7 .7 .7],'edgecolor','k'); 
    c = colorbar; ylabel(c,'Annual PP 20-60 m');
    title(tsrt,'FontWeight','bold')
    
    % annual XCORR between surface and depth integrated PP
    subplot(1,3,3);
    m_proj('stereographic','lat',90,'long',30,'radius',30)
    m_elev('contour',[-500 -500],'edgecolor','r'); hold on
    m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
    m_pcolor(xx,yy,corAN); shading flat; hold on; % XCORR 
    m_coast('patch',[.7 .7 .7],'edgecolor','k'); 
    c = colorbar; ylabel(c,'Annual XCORR (surface vs depth integrated');
    title(tsrt,'FontWeight','bold')
    
    
    
    
% ======
% PLOT 2
% ======



% ==
% D
% ==

% model surface chlorophyll
modsurfchl = squeeze(nansum(modchl(:,:,1:4,:)),3); % taking depth = 1-20 m, then integrating over this depth and squeezing this dimension. i.e. 'surface' = 0-20 m not 0-3 m
ANmodsurfchl = nansum(modsurfchl,3); % Annual mod surface Chl (integrating over all months)


% ==
% E
% ==


% Reading in satellite chlorophyll
    yr = 2005;
    
    dnomP = sprintf('/noc/users/jl2g13/MATLAB/seawifsCHL/%d/*.hdf', yr);
    fnamesP = dir (dnomP); 
    [num, x] = size(fnamesP);
    satsurfchl = zeros(2160, 1080, num);  % preallocation: lon x lat x num of 8-day averages reading in (note lon and lat both at 1/6�)
    
    pos = 1;
    for f = 1:num
        
        fnameP = sprintf('/noc/users/jl2g13/MATLAB/seawifsCHL/%d/%s', yr, fnamesP(f).name); % for when you start using multiple yrs of satellite data
        fprintf('Reading %s\n', fnameP);
        
        
        % For each file: extracting Chl from hdf file into matlab array (verbose documentation see HDF_to_array.m)
        fileID = hdfsd('start',fnameP, 'read');      
        [numdata, numdescr] = hdfsd('fileinfo', fileID);
        chlid = hdfsd('select', fileID, 0);
        [name,numdim,dimvector,type,numdescr] = hdfsd('getinfo', chlid);
        
        startvector = [0 0]; % Start reading HDF file at the beginning of Chl
        endvector = dimvector; % Finish reading HDF file at end of Chl
        stride = []; % specifies number of data steps to read at (n=1 here, i.e. read all)       
        Chl = hdfsd('readdata', chlid, startvector, stride, endvector);
        Chl = double(Chl);
        ii=find(Chl < 0); Chl(ii)=NaN; % removing missing data flags
        
        

        % Satellite data by month (dim 4 = month)
        % uses reg exp to split file names by month: month given by 'filter'.


        filter = regexp(fnameP,'2005','split'); % hard code: 2005 - has to be string to 
        filter = textscan(filter{3}, '%3s', 1);
        filter = str2num(cell2mat(filter{1}))
        
        satsurfchl(:,:,pos) = Chl;
        
        pos = pos + 1;
        
    end

ochl = flipdim(permute(satsurfchl(:,:,1:46), [2 1 3]), 1);
ochl (ochl == 0) = NaN;

Ansatchl = nansum(ochl,3); % Annual sat Chl (integrating over all months)

% creating satellite grid for M_Map
x1 = (-180 + (1/12)):(1/6):(180 - (1/12));
y1 = (-90 + (1/12)):(1/6):(90 - (1/12));
[x2,y2] = meshgrid(x1, y1);


% ==
% F
% ==

% Calculating vgpm primary productivity from model fields
% vgpm.m (Behrenfeld and Falkowski, 1997)

    for n=1:size(modsurfchl,3) % for each month

        CHL=modsurfchl(:,:,n);   
        SST=sst(:,:,n);
        PAR=par_surface(:,:,n);

        if CHL<1
            Ctot=38*(CHL.^0.425); % integrated water column chl
        else
            Ctot=40.2*(CHL.^0.507);
        end

        Zeu=568.2*(Ctot.^-0.746); % euphotic zone depth, Morel 1989, eq 1a

        if Zeu>102
            Zeu=200*(Ctot.^-0.293); % Morel 1989, eq 1b
        end

        if SST<-1   %% the simple (!) estimate of maximum daily PP
            PBopt=1.13;
        elseif SST>28.5
            PBopt=4;
        else
            PBopt=1.2956 + SST.*2.749e-1 + (SST.^2)*6.17e-2 - (SST.^3)*2.05e-2...
                + (SST.^4)*2.462e-3 - (SST.^5)*1.348e-4 + (SST.^6)*3.4132e-6...
                - (SST.^7)*3.27e-8;
        end

        x=sym(n/12);
        xx=frac(x);
        if xx==0
            N=12;
        else
            N=double(xx)*12;
        end

        DLen= double(24);        % Dirr(:,:,N);  % decimal day length % ARBITRARILY SET, NEED TO CHANGE THIS

        NPPvgpm(:,:,n)=0.66125*(PBopt.*(PAR./(PAR+4.1)).*CHL.*Zeu.*DLen);
        fprintf('NPPvgpm calculated for month %d\n', n)
    end

ANvgpm = nansum(NPPvgpm,3); % integrating over all months

% ==
% G
% ==

% model annual depth integrated PP
ANmodPPz = nansum(nansum(modppC,3),4); % Annual, depth integrated mod PP (integrating over depth and all months)







% Plotting

    figure(2);
    
    % D - Annual model surface chlorophyll
    subplot(2,2,1);
    m_proj('stereographic','lat',90,'long',30,'radius',30)
    m_elev('contour',[-500 -500],'edgecolor','r'); hold on
    m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
    m_pcolor(xx,yy,ANmodsurfchl); shading flat; hold on;
    m_coast('patch',[.7 .7 .7],'edgecolor','k');
    ylabel(cb, 'Model surface chlorophyll');
    title(tsrt,'FontWeight','bold')
    
    
    % E - Annual satellite surface chl
    subplot(2,2,2);
    m_proj('stereographic','lat',90,'long',30,'radius',30)
    m_elev('contour',[-500 -500],'edgecolor','r'); hold on
    m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
    m_pcolor(x2,y2,Ansatchl); shading flat; hold on;
    m_coast('patch',[.7 .7 .7],'edgecolor','k');
    ylabel(cb, 'Annual satellite surface chl');
    title(tsrt,'FontWeight','bold')
    
    % F - Annual depth-int satellite PP (forced by model)
    subplot(2,2,3);
    m_proj('stereographic','lat',90,'long',30,'radius',30)
    m_elev('contour',[-500 -500],'edgecolor','r'); hold on
    m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
    m_pcolor(xx,yy,ANvgpm); shading flat; hold on;
    m_coast('patch',[.7 .7 .7],'edgecolor','k');
    cb = colorbar; ylabel(cb, 'Annual depth-int satellite PP');
    title(tsrt,'FontWeight','bold')
    
    % G - Annual depth-int mod PP
    subplot(2,2,4);
    m_proj('stereographic','lat',90,'long',30,'radius',30)
    m_elev('contour',[-500 -500],'edgecolor','r'); hold on
    m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
    m_pcolor(x2,y2,ANmodPPz); shading flat; hold on;
    m_coast('patch',[.7 .7 .7],'edgecolor','k');
    cb = colorbar; ylabel(cb, 'Annual depth-int satellite PP');
    title(tsrt,'FontWeight','bold')

    