% Model vs Satellite PP



% (1) Reading in model variables

pos = 1; % Initialsing index for where to put data in array


for y = 1:1:1 % Looping over one year -> all years (full possible interval is 1990-2006)
    yr = 2004 + y;
    
    lattotal = double(1021);
    latstart = double(1);
    lontotal = double(1442);
    depth = 1; % 3 m.
    
    
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
    
%     % (4) get a list of I (Ice) files for this year
%     dnomI = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/*d05I.nc', yr);
%     fnamesI = dir (dnomI);
    
    
    
    % Finding number of  files to loop over in each yearly folder (73 for 5 day means; n.b. needs to be same for all file types)
    [num, x] = size(fnamesP);
    
    
    % Preallocating arrays for model output to be read into
    
      modsst = zeros(lontotal,lattotal,num);
      modndpp = zeros(lontotal,lattotal,num);
      moddpp = zeros(lontotal,lattotal,num);
      modsurfchn = zeros(lontotal,lattotal,depth,num);
      modsurfchd = zeros(lontotal,lattotal,depth,num);
      modsurfpar = zeros(lontotal,lattotal,num);
%       modice = zeros(lontotal,lattotal,num);


    
    % Initiating extraction loop
    
    for f = 1:1:num % looping from file 1 to num (73) - for each year
        
        % (1) For Chlorophyll - saved modchl variable is sum of top 13 boxes (~95 m)
        
        % Extracting model surface chlorophyll
        fnameP = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/%s', yr, fnamesP(f).name);
        fprintf('Reading %s\n', fnameP);
        
        % Reading in Surface Chl model field
        % (i) For non-diatom chlorophyll
        
        t1 = ncread(fnameP,'CHN',[1 latstart 1 1], [lontotal lattotal depth 1]); % can set this to read in just N of Arctic Circle
        t1 = squeeze(t1);
        t1 (t1 == 0) = NaN; % removing zero values (land)
        
        modsurfchn(:,:,:,pos) = t1; % lat, lon, no of 5-day means (73/yr)
        
        % (ii) For diatom chlorophyll
        
        t1 = ncread(fnameP,'CHD',[1 latstart 1 1], [lontotal lattotal depth 1]);
        t1 = squeeze(t1);
        t1 (t1 == 0) = NaN; % removing zero values (land)
        
        modsurfchd(:,:,:,pos) = t1;
          
        
        
        % (2) Primary Production 
        
        fnameD = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/%s', yr, fnamesD(f).name);
        fprintf('Reading %s\n', fnameD);

        % (a) Model depth-integrated
                
        % (i) For non-diatom production
        
        t1 = ncread(fnameD,'PRN',[1 latstart 1], [lontotal lattotal 1]); % PRN is depth-integrated so a 2D model field
        t1 = squeeze(t1);
        t1 (t1 == 0) = NaN; % removing zero values (land)
        
        modndpp(:,:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)
        
        % (ii) For diatom production
        t1 = ncread(fnameD,'PRD',[1 latstart 1], [lontotal lattotal 1]); % PRD is depth-integrated so a 2D model field
        t1 = squeeze(t1);
        t1 (t1 == 0) = NaN; % removing zero values (land)
        
        moddpp(:,:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)
        
        
        
        % (4) For Sea Surface Temperature
        
        % Extracting model sea surface temperature
        fnameT = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/%s', yr, fnamesT(f).name);
        fprintf('Reading %s\n', fnameT);
        
        % Reading in sea surface temperature model field
        t1 = ncread(fnameT,'sosstsst',[1 latstart 1], [lontotal lattotal 1]); % sosstsst is a 2D model field, has no depth dimension
        t1 = squeeze(t1);
        t1 (t1 == 0) = NaN; % removing zero values (land)
        
        modsst(:,:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)
        
        
        % (5) For Photosynthetically Active Radiation
        
        % Reading in surface PAR model field
        t1 = ncread(fnameD,'MED_QSR',[1 latstart 1], [lontotal lattotal 1]); % MED_QSR ('top level radiation) is a 2D model field, has no depth dimension
        t1 = squeeze(t1);
        t1 (t1 == 0) = NaN; % removing zero values (land)
        
        totalshortwave(:,:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)
        
        
        
        
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
        
                 
        pos = pos + 1; % resetting position counter for next iteration
         
    end
end


modsurfchl = squeeze(modsurfchd + modsurfchn); % adding diatom and non-diatom chl components
modpp = modndpp + moddpp; % adding diatom and non-diatom primary production components
modppC = modpp * 12.011 * 6.625; % conversion from mmolN (model) to mgC (satellites) is * 12.011 (g/mol) * 6.625 (assuming Redfieldian)
modsurfpar = totalshortwave .* 0.43; % Fudge conversion from total shortwave radiation to PAR (used previously by Andrew)


% Emptying memory
clearvars totalshortwave modsurfchd modsurfchn modndpp moddpp modpp










% (2) Reading in satellite variables


pos = 1;

 for y = 1:1:1         % for using multiple years of satellite data, seawifs data for end1997-mid2010
     yr = 2004 + y
    
    
    dnomP = sprintf('/noc/users/jl2g13/MATLAB/seawifsCHL/%d/*.hdf', yr);
    fnamesP = dir (dnomP); 
    [num, x] = size(fnamesP);
    satsurfchl = zeros(2160, 1080, num);  % preallocation: lon x lat x num of 8-day averages reading in (note lon and lat both at 1/6°)
    
    
    for f = 1:num
        
      fnameP = sprintf('/noc/users/jl2g13/MATLAB/seawifsCHL/%d/%s', yr, fnamesP(f).name); % for when you start using multiple yrs of satellite data
      % fnameP = fnamesP(f).name;
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
        
        satsurfchl(:,:,pos) = Chl;
        
        pos = pos + 1;
        
    end
 end

ochl = flipdim(permute(satsurfchl, [2 1 3]), 1);
ochl (ochl == 0) = NaN; % setting land positions to zero
% figure; h = pcolor(log(ochl(:,:,20))); shading flat; colorbar; title('ochl t = 20') % checking plot








% (3) Plotting Chlorophyll

% (i) Preamble

    ft='/noc/altix2/scratch/omfman/ORCA025-N201/means/1990/ORCA025-N201_1990m09P.nc';
    xx=ncread(ft,'nav_lon', [1 1], [lontotal lattotal]);
    yy=ncread(ft,'nav_lat', [1 1], [lontotal lattotal]);


    % removing spurious spikes in N Atl when plotting (beware true zeros, quick fix)
    land = xx(1,1);
    xx(xx == land) = NaN; xx(xx == 0) = NaN;
    land = yy(1,1);
    yy(yy == land) = NaN; yy(yy == 0) = NaN;


    % creating satellite grid for M_Map
    x1 = (-180 + (1/12)):(1/6):(180 - (1/12));
    y1 = (-90 + (1/12)):(1/6):(90 - (1/12));
    [x2,y2] = meshgrid(x1, y1);



% (ii) Model
    
    t1 = modsurfchl(:,:,37:42);
    t2 = nansum(t1,3) ./6; % division is by number of files summed (37-42) because files are in units mgChl/m2/day - adding 6 files of daily rate together won't give you the daily rate (averages are sum(individuals)/number(individuals)!)
    
    figure(1); odvpal(10);
    subplot(1,2,1);
    m_proj('stereographic','lat',90,'long',30,'radius',30)
    m_elev('contour',[-500 -500],'edgecolor','r'); hold on
    m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
    m_pcolor(xx,yy,t2); shading flat; hold on;
    m_coast('patch',[.7 .7 .7],'edgecolor','k');
    caxis([0 2]); cb = colorbar; ylabel(cb, 'Chl (mg Chl/m2/day)');
    title ('Model Chlorophyll (to 3m), Jul 2005','FontWeight','bold');
     
       
    
 % (iii) Satellite
        
    t1 = ochl(:,:,24:27); % Mar t = 9-12, July t = 24-27, Sep t = 32-35 // this line reverses setting ochl (ochl == 0) = NaN; because at each land position summing along n nans - all of which are ignored, so for that position nansum = 0
    t2 = nansum(t1,3) ./ 4; % division is by number of files summed (24-27), to keep units as /day 
    t2 (t2 == 0) = NaN;

    subplot(1,2,2);
    m_proj('stereographic','lat',90,'long',30,'radius',30)
    m_elev('contour',[-500 -500],'edgecolor','r'); hold on
    m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
    m_pcolor(x2,y2,t2); shading flat; hold on;
    m_coast('patch',[.7 .7 .7],'edgecolor','k'); 
    caxis ([0 5]); c = colorbar; ylabel(c,'Chl (mg Chl/m2/day)');
    title('Satellite Chlorophyll (z90), Jul 2005','FontWeight','bold')
    
    
  % Emptying memory
    clearvars satsurfchl ochl
    
    
    
    
    
    
    
    
    
    
        
% (4) Satellite primary production


% (i) vgpm.m (Behrenfeld and Falkowski, 1997)


    % calculate


    for n=1:size(modsurfchl,3)

        CHL=modsurfchl(:,:,n);   
        SST=modsst(:,:,n);
        PAR=modsurfpar(:,:,n);

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
        fprintf('NPPvgpm calculated for 5-day mean %d\n', n)
    end

    

   
    
% (ii) carr.m (Carr 2002; dev. eastern boundary currents)

alpha=2.64;

for n=1:size(modsurfchl,3)
    
    CHL=modsurfchl(:,:,n);
    SST=modsst(:,:,n);
    % PAR=(modsurfpar(:,:,n)*6e23)/(86400*2.77e18);    % converts PAR from Einsteins to Watts/m2
    PAR = modsurfpar(:,:,n); %  model par already in W/m2 so no need for conversion
    Kpar=0.04+(0.0088*CHL)+(0.054*(CHL.^0.66)); % attenuation coefficient of PAR, based on chl content of water
    
    if CHL<1
        Ctot=38*(CHL.^0.425);  % integrated water column chl
    else
        Ctot=40.2*(CHL.^0.507);
    end
    Zeu=568.2*(Ctot.^-0.746);  % euphotic depth
    if Zeu>102
        Zeu=200*(Ctot.^-0.293);
    end
    D=Zeu;
    Ed=(PAR.*(1-exp(-Kpar.*D)))./(Kpar.*D);  % downwelling irradiance
    Pmax=24*exp(0.09*SST);   % maximum photosynthetic rate (Eppley, 1972)
    P=CHL.*(((Pmax.*Ed)./((Pmax./alpha)+Ed))).*D;
    
    NPPcarr(:,:,n)=P;
    fprintf('NPPcarr calculated for 5-day mean %d\n', n)
end

NPPcarr(find(NPPcarr==0))=NaN;
   
    
    
    
  
    
% (iii) marra.m (Marra 2003)    

NPPmarra = zeros(lontotal,lattotal,num); % preallocating


for n=1:size(modsurfchl,3) % note is the number of 5 day means (73/yr), equal to 'num'
    
    SST = modsst(:,:,n);
    CHL = modsurfchl(:,:,n);
    PAR = modsurfpar(:,:,n);
    
    
    kc=0.00433*exp(0.08249*SST);  % various light attenuation coefficients
    k4=find(SST<12);
    kc(k4)=0.0105+SST(k4)*0.0001;
    kc_p=0.00433*exp(0.85*0.08249*SST);
    k4=find(SST<14.475);
    kc_p(k4)=0.0105+SST(k4)*0.0001;
    kw=0.04;  kx=0.02;
    
    z=0:2:98;  %% depths
    
    for j=1:size(CHL,1)
        for k=1:size(CHL,2)     %  only calculating PP for Arctic to save computation time (825 = 65°N) 
            cc=CHL(j,k);
            
            if cc>0.4                   % HERE is where you would adjust sensu ardyna
                h=42;sig=50;zm=5;y=20;  %% parameters defining shape of curve of chl vs depth
            else
                h=20;sig=18;zm=75;y=1;
            end
            
            cz=cc+(h/sig/sqrt(2*pi))*exp(-(z-zm).^2/2/sig^2)/y;
            
            kpar=kw+kc(j,k)*cz+kx;  %% attenuation of PAR
            
            E=PAR(j,k);
            for ii=2:50
                E(ii)=E(ii-1)*exp(-kpar(ii)*2);
            end
            
            EE=E/PAR(j,k);
            dp=min(find(EE<0.01));
            
            prod=12*2*0.03*(10./(10+E))*kc_p(j,k).*cz.*E;
            
            NPPmarra(j,k)=sum(prod(1:dp));  %% production integrated to the 1% light level
        end
        fprintf('%d\n', j); 
    end   
    fprintf('Read 5day file %d\n', n)
end
    NPPmarra(find(NPPmarra==0))=NaN;
   

    
    
    
    
% Plotting primary production    

figure(2);

% (i) Model

    t1 = modppC(:,:,37:42); % go to the directory and count which 5 day av are within July.
    t2 = nansum(t1,3) ./ 6; % % division is by number of files summed (37-42), to keep units as /day
    % test = t1(:,:,1) ./ (t2 ./ 6); % check it looks alright - values should be around 1 (plot it using m_map bundle below)

    subplot(2,2,1)
    m_proj('stereographic','lat',90,'long',30,'radius',30)
    m_elev('contour',[-500 -500],'edgecolor','r'); hold on % Arctic shelf break is at 220 m. 
    m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
    m_pcolor(xx,yy,t2); shading flat; hold on;
    m_coast('patch',[.7 .7 .7],'edgecolor','k');
    caxis([0 1500]); cb = colorbar; ylabel(cb,'PP (mgC/m2/day)')
    title('Model depth-integrated PP, Jul 2005','FontWeight','bold')
    
% (ii) vgpm
    
    NPPvgpm(find(NPPvgpm==0)) = NaN;
    % figure; pcolor(permute(NPPvgpm(:,:,40), [2 1])); shading flat; colorbar; title('vgpm depth-integrated pp, t = 40');
    
    t1 = NPPvgpm(:,:,37:42); % Mar t = 19-24, Jul t = 37-42, Sep t = 49-54
    t2 = nansum(t1,3) ./ 6; % division is by number of files summed (37-42), to keep units as /day

    subplot(2,2,2)
    m_proj('stereographic','lat',90,'long',30,'radius',30)
    m_elev('contour',[-500 -500],'edgecolor','r'); hold on
    m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
    m_pcolor(xx,yy,t2); shading flat; hold on;
    m_coast('patch',[.7 .7 .7],'edgecolor','k'); caxis ([0 5000]); cb = colorbar; ylabel(cb, 'PP (mgC/m2/day)');
    title('VGPM depth-integrated PP, Jul 2005','FontWeight','bold')
    
% (iii) carr

    t1 = NPPcarr(:,:,37:42); % Mar t = 19-24, Jul t = 37-42, Sep t = 49-54
    t2 = nansum(t1,3) ./6; % division is by number of files summed (37-42), to keep units as /day

    subplot(2,2,3)
    m_proj('stereographic','lat',90,'long',30,'radius',30)
    m_elev('contour',[-500 -500],'edgecolor','r'); hold on
    m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
    m_pcolor(xx,yy,t2); shading flat; hold on;
    m_coast('patch',[.7 .7 .7],'edgecolor','k'); caxis ([0 5000]); cb = colorbar; ylabel(cb, 'PP (mgC/m2/day)');
    title('Carr depth-integrated PP, Jul 2005','FontWeight','bold')
    
% (iv) marra

    t1 = NPPmarra(:,:,37:42); % Mar t = 19-24, Jul t = 37-42, Sep t = 49-54
    t2 = nansum(t1,3) ./ 6; % division is by number of files summed (37-42), to keep units as /day

    subplot(2,2,4)
    m_proj('stereographic','lat',90,'long',30,'radius',30)
    m_elev('contour',[-500 -500],'edgecolor','r'); hold on
    m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
    m_pcolor(xx,yy,t2); shading flat; hold on;
    m_coast('patch',[.7 .7 .7],'edgecolor','k'); caxis ([0 5000]); cb = colorbar; ylabel(cb, 'PP (mgC/m2/day)');
    title('Marra depth-integrated PP, Jul 2005','FontWeight','bold')        
    
    
    
    
    
     
    
    
% ============
% NOTES

% (1) to save matlab fig to eps and print in colour
    %1) save matlab figure
    print -depsc -painters satchljul % painters = save as vector image (not bitmap). depsc = eps_colour
    % 2) check in ghostscript and print
    gs test_01.eps % view image in ghostscript (i.e. execute outside matlab enviro)
    lpr -Pspartacus % print to spartacus
    lpq -Pspartacus % just to check it's printing
    
    
    
% (2) Daylength for all latitudes for whole year

% calculates daylength for 180 degrees of longitude for 365 days of the
% year and plots daylength on a (y = lat,x = year day) grid.


for l = 1:1:180
    lat = -91 + l;
    for n = 1:1:365
        dlen(l, n) = daylength(lat, n); % using function daylength
    end
end

figure; pcolor(dlen); shading flat; cb = colorbar;
ylabel(cb, 'Daylength (hr)'); title('Daylength(lat,yDay)','FontWeight','bold')




% (3) Additional algorithms


% (i) vgpm_arc (Hill and Zimmerman, 2010)

% CURRENTLY A PROBLEM BETWEEN DAYS 30 and 60 WHERE NPP IS UNIFORMLY LARGE
% -VE NUMBERS (ORDER -10^4)
% Bear in mind Hill and Zimmerman never tested it for pan-arctic.

% Hill and Zimmerman make two changes to the original vgpm to tune it for
% the Arctic. (1) use R -> Chl relationship of Cota et al, 2004 and (2)
% change the temperature dependency of PBopt
%   Since I've taken Chl from the model I'm only changing the temperature
%   dependency of PBopt. Note this new dependecy is derived by fitting to
%   Chukchi Sea data only.


for n=1:size(modsurfchl,3)
    
    CHL=modsurfchl(:,:,n);
    SST=modsst(:,:,n);
    PAR=modsurfpar(:,:,n);
    
    if CHL<1
        Ctot=38*(CHL.^0.425); % integrated water column chl
    else
        Ctot=40.2*(CHL.^0.507);
    end
    
    Zeu=568.2*(Ctot.^-0.746); % euphotic zone depth
    
    if Zeu>102
        Zeu=200*(Ctot.^-0.293);
    end
    
    if SST<-1   %% the simple (!) estimate of maximum daily PP
        PBopt=1.13;
    elseif SST>28.5
        PBopt=4;
    else
        PBopt=1.44983 + SST.*5.3733e-1 - (SST.^2)*1.051e-1 - (SST.^3)*3.877e-2...
            + (SST.^4)*1.036e-2 - (SST.^5)*5.9e-4; % Hill,2010 eq 9. NB NOT entirely clear in H2010 whether they apply the polynomial to SST between 1 and 28.5 (as B97) or to all SSTs
    end
    
    x=sym(n/12);
    xx=frac(x);
    if xx==0
        N=12;
    else
        N=double(xx)*12;
    end
    
    DLen= double(12);        % Dirr(:,:,N);  % decimal day length % ARBITRARILY SET, NEED TO CHANGE THIS
    
    NPPvgpmA(:,:,n)=0.66125*(PBopt.*(PAR./(PAR+4.1)).*CHL.*Zeu.*DLen);
    fprintf('NPPvgpm calculated for 5-day mean %s\n', n)
end

NPPvgpmA(find(NPP==0))=NaN;

% checking it looks alright

figure; pcolor(NPPvgpmA(:,:,40)); shading flat; colorbar; title('vgpm depth-integrated pp');




% (ii) Simple Biomass Model (Hill & Zimmerman 2010, ammended by Hill 2013)



% for chl(z) constant

z = 100;

for n = 1: size(modsurfchl,3) % time intervals
    for i = 1:size(modsurfchl,1)
        for j = 1:size(modsurfchl,2)
        ipp = 1.04 + 0.86 * log(modsurfchl(:,:,n)) * z; % Hill 2013 p 111
        end
        fprintf('%d\n', i)
    end
end  





% (4) Additional code

    % Satellite chl on logarithmic scale
% 
%     % Calculating Chl
%     t1 = log10(ochl(:,:,24:27)); % Mar t = 9-12, July t = 24-27, Sep t = 32-35
%     t2 = nansum(t1,3);
% 
%     %t1 = log10(ochl(:,:,20));
%     % setting colorbar spacing
%     yt = [0.01 0.1 1 10];
%     lyt = log10(yt);
% 
%     % plotting (4) of Katya monthly plots
%     figure; clf; odvpal(10);
%     m_proj('stereographic','lat',90,'long',30,'radius',30)
%     m_elev('contour',[-500 -500],'edgecolor','r'); hold on
%     m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
%     m_pcolor(x2,y2,t2); shading flat; hold on;
%     m_coast('patch',[.7 .7 .7],'edgecolor','k'); 
%     caxis ([lyt(1) lyt(end)]); cb = colorbar; set(cb, 'YTick', lyt, 'YTickLabel', yt); ylabel(cb,'log Chl (mg Chl/m2)');
%     title('Satellite Chlorophyll (z90), Jul 2005','FontWeight','bold')   
