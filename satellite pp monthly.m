% Monthly satellite PP


% (1) Reading in model variables

yr = 2005;
    
lattotal = double(1021);
latstart = double(1);
lontotal = double(1442);
depth = double(1); % 3 m
pos = 1;
    
    % combined
    for m = 1:12
        if m < 10
            dnomP = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/ORCA025-N201_%d0%d*d05P.nc', yr, yr, m);
            dnomD = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/ORCA025-N201_%d0%d*d05D.nc', yr, yr, m);
            dnomT = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/ORCA025-N201_%d0%d*d05T.nc', yr, yr, m);
        else
            dnomP = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/ORCA025-N201_%d%d*d05P.nc', yr, yr, m);
            dnomD = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/ORCA025-N201_%d%d*d05D.nc', yr, yr, m);
            dnomT = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/ORCA025-N201_%d%d*d05T.nc', yr, yr, m);
        end

        fnamesP = dir(dnomP);
        fnamesD = dir(dnomD);
        fnamesT = dir(dnomT);
        [num,x] = size(fnamesP);
        
        data(:,m) = size(fnamesP,1); % for plotting need to know number of files with data (non-padded) to maintain /day units


        fprintf('Reading month %d\n', m);

            for f = 1:num % ie for each of the files within that given month - read in modpp
                
                % (1) Chlorophyll (surface)
                
                fnameP = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/%s', yr, fnamesP(f).name);
                fprintf('Reading %s\n', fnameP);

                % (i) For non-diatom chlorophyll
                t1 = ncread(fnameP,'CHN',[1 latstart 1 1], [lontotal lattotal depth 1]); % PRN is depth-integrated so a 2D model field                
                t1 (t1 == 0) = NaN; % removing zero values (land)
                chn(:,:,:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)

                % (ii) For diatom chlorophyll
                t1 = ncread(fnameP,'CHD',[1 latstart 1 1], [lontotal (lattotal) depth 1]); % PRD is depth-integrated so a 2D model field
                t1 (t1 == 0) = NaN; % removing zero values (land)
                chd(:,:,:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)       
                
                
                % (2) Primary Production (depth integrated)
                
                fnameD = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/%s', yr, fnamesD(f).name);
                fprintf('Reading %s\n', fnameD);
                
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
                
                
                % (3) Sea Surface Temperature
        
                % Extracting model sea surface temperature
                fnameT = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/%s', yr, fnamesT(f).name);
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

            
            % Padding - making arrays for each individual month the same size
            % so that they can be sequentially appended to form a single 4D
            % array. Pad value = NaN
            chd = padarray(chd, [0 0 0 (7-num)], NaN, 'post'); % to be able to append each month to the same array requires that each appended array is the same size as the matrix it's being appended to - 7 is the most files that appear in a given month (December)
            chn = padarray(chn, [0 0 0 (7-num)], NaN, 'post');
            modndpp = padarray(modndpp, [0 0 (7-num)], NaN, 'post');
            moddpp = padarray(moddpp, [0 0 (7-num)], NaN, 'post');
            sst = padarray(sst, [0 0 (7-num)], NaN, 'post');
            totalshortwave = padarray(totalshortwave, [0 0 (7-num)], NaN, 'post');
           
            
            % assign current month to array that will store all months for a given var
            modsurfchl(:,:,:,:,m) = chn + chd; % remove depth dimension (=1)
            modppC(:,:,:,m) = (modndpp + moddpp) * 12.011 * 6.625; % conversion from mmolN (model) to mgC (satellites) is * 12.011 (g/mol) * 6.625 (assuming Redfieldian)
            par_surface(:,:,:,m) = totalshortwave .* 0.43;
            modsst(:,:,:,m) = sst;

            clearvars chn chd moddpp modndpp totalshortwave sst % clear all variables read in within f-loop (otherwise when come to assign modice on second,third etc loops it can't assign if the num of files in the previous month isn't the same as the month it's overwriting with
            pos = 1; % must reset for the next month (loop): modndpp and modpp then overwrite preceeding months matricies (but only if num(monthN) >= num(month(N-1)) - this is why all model vars read in within f loop must be cleared

    end
    
modsurfchl = squeeze(modsurfchl); % removing depth dimension (=1)

    
    
    
% (2) Reading in satellite chlorophyll
    yr = 2005;
    
    dnomP = sprintf('/noc/users/jl2g13/MATLAB/seawifsCHL/%d/*.hdf', yr);
    fnamesP = dir (dnomP); 
    [num, x] = size(fnamesP);
    satsurfchl = zeros(2160, 1080, num);  % preallocation: lon x lat x num of 8-day averages reading in (note lon and lat both at 1/6°)
    
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
        filter = str2num(cell2mat(filter{1})) % converts cell to num 

            % PLAN A
%             fprintf('Jan file\n');
%             chl_Jan(:,:,pos) = Chl;
%             elseif filter > 31 & filter <= 59
%                 fprintf('Feb file\n');
%                 chl_Feb(:,:,pos) = Chl;
%             elseif filter > 59 & filter <= 90
%                 fprintf('Mar file\n');
%                 chl_Mar(:,:,pos) = Chl;
%             elseif filter > 90 & filter <= 120
%                 fprintf('Apr file\n');
%                 chl_Apr(:,:,pos) = Chl;
%             elseif filter > 120 & filter <= 151
%                 fprintf('May file\n');
%                 chl_May(:,:,pos) = Chl;
%             elseif filter > 151 & filter <= 181
%                 fprintf('Jun file\n');
%                 chl_Jun(:,:,pos) = Chl;
%             elseif filter > 181 & filter <= 212
%                 fprintf('Jul file\n');
%                 chl_Jul(:,:,pos) = Chl;
%             elseif filter > 212 & filter <= 243
%                 fprintf('Aug file\n');
%                 chl_Aug(:,:,pos) = Chl;
%             elseif filter > 243 & filter <= 273
%                 fprintf('Sep file\n');
%                 chl_Sep(:,:,pos) = Chl;
%             elseif filter > 273 & filter <= 304
%                 fprintf('Oct file\n');
%                 chl_Oct(:,:,pos) = Chl;
%             elseif filter > 304 & filter <= 334
%                 fprintf('Nov file\n');
%                 chl_Nov(:,:,pos) = Chl;
%             else
%                 fprintf('Dec file\n');
%                 chl_Dec(:,:,pos) = Chl;
%         end
%         
        satsurfchl(:,:,pos) = Chl;
        
        % PLAN B
        % satsurfchl(:,:,pos,filter) = Chl;
        % months_filters = [0-31, 32-.. etc]  
        
        %for i = 1:length(month_filters) % ie for each month
        % temp = find(satsurfchl(:,:,:,filter(i))
        % satsurfchl(:,:,?,filer(i) = temp
        %end
         
        
        pos = pos + 1;
        
    end

ochl = flipdim(permute(satsurfchl(:,:,1:46), [2 1 3]), 1);
ochl (ochl == 0) = NaN;




% (3) General plotting preamble

% Setting up M_Map model and satellite grids
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
    
% setting plot titles
    mnom = char('May','Jun','Jul','Aug','Sep'); % plot titles
    pnom = char('Pre','Post','Winter');

    
    

% (4) Plotting Chlorophyll

% (a) By phenology bins


% (i) Model    
    preM = nansum(nansum(modsurfchl(:,:,:,2:4),3),4) / sum(data(2:4)); % feb-apr. data: number of files with stuff in, i.e. non-padded
    postM = nansum(nansum(modsurfchl(:,:,:,5:9),3),4) / sum(data(5:9)); % may-sept
    winterM = nansum(nansum(modsurfchl(:,:,:,10:12),3),4) / sum(data(10:12)); % oct-dec
    
    PhenM_cell = [{preM} {postM} {winterM}];
    
% (ii) Satellite

    preS = nansum(ochl(:,:,5:15),3) / 11; % 5-15 covers Feb-Apr, which has 11 data files
    preS (preS ==0) = NaN;

    postS = nansum(ochl(:,:,16:35),3) / 20; % 16-35 is May-Sep, which has 20 data files
    postS (postS ==0) = NaN;
    
    winterS = nansum(ochl(:,:,36:46),3) / 11; % 36-36 is Oct-Dec, which has 11 data files
    winterS (winterS == 0) = NaN;

    PhenS_cell = [{preS} {postS} {winterS}];

    % converting 'PhenM' and 'PhenS' from cell array to normal array        
    pos = 1;
    for i = 1:length(PhenS_cell)
        
        % model
        PhenM(:,:,pos) = cell2mat(PhenM_cell(i));
        
        % satellite
        PhenS(:,:,pos) = cell2mat(PhenS_cell(i)); % matlab can't directly convert cell to double - have to 'force' it to with cell2mat
        pos = pos + 1;
    end
    
    
% (iii) Plotting

        % want less verbose code so to loop over plotting in M_Map    
    for i = 1:size(PhenM,3)
        
        plot_mod = squeeze(PhenM(:,:,i));
        plot_sat = squeeze(PhenS(:,:,i));
                
        figure(i); odvpal(10);
        subplot(1,2,1);
        suptitle('Chlorophyll');
        m_proj('stereographic','lat',90,'long',30,'radius',30)
        m_elev('contour',[-500 -500],'edgecolor','r'); hold on
        m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
        m_pcolor(xx,yy,plot_mod); shading flat; hold on;
        m_coast('patch',[.7 .7 .7],'edgecolor','k');
        caxis([0 2]); cb = colorbar; ylabel(cb, 'Chl (mg Chl/m2/day)');
        tsrt = sprintf('Model %s', pnom(i,:));
        title(tsrt,'FontWeight','bold')
        
        subplot(1,2,2);
        m_proj('stereographic','lat',90,'long',30,'radius',30)
        m_elev('contour',[-500 -500],'edgecolor','r'); hold on
        m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
        m_pcolor(x2,y2,plot_sat); shading flat; hold on;
        m_coast('patch',[.7 .7 .7],'edgecolor','k');
        caxis([0 2]); cb = colorbar; ylabel(cb, 'Chl (mg Chl/m2/day)');
        tsrt = sprintf('Satellite %s', pnom(i,:));
        title(tsrt,'FontWeight','bold')
    end
    
    
    
    
% (b) By month bins    

% (i) Model
    
    pos = 1;
    for m = 1:12
        MonthM(:,:,pos) = nansum(modsurfchl(:,:,:,m),3) / data(m);   % model post-bloom monthly means
        pos = pos + 1;
    end
    
    
% (ii) Satellite
        
        JanS = nansum(ochl(:,:,1:4),3) / 4;
        JanS (JanS ==0) = NaN;
        
        FebS = nansum(ochl(:,:,5:8),3) / 4;
        FebS (FebS ==0) = NaN;
        
        MarS = nansum(ochl(:,:,9:12),3) / 4;
        MarS (MarS ==0) = NaN;
        
        AprS = nansum(ochl(:,:,13:15),3) / 3;
        AprS (AprS ==0) = NaN;

        MayS = nansum(ochl(:,:,16:19),3) / 4;
        MayS (MayS ==0) = NaN;
        
        JunS = nansum(ochl(:,:,20:23),3) / 4;
        JunS (JunS ==0) = NaN;
        
        JulS = nansum(ochl(:,:,24:27),3) / 4;
        JulS (JulS ==0) = NaN;
        
        AugS = nansum(ochl(:,:,28:31),3) / 4;
        AugS (AugS ==0) = NaN;
        
        SepS = nansum(ochl(:,:,32:35),3) / 4;
        SepS (SepS ==0) = NaN;
        
        OctS = nansum(ochl(:,:,36:38),3) / 3;
        OctS (OctS ==0) = NaN;
        
        NovS = nansum(ochl(:,:,39:42),3) / 4;
        NovS (NovS ==0) = NaN;
        
        DecS = nansum(ochl(:,:,43:46),3) / 4;
        DecS (DecS ==0) = NaN;
        
        
        MonthS_cell = [{JanS} {FebS} {MarS} {AprS} {MayS} {JunS} {JulS} {AugS} {SepS} {OctS} {NovS} {DecS}]; % Monthly sat chl
        
        % converting 'MonthS' from cell array to normal array        
        pos = 1;
        for i = 1:length(MonthS_cell)
            MonthS(:,:,pos) = cell2mat(MonthS_cell(i)); % matlab can't directly convert cell to double - have to 'force' it to with cell2mat
            pos = pos + 1;
        end


        
% (iii) plotting
%    each month on a seperate plot as model-satellite pairs
        
        for i = 1:size(MonthS,3)
            figure(i); odvpal(10);
            
            
            plot_sat = squeeze(MonthS(:,:,i));
            plot_mod = squeeze(MonthM(:,:,i));
            
            % Model
            subplot(1,2,1);
            suptitle('Chlorophyll'); % plots title over 
            m_proj('stereographic','lat',90,'long',30,'radius',30)
            m_elev('contour',[-500 -500],'edgecolor','r'); hold on
            m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
            m_pcolor(xx,yy,MonthM(:,:,i)); shading flat; hold on;
            m_coast('patch',[.7 .7 .7],'edgecolor','k'); 
            caxis ([0 2]); c = colorbar; ylabel(c,'Chl (mg Chl/m2/day)');
            tsrt = sprintf('Model %s', mnom(i,:));
            title(tsrt,'FontWeight','bold')
            
            % Satellite
            subplot(1,2,2);
            m_proj('stereographic','lat',90,'long',30,'radius',30)
            m_elev('contour',[-500 -500],'edgecolor','r'); hold on
            m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
            m_pcolor(x2,y2,MonthS(:,:,i)); shading flat; hold on;
            m_coast('patch',[.7 .7 .7],'edgecolor','k'); 
            caxis ([0 2]); c = colorbar; ylabel(c,'Chl (mg Chl/m2/day)');
            tsrt = sprintf('Satellite %s', mnom(i,:));
            title(tsrt,'FontWeight','bold')
        end
             
   
% (iv) Bar plot
        
        % averaging Pan_Arctic
        pos = 1;
        for i= 1:size(MonthM,3)
            
            % model
            mod_bar(pos) = nansum(nansum(MonthM(:,:,i),1),2) / (size(MonthM,1) * size(MonthM,2));

            % satellite
            sat_bar(pos) = nansum(nansum(MonthS(:,:,i),1),2) / (size(MonthS,1) * size(MonthS,2));
            pos = pos + 1;        
        end
        
        sat
        colormap(jet)
        plot_bar = permute([mod_bar; sat_bar], [2 1]);
        figure(1); bar(plot_bar); legend('mod chl','sat chl'); xlabel('Month (post = May-Sep)'); ylabel ('Chl (mgChl/m3)')
        title('Pan-Arctic average Chl concentrations');
    
  % Emptying memory
    clearvars satsurfchl ochl
    
    
    
    
    
    
    
    
    
% (5) Calculating primary productivity   



% (i) vgpm.m (Behrenfeld and Falkowski, 1997)


    for n=1:size(modsurfchl,4)

        CHL=modsurfchl(:,:,:,n);   
        SST=modsst(:,:,:,n);
        PAR=par_surface(:,:,:,n);

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

        NPPvgpm(:,:,:,n)=0.66125*(PBopt.*(PAR./(PAR+4.1)).*CHL.*Zeu.*DLen);
        fprintf('NPPvgpm calculated for month %d\n', n)
    end
    
    
    % (ii) carr
    
    alpha=2.64;

for n=1:size(modsurfchl,4)
    
    CHL=modsurfchl(:,:,:,n);
    SST=modsst(:,:,:,n);
    % PAR=(modsurfpar(:,:,n)*6e23)/(86400*2.77e18);    % converts PAR from Einsteins to Watts/m2
    PAR = par_surface(:,:,:,n); %  model par already in W/m2 so no need for conversion
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
    
    NPPcarr(:,:,:,n)=P;
    fprintf('NPPcarr calculated for month %d\n', n)
end

NPPcarr(find(NPPcarr==0))=NaN;

    
    
    
    
    
% (6) Plotting primary production


% (a) by month

% (i) model

    pos = 1;
    for m = 1:12
        MonthPPM(:,:,pos) = nansum(modppC(:,:,:,m),3) / data(m);   % model post-bloom monthly means
        pos = pos + 1;
    end   


% (ii) satellite
    
    pos = 1;
    for m = 1:12
        MonthPPS(:,:,pos) = nansum(NPPvgpm(:,:,:,m),3) / data(m);   % model post-bloom monthly means
        MonthPPS2(:,:,pos) = nansum(NPPcarr(:,:,:,m),3) / data(m); 
        pos = pos + 1;
    end
    
    
% (iii) plotting

    for i = 1:5
 
        plot_sat = squeeze(MonthPPS(:,:,i));
        plot_mod = squeeze(MonthPPM(:,:,i));
        
        figure(i); odvpal(10);

        subplot(1,2,1);
        suptitle('Primary Prodution'); % plots title over 
        m_proj('stereographic','lat',90,'long',30,'radius',30)
        m_elev('contour',[-500 -500],'edgecolor','r'); hold on
        m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
        m_pcolor(xx,yy,MonthPPM(:,:,i)); shading flat; hold on;
        m_coast('patch',[.7 .7 .7],'edgecolor','k'); 
        caxis ([0 1000]); c = colorbar; ylabel(c,'Chl (mg Chl/m2/day)');
        tsrt = sprintf('Model %s', mnom(i,:));
        title(tsrt,'FontWeight','bold')
        
        subplot(1,2,2);
        m_proj('stereographic','lat',90,'long',30,'radius',30)
        m_elev('contour',[-500 -500],'edgecolor','r'); hold on
        m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
        m_pcolor(xx,yy,MonthPPS(:,:,i)); shading flat; hold on;
        m_coast('patch',[.7 .7 .7],'edgecolor','k'); 
        caxis ([0 1000]); c = colorbar; ylabel(c,'Chl (mg Chl/m2/day)');
        tsrt = sprintf('Satellite %s', mnom(i,:));
        title(tsrt,'FontWeight','bold')
    end

    
% (iv) bar plot    
    
     % averaging Pan_Arctic
        pos = 1;
        for i= 1:size(MonthPPM,3)
            
            % model
            mod_bar(pos) = nansum(nansum(MonthPPM(:,831:end,i),1),2) / (size(MonthPPM,1) * (size(MonthPPM,2)-831)); % 831 = 65° N. ie. Arctic cicle

            % satellite
            sat_bar(pos) = nansum(nansum(MonthPPS(:,831:end,i),1),2) / (size(MonthPPS,1) * (size(MonthPPS,2)-831));
            sat_bar2(pos) = nansum(nansum(MonthPPS2(:,831:end,i),1),2) / (size(MonthPPS2,1) * (size(MonthPPS2,2)-831));
            pos = pos + 1;        
        end
        
        sat_res = sat_bar - mod_bar;
        sat_res2 = sat_bar2 - mod_bar;
        
        colormap(jet);
        plot_bar = permute([mod_bar; sat_bar; sat_res; sat_res2], [2 1]);
        figure(1); bar(plot_bar); legend('mod PP','sat PP','sat res:vgpm', 'sat res:car'); xlabel('Month (post = May-Sep)'); ylabel ('PP (mgC/m2/d)')
        title('Pan-Arctic PP');
    
        % print -depsc -painters bar_sat_resPP.eps
    