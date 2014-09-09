% Hovmuller for plotting each station for a single year individually


% Notes
% (1) Toggle between transect by selecting 'coords' vector to initialise
% (2) Matlab struggles to plot more than 3 of the subplots in the same
% figure window



% Transect points for Hovmuller (Google Earth)

    % Bering to Fram
    %  70°46'23.37"N            
    % 169°45'22.52"W
    % 
    %  75°50'5.98"N
    % 167° 9'40.50"W
    % 
    %  79°36'58.30"N
    % 164°44'44.00"W
    % 
    %  83° 4'14.85"N
    % 162°45'16.88"W
    % 
    %  87°57'25.73"N
    % 149°25'23.51"W
    % 
    %  87°24'5.17"N
    %   6°36'12.22"W
    % 
    %  82°44'23.33"N
    %   2°37'48.25"E
    % 
    %  76°54'6.69"N
    %   0°50'48.84"W

    % Beaufort to Barents
    %  71°17'51.32"N
    % 140°14'49.39"W  [140 71] land in the model (too near the land for 1/4
    % degree model)
    % 
    %  76°34'51.25"N
    % 140°10'46.77"W
    % 
    %  82° 9'52.66"N
    % 139°39'51.68"W
    % 
    %  86°30'20.89"N
    % 141°31'17.84"W
    % 
    %  88°43'9.77"N
    %  65°17'6.99"E
    % 
    %  83° 2'27.66"N
    %  43°55'36.09"E
    % 
    %  77°59'4.95"N
    %  39°18'31.13"E
    % 
    %  73° 2'45.82"N
    %  38°53'26.67"E

    % SPOTTED ERROR - ALL W COORDS SHOULD BE -VE FOR COORD2INDEX - ALL
    % HOVMULLERS WEST OF MERIDIAN AND PLOTTED BEFORE THIS DATE WILL ACTUALLY BE HOVMULLERS OF THE MIRROR LOCATION (IN THE MEDRIDIAN) I.E. TAKEN E NOT W LON. (16/5/2014) 
% coords = [169 70; 167 76; 164 79; 162 83; 149 88; 6 87; 2 82; 0 77]; % Bering to Fram
coords = [140 71; 140 76; 139 83; 141 86; 65 88; 44 83; 39 78; 39 73]; % Beaufort to Barents

for i = 1:length(coords) % this loop is so ugly but I can't get Matlab to let me index directly into cell array (creates ?[1x2]? structure rather than ?[a] [d]?)
    current_loc = coords(i,:);
    [lon lat] = coord2index(current_loc(1),current_loc(2));
   
    
    % index for where data should be put
    pos = 1;

    for y = 1:1:1   %% start y=1,count=1,end=1 
      % set up year, i.e. 1990 onwards
      yr = 1994 + y;


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


        [num, x] = size(fnamesP); % how many files are there?

        % Preallocate arrays
        chn = zeros(64,73);
        chd = zeros(64,73);
        dpp = zeros(64,73);
        ndpp = zeros(64,73);
        DIN = zeros(64,73);
        phn = zeros(64,73);
        phd = zeros(64,73);


      for f = 1:1:num


        % (1) Chlorophyll

            % Extracting 
            fnameP = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/%s', yr, fnamesP(f).name);
            fprintf('Reading %s\n', fnameP)

            % Non-diatom
            t1 = ncread(fnameP, 'CHN', [lon lat 1 1], [1 1 64 1]); % reading single point at (414,1010)
            t2 = squeeze(t1);

            chn(:,pos) = t2;

            % Diatom
            t1 = ncread(fnameP, 'CHD', [lon lat 1 1], [1 1 64 1]);
            t2 = squeeze(t1);

            chd(:,pos) = t2;



        % (2) Biomass

            % (i) For non-diatom biomass

            t1 = ncread(fnameP,'PHN',[lon lat 1 1], [1 1 64 1]);
            t1 = squeeze(t1);

            phn(:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)

            % (ii) For diatom biomass
            t1 = ncread(fnameP,'PHD',[lon lat 1 1], [1 1 64 1]);
            t1 = squeeze(t1);

            phd(:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)



        % (3) Primary Production 

            % Extracting D files
            fnameD = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/%s', yr, fnamesD(f).name);
            fprintf('Reading %s\n', fnameD); 

            % (i) For non-diatom production

            t1 = ncread(fnameD,'PRN3',[lon lat 1 1], [1 1 64 1]); % PRN3 is 3D model field
            t1 = squeeze(t1);

            ndpp(:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)

            % (ii) For diatom production
            t1 = ncread(fnameD,'PRD3',[lon lat 1 1], [1 1 64 1]);
            t1 = squeeze(t1);

            dpp(:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)    



         % (4) For DIN        

            % Reading in surface nitrogen model field - DIN
            t1 = ncread(fnameP,'DIN',[lon lat 1 1], [1 1 64 1]);
            t1 = squeeze(t1);

            DIN(:,pos) = t1;     



         % (5) Mixed Layer Depth

            % Extracting D files
            fnameT = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/%s', yr, fnamesT(f).name);
            fprintf('Reading %s\n', fnameT);

            % Reading in mixed layer depth model field
            t1 = ncread(fnameT,'somxl010',[lon lat 1], [1 1 1]); % mixed layer depth is 2D field
            t1 = squeeze(t1);

            mld(pos) = t1; % lattotal, lon, no of 5-day means (73/yr)


         % (6) For Sea Ice Cover

            % Reading in sea ice cover model field
            t1 = ncread(fnameT,'soicecov',[lon lat 1], [1 1 1]); % sea ice cover is a 2D model field, has no depth dimension
            t1 = squeeze(t1);

            icecover(pos) = t1; % lattotal, lon, no of 5-day means (73/yr)


          % (7) For Sea Ice Thickness 

            % Extracting I files
            fnameI = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/%s', yr, fnamesI(f).name);
            fprintf('Reading %s\n', fnameI);

            % Reading in sea ice thickness model field
            t1 = ncread(fnameI,'iicethic',[lon lat 1], [1 1 1]); % sea ice thickness is a 2D model field, has no depth dimension
            t1 = squeeze(t1);

            icethick(pos) = t1;


           pos = pos + 1;

      end
    end

    chl = chn + chd;
    pp = (dpp + ndpp) * 12.011 * 6.625; % converting mmolN/m3 to mgc/m3
    biomass = (phn + phd) * 12.011 * 6.625;

    % this empties Matlab's buffers of file information that can cause trouble
    % down the line (or so we have found)
    clear functions


    % Plotting hovmullers  

          %  Depth grid -> meters, time -> days and adding row and column of
          %  NaNs to each variable because pcolor cuts that last row and column
          %  off.

              dep_nemo; % script that produces nemo_dep var which is depth interfaces of nemo grid boxes
              time = 0:5:(size(chl,2))*5; % changing time axis units from 5day files to days   

              % Adding a row of NaNs z = 65, and column of NaNs t = 74 (two alt
              % ways)

                % row of NaNs
                   z65 = NaN(1,size(chl,2));

                   chlplot = [chl; z65]; % cocantenating an extra depth level to all time slices of chl
                   ppplot = [pp; z65];
                   DINplot = [DIN; z65];
                   biomassplot = [biomass; z65];


                % column of NaNs    
                   t1 = chlplot; t1(:,end+1) = NaN;
                   t2 = ppplot; t2(:,end+1) = NaN;
                   t3 = DINplot; t3(:,end+1) = NaN;
                   t4 = biomassplot; t4(:,end+1) = NaN;

                % Scaling thickness of ice * 10 so that it can be seen on plot,
                % and ice cover to a percent
                 icethick = icethick * 10;
                 icecover = icecover * 100;

                % plotting Carbon to chl ratio
                 CtoChl = t4 ./ t1;



         % just plotting the top 13 boxes (n+1 interfaces) ~ 100 m
           figure(i); odvpal(30); % i is the x-y coord-position index          
             subplot(5,1,1)
             %suptitle('Basin transect, station 7, 1995')
             pcolor(time, nemo_dep(1:14),t1(1:14,:)); shading flat; hold on;
             plot(2.5:5:362.5, mld,'w','LineWidth',2); % overlaying mld. 5day mean files taken at their midpoint
             plot(2.5:5:362.5, -(icethick),'k','LineWidth',2); % plotting icethickness
             plot(2.5:5:362.5, -(icecover),'r','LineWidth',2); % plotting icethickness
             axis([0 365 -100 100]) % manually set y axis to go to -100 to accomodate ice thickness variable on same axis
             caxis auto; cb = colorbar; ylabel(cb, 'Chl (mg Chl/m3)')
             set(gca, 'YDir', 'reverse');
             ylabel(''); xlabel('time (5day means)');
             title ('Chlorophyll','FontWeight','bold');

             subplot(5,1,2)
             pcolor(time, nemo_dep(1:14),t2(1:14,:)); shading flat; hold on;
             plot(2.5:5:362.5, mld,'w','LineWidth',2);
             plot(2.5:5:362.5, -(icethick),'k','LineWidth',2);
             plot(2.5:5:362.5, -(icecover),'r','LineWidth',2); % plotting icethickness
             axis([0 365 -100 100])
             caxis auto; cb = colorbar; ylabel(cb, 'PP (mgC/m3/day)')
             set(gca, 'YDir', 'reverse');
             ylabel('Depth, ice thickness (dm) and ice cover (%)'); xlabel('time (5day means)');
             title ('Primary production','FontWeight','bold');

             subplot(5,1,3)
             pcolor(time, nemo_dep(1:14),(t4(1:14,:))); shading flat; hold on % t4 = biomass
             plot(2.5:5:362.5, mld,'w','LineWidth',2);
             plot(2.5:5:362.5, -(icethick),'k','LineWidth',2);
             plot(2.5:5:362.5, -(icecover),'r','LineWidth',2); % plotting icethickness
             axis([0 365 -100 100])
             caxis auto; cb = colorbar; ylabel(cb, 'Biomass (mgC/m3)')
             set(gca, 'YDir', 'reverse');
             ylabel(''); xlabel('time (5day means)');
             title ('Biomass','FontWeight','bold');

             subplot(5,1,4)
             pcolor(time, nemo_dep(1:14),log(t3(1:14,:))); shading flat; hold on % NB log(DIN)
             plot(2.5:5:362.5, mld,'w','LineWidth',2);
             plot(2.5:5:362.5, -(icethick),'k','LineWidth',2);
             plot(2.5:5:362.5, -(icecover),'r','LineWidth',2); % plotting icethickness
             axis([0 365 -100 100])
             caxis auto; cb = colorbar; ylabel(cb, 'log(DIN) (mmol-N/m3)')
             set(gca, 'YDir', 'reverse');
             ylabel('"'); xlabel('time (5day means)');
             title ('DIN','FontWeight','bold');

             subplot(5,1,5)
             pcolor(time, nemo_dep(1:14),CtoChl(1:14,:)); shading flat; hold on
             plot(2.5:5:362.5, mld,'w','LineWidth',2);
             plot(2.5:5:362.5, -(icethick),'k','LineWidth',2);
             plot(2.5:5:362.5, -(icecover),'r','LineWidth',2); % plotting icethickness
             axis([0 365 -100 100])
             caxis auto; cb = colorbar; ylabel(cb, 'C/Chl (g/g)')
             set(gca, 'YDir', 'reverse');
             ylabel('"'); xlabel('time (5day means)');
             title ('Chl/C','FontWeight','bold');
             


end

            print -depsc -painters basin7b.eps

