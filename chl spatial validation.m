% Ardyna-model temporal comparison

% (1) Ardyna

 z = 1:10:200; % 200 m because that's what Ardyna13 uses (Figure 5)
 
 % Prebloom (Feb-Apr)
        prec1 = 0.8356 - (0.0026 * z) + 0.9450 * exp (-(((z - 3.83) / 22.21).^2)); % C1
        prec2 = 0.7272 - (0.0009 * z) + 0.8371 * exp (-(((z - 0.00) / 36.20).^2)); % C2
        prec3 = 0.4542 - (0.0007 * z) + 0.8127 * exp (-(((z - 1.91) / 80.52).^2)); % C3
        prec4 = 0.4751 - (0.0013 * z) + 0.9337 * exp (-(((z - 0.00) / 68.35).^2)); % C4
        Apre = permute([(prec1); (prec2); (prec3); (prec4)], [2,1]);

 % Postbloom (May-Sep)
        postc1 = 0.4908 - (0.0019 * z) + 1.2039 * exp (-(((z - 48.07) / 26.43).^2)); % C1
        postc2 = 0.6087 - (0.0026 * z) + 0.9656 * exp (-(((z - 36.05) / 27.27).^2)); % C2
        postc3 = 0.5461 - (0.0016 * z) + 1.0198 * exp (-(((z - 23.81) / 28.47).^2)); % C3
        postc4 = 0.5093 - (0.0017 * z) + 1.1552 * exp (-(((z - 17.77) / 30.12).^2)); % C4
        Apost = permute([(postc1); (postc2); (postc3); (postc4)], [2,1]);
 
 % Winter (Oct-Dec)
        wintc1 = 1.1696 - (0.0045 * z) + 0.1130 * exp (-(((z - 83.42) / 24.99).^2)); % C1
        wintc2 = 0.6519 - (0.0030 * z) + 0.7873 * exp (-(((z - 2.37) / 63.03).^2)); % C2
        wintc3 = 0.0939 - (0.0001 * z) + 1.4592 * exp (-(((z - 1.34) / 66.32).^2)); % C3
        wintc4 = 0.3126 - (0.0013 * z) + 1.3075 * exp (-(((z - 0.00) / 54.03).^2)); % C4     
        Awinter = permute([(wintc1); (wintc2); (wintc3); (wintc4)], [2,1]);
                
        
% (2) Model 

yr = 2005;
    
    lattotal = double(1021);
    latstart = double(820); % only want N of 60°N - note grid is not linear so this is an approximation of 60°N latt line (ie. grid_lat f(grid_lon))
    lontotal = double(1442);
    depth = double(20); % first 20 grid levels = 200 m (to match Ardyna13)
    pos = 1;
    
    % combined
    for m = 1:12
        if m < 10
            dnomP = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/ORCA025-N201_%d0%d*d05P.nc', yr, yr, m);
        else
            dnomP = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/ORCA025-N201_%d%d*d05P.nc', yr, yr, m);
        end

        fnamesP = dir(dnomP);
        [num,x] = size(fnamesP);

        fprintf('Reading month %d\n', m);

            for f = 1:num % ie for each of the files within that given month - read in modpp
                
                % (1) Chlorophyll
                
                fnameP = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/%s', yr, fnamesP(f).name);
                fprintf('Reading %s\n', fnameP);

                % (i) For non-diatom chlorophyll
                t1 = ncread(fnameP,'CHN',[1 latstart 1 1], [lontotal (lattotal-latstart) depth 1]); % PRN is depth-integrated so a 2D model field                
                t1 (t1 == 0) = NaN; % removing zero values (land)
                chn(:,:,:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)

                % (ii) For diatom chlorophyll
                t1 = ncread(fnameP,'CHD',[1 latstart 1 1], [lontotal (lattotal-latstart) depth 1]); % PRD is depth-integrated so a 2D model field
                t1 (t1 == 0) = NaN; % removing zero values (land)
                chd(:,:,:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)       


                pos = pos + 1;

            end

            % Padding - making arrays for each individual month the same size
            % so that they can be sequentially appended to form a single 4D
            % array. Pad value = NaN
            chd = padarray(chd, [0 0 0 (7-num)], NaN, 'post'); % to be able to append each month to the same array requires that each appended array is the same size as the matrix it's being appended to - 7 is the most files that appear in a given month (December)
            chn = padarray(chn, [0 0 0 (7-num)], NaN, 'post');
           
            
            % assign current month to array that will store all months for a given var
            modchl(:,:,:,:,m) = chn + chd; % t1 is non-diatom pp, t2 is diatom pp

            clearvars chn chd % clear all variables read in within f-loop (otherwise when come to assign modice on second,third etc loops it can't assign if the num of files in the previous month isn't the same as the month it's overwriting with
            pos = 1; % must reset for the next month (loop): modndpp and modpp then overwrite preceeding months matricies (but only if num(monthN) >= num(month(N-1)) - this is why all model vars read in within f loop must be cleared

    end
    
    
% Applying water depth filter

% If water depth is < 50 m then ignore that water column, following
% Ardyna13

                                       % anywhere water depth < 50 m -> NaN
        t1 = modchl(:,:,9,:,:) * 0; % only want filter at top 50 m -> not over whole grid. nb. n=9 is z~54m. 
        t2 = sum(isfinite(t1), 3); % chl = finite in water, = NaN below seafloor
        kbathy = t2; % model water depth in terms of k model levels
        kbathy (kbathy == 0) = NaN;
       % figure; pcolor(kbathy(:,:,1,1,1)); shading flat; colorbar % first time slice: test figure only
            % this gives (water depth > 50 m) = 1, water depth < 50 m (inc
            % land) = 0
         
        
        % rearranging dimensions so that loop below appends to last dim  
        modchl = permute(modchl, [1 2 4 5 3]);
        kbathy = permute(kbathy, [1 2 4 5 3]);
        kbathy = reshape(kbathy,1442,201,7,12,1); % need 5th dimension to multiply by modchl in loop
        
        pos = 1;
        for q = 1:size(modchl,5)
            tempChl(:,:,:,:,pos) = modchl(:,:,:,:,q) .* kbathy; % kbathy = bathy mask
            tempChl (tempChl == 0) = NaN;

            fprintf('Water depth filter: read level %d\n', pos)
            pos = pos + 1;
        end

        tempChl = permute(tempChl, [1 2 5 3 4]); % 1442 201 20 7 12

        

% Binning by season

prebloom = tempChl(:,:,:,:,2:4); % months feb-apr etc
postbloom = tempChl(:,:,:,:,5:9);
winter = tempChl(:,:,:,:,10:12);


allseasons = [{prebloom} {postbloom} {winter}];


for i = 1:size(allseasons,2)

    % Binning by surface chlorophyll concentration
        season = permute(allseasons{i}, [1 2 5 4 3]);
        c1 = season; c2 = season; c3 = season; c4 = season;
        

        % setting c1 filter (< 0.1)
            c1 (c1(:,:,:,:,1) > 0.1) = NaN;         
            mask = c1(:,:,:,:,1);
            mask = ~isnan(mask);

            pos = 1; % resetting because used to read in model Chl
            for q = 1:size(c1,5)
                c1(:,:,:,:,pos) = c1(:,:,:,:,q) .* mask;
                c1 (c1 == 0) = NaN; % be careful with these lines - chl values = 0 in 'legitimate' columns (chl(z=0) < 0.1) will also be set to NaN;  

                pos = pos + 1;
            end

            c1 (c1 == 0) = NaN;
            c1 = permute(c1, [1 2 5 4 3]);
            fprintf('Read c1 \n');


        % setting c2 filter (0.1-0.3)
            c2 (c2(:,:,:,:,1) <= 0.1) = NaN; % assign range in one line
            c2 (c2(:,:,:,:,1) > 0.3) = NaN;
            mask = c2(:,:,:,:,1);
            mask = ~isnan(mask);

            pos = 1;
            for q = 1:size(c2,5)
                c2(:,:,:,:,pos) = c2(:,:,:,:,q) .* mask;
                c2 (c2 == 0) = NaN; 

                pos = pos + 1;
            end

            c2 (c2 == 0) = NaN;
            c2 = permute(c2, [1 2 5 4 3]);
            fprintf('Read c2 \n');


         % setting c3 filter (0.3-0.5)
            c3 (c3(:,:,:,:,1) <= 0.3) = NaN;
            c3 (c3(:,:,:,:,1) > 0.5) = NaN;                    
            mask = c3(:,:,:,:,1);
            mask = ~isnan(mask);

            pos = 1;
            for q = 1:size(c3,5)
                c3(:,:,:,:,pos) = c3(:,:,:,:,q) .* mask;
                c3 (c3 == 0) = NaN;

                pos = pos + 1;
            end

            c3 (c3 == 0) = NaN;
            c3 = permute(c3, [1 2 5 4 3]);
            fprintf('Read c3 \n');

         % setting c4 filter (0.5-0.7)
            c4 (c4(:,:,:,:,1) <= 0.5) = NaN;
            c4 (c4(:,:,:,:,1) > 0.7) = NaN;
            mask = c4(:,:,:,:,1);
            mask = ~isnan(mask);

            pos = 1;
            for q = 1:size(c4,5)
                c4(:,:,:,:,pos) = c4(:,:,:,:,q) .* mask;
                c4 (c4 == 0) = NaN;

                pos = pos + 1;
            end

            c4 (c4 == 0) = NaN;
            c4 = permute(c4, [1 2 5 4 3]);
            fprintf('Read c4 \n');


    seasonbins = [{c1} {c2} {c3} {c4}];
    
    pos = 1; % resetting because used to apply chl filter

    for j = 1:size(seasonbins,2) % seasonbins is 1x4 cell array
        chltplot = nansum(nansum(seasonbins{j},1),2) / (size(seasonbins{j},1) * size(seasonbins{j},2)); % x-y averaging
        chltplot = squeeze(chltplot); % t1(z,t1,t2) = t1(20,7,3)

        chltplot = nansum(nansum(chltplot,2),3) / (size(chltplot,2) * size(chltplot,3)); % t averaging // t1(z) = t1(20)
        chlav = sum(chltplot) / length(chltplot); % average of chl profile to normalise it

        nchl(:,pos) = chltplot / chlav;
        fprintf('Calculated normalised Chl for C%d\n', pos);
        
        pos = pos + 1;
    end
    
    Mchl(:,:,i) = nchl;
   
    fprintf('Finished iteration for i = %d\n', i);
 
end

% (3) Plotting

% method 2: plot ardyna vs model out here

   dep_nemo; % grid to true depth function - for plotting
   z = nemo_dep(1:20); % plot top 20 levels in metres

   
% (i) Model Chl profiles for each time-period (pre,post,winter)

% Checking temp: plotting model curves for c1-c4 on same axes
   figure(1); plot(z,Mchl(:,1,1)); hold on; plot(z,Mchl(:,2,1),'r'); plot(z,Mchl(:,3,1),'g'); plot(z,Mchl(:,4,1),'k'); view(90,90); ylim([0 4]); legend('C1','C2','C3','C4') % c1-c4, prebloom (t=1)
   figure(2); plot(z,Mchl(:,1,2)); hold on; plot(z,Mchl(:,2,2),'r'); plot(z,Mchl(:,3,2),'g'); plot(z,Mchl(:,4,2),'k'); view(90,90); ylim([0 4]); legend('C1','C2','C3','C4') % c1-c4, postbloom (t=2)
   figure(3); plot(z,Mchl(:,1,3)); hold on; plot(z,Mchl(:,2,3),'r'); plot(z,Mchl(:,3,3),'g'); plot(z,Mchl(:,4,3),'k'); view(90,90); ylim([0 4]); legend('C1','C2','C3','C4') % c1-c4, winter (t=3)


% (ii) Model vs Ardyna Chl profiles   
figure(4); % Apre/Apost/Awint(z,c) = (20,4) // MnormChl(z,c,t) = (20,4,3)
% prebloom row
subplot(3,4,1); plot(z,Apre(:,1),'r'); hold on; plot(z,(Mchl(:,1,1))); view(90,90); ylim([0 4]); title('Prebloom C1','FontWeight','bold'); % legend('Ardyna2013','model'); xlabel('Depth (m)'); ylabel('normalised Chl');
subplot(3,4,2); plot(z,Apre(:,2),'r'); hold on; plot(z,(Mchl(:,2,1))); view(90,90); ylim([0 4]); title('Prebloom C2','FontWeight','bold')
subplot(3,4,3); plot(z,Apre(:,3),'r'); hold on; plot(z,(Mchl(:,3,1))); view(90,90); ylim([0 4]); title('Prebloom C3','FontWeight','bold')
subplot(3,4,4); plot(z,Apre(:,4),'r'); hold on; plot(z,(Mchl(:,4,1))); view(90,90); ylim([0 4]); title('Prebloom C4','FontWeight','bold')
% postbloom row
subplot(3,4,5); plot(z,Apost(:,1),'r'); hold on; plot(z,(Mchl(:,1,2))); view(90,90); ylim([0 4]); title('Postbloom C1','FontWeight','bold')
subplot(3,4,6); plot(z,Apost(:,2),'r'); hold on; plot(z,(Mchl(:,2,2))); view(90,90); ylim([0 4]); title('Postbloom C2','FontWeight','bold')
subplot(3,4,7); plot(z,Apost(:,3),'r'); hold on; plot(z,(Mchl(:,3,2))); view(90,90); ylim([0 4]); title('Postbloom C3','FontWeight','bold')
subplot(3,4,8); plot(z,Apost(:,4),'r'); hold on; plot(z,(Mchl(:,4,2))); view(90,90); ylim([0 4]); title('Postbloom C4','FontWeight','bold')
% winter row
subplot(3,4,9); plot(z,Awinter(:,1),'r'); hold on; plot(z,(Mchl(:,1,3))); view(90,90); ylim([0 4]); title('Winter C1','FontWeight','bold')
subplot(3,4,10); plot(z,Awinter(:,2),'r'); hold on; plot(z,(Mchl(:,2,3))); view(90,90); ylim([0 4]); title('Winter C2','FontWeight','bold')
subplot(3,4,11); plot(z,Awinter(:,3),'r'); hold on; plot(z,(Mchl(:,3,3))); view(90,90); ylim([0 4]); title('Winter C3','FontWeight','bold')
subplot(3,4,12); plot(z,Awinter(:,4),'r'); hold on; plot(z,(Mchl(:,4,3))); view(90,90); ylim([0 4]); title('Winter C4','FontWeight','bold')
        

        
        
        
        

% ==================        
% Surplus Code

% (1) Entire Open water period mask

%   % Applying additional ice filter to open water area (ice < 10%)  
%     icemask = modicecov < 0.10; % if true = 1, false = 0
%     % t1 = modchl(:,:,1,:,:); % matching modchl dimensions to icemask to do element-wise multiplication
%     % t1 = squeeze(t1)  
%     surfacemask = modchl(:,:,1,:,:) .* icemask; % keep modchl values that are under ice < 10%
%     Eopen = modchl;
%     
%     pos = 1;
%     %for q = 1:size(Eopen,5)
%     for q = 1:size(Eopen,3)    
%         Eopen(:,:,pos,:,:) = Eopen(:,:,q,:,:) .* surfacemask;
%         % Eopen(:,:,:,:,pos) = Eopen(:,:,:,:,q) .* mask;
%         % Eopen (Eopen == 0) = NaN;
%         pos = pos + 1;
%     end
% 
%     Eopen (Eopen == 0) = NaN;


% ice from model not great so will introduce errors into the analysis that
% aren't part of the accuracy of the chl field itself

% (2) Converting bathymetry from levels to metres
% 
%         % Bathymetry in metres
%         bathy = t2 * 0;
%         for k = 1:1:63
%            t3 = find(t2 == k);
%            bathy(t3) = nemo_dep(k+1);
%         end
%         
%         t4 (t4 == 0) = NaN; % removing land
%         figure; pcolor(t4); shading flat; colorbar;