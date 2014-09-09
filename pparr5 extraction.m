% PPARR5 Extraction

% Extract model output at lat,lon,t specified in text file
% PPARR5_GCM_CASES.txt, for PPARR5 intercomparison project.


fileID = fopen('/noc/users/jl2g13/MATLAB/PPARR5/PPARR5_GCM_CASES.txt');
N = 8; % number of columns in text file
C_titles = textscan(fileID,'%s',N);
C_data = textscan(fileID,'%d %d %f %f %d %d %d %d','CollectOutput',1); % produces 3 cell arrays (collect output tries to produce 1 cell array, but gets split where  type (integer, floating-point) changes.
loc = C_data{2}; % text-file station lattitudes and longditudes
year = C_data{3}(:,1); % corresponding year station was sampled
year_day = C_data{3}(:,4); % and year day station was sampled
station = C_data{1}(:,1);
case_no = C_data{1}(:,2);
month_day = C_data{3}(:,2:3); % month and day columns into 1005 x 2 matrix

% titles = char(C_titles{1}); % not sure exactly whether this is a good idea, but it works (cell -> character array conversion)
fclose(fileID); % finished reading in text file to Matlab

pos = 1; % iterator for each location

for i = 1:size(loc,1) % for every station in text file
    obs_loc = loc(i,:);
    
    yr = year(i);
    
    
    % converting eastings to regular lon
    if obs_loc(2) > 180 % coord2index requires lon in range (-180,180), text file gives eastings (0-360)
        obs_loc(2) = (360-obs_loc(2)) * -1;
        fprintf('regular longditude is %d\n', obs_loc(2))
    else
        fprintf('unchanged-format longditude is %d\n', obs_loc(2))
    end
    
    
    % converting true locations to model grid cell locations
    % [lon, lat, lon2, lat2] = coord2index(obs_loc(1),obs_loc(2)); %produces lon2, lat2 which are nearest grid points that are ocean to 'lon' and 'lat'.
    [lon, lat, lon2, lat2, actlon1, actlat1, actlon2, actlat2] = coord2index(obs_loc(1),obs_loc(2)); %produces lon2, lat2 which are nearest grid points that are ocean to 'lon' and 'lat'.
   % fprintf('lon = %d, lat = %d, lon2 = %d, lat2 = %d\n', actlon, actlat, actlon2, actlat2);
    if (lon == lon2) & (lat == lat2) % if statments check whether station coordinates (ship obs) are a land or ocean cell in the model
        % they're the same!
        fprintf('- ocean cell\n');
        column14(pos) = 1;
        actlon = actlon1; actlat = actlat1;
    else
        % they're different!
        fprintf('- land cell\n');
        column14(pos) = 0;
        lon = lon2; % if station coordinates are a land cell in model then take nearest ocean cell
        lat = lat2;
        actlon = actlon2; actlat = actlat2;
    end
    fprintf('grid lon = %d, grid lat = %d\n', lon, lat); 
    fprintf('lon = %d, lat = %d, lon2 = %d, lat2 = %d\n', actlon, actlat, actlon2, actlat2);
    
    
    dist_lat = [actlat2, obs_loc(1)];
    dist_lon = [actlon2, obs_loc(2)]; % note if column14(i) = 1, actlon = actlat2, else want actlat2 anyway.
    
    [dist] = dist_fischer(dist_lat,dist_lon,'km');
    distance(pos) = dist;
    
    % calculating distance between station coordinates (ship obs) and coordinates model data extracted from
    
    
    
    
    
    
    
    
                        % convert lon and lat that model will use to read in 

                        % if column14(i) == 1, mod_lat = lat; mod_lon = lon; end
                        %else, % mod_lat = lattitude taken by model, converted back into a true lattitude. 

                        % syntax for calculating distance between two locations (both in
                        % geographic lat-lon should look like this
                    %     obs_lat = obs_loc(1);
                    %     obs_lon = obs_loc(2);
                    %     
                    %     dist_lat = [mod_lat,obs_lat];
                    %     dist_lon = [mod_lon,obs_lon];
                    %     [dist] = dist_fischer(dist_lat,dist_lon,'km');
                    %     dist_n(pos) = dist;

                    %       --> then replace dist_NaNs variable with dist_n, below
    
    
    yr = year(i);
    yrday = year_day(i);
    
    % read in top 100 m (top 14 levels = 13 boxes?)
    depth_res = 14;
    
    % extracting PP, Chl and ice cover -- for now just for 2005
    
    % Getting list of all files for this year
    dnomP = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/*d05P.nc', yr);
    fnamesP = dir (dnomP);
    
    dnomD = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/*d05D.nc', yr);
    fnamesD = dir (dnomD);
    
    dnomT = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/*d05T.nc', yr);
    fnamesT = dir (dnomT);
        
    % how many files are there
    [num, x] = size(fnamesP);
    
        
    % reading in
    for f = 1:1:num
        
        % extract year day from output file name - take date from file name, convert date to year day, if yrday (text file) is on day or within 5 days before filename then read in variables, else skip this output file
        str = fnamesP(f).name;
        str = str (14:end-7); % hard code
        dash = '/';
        d = strcat(str(1:4),dash,str(5:6),dash,str(7:8));
        v = datevec(d); % datevec calculates year day using leap years, but model doesn't have leap years, for alternative work-around function see: http://www.mathworks.co.uk/matlabcentral/answers/81453-datevec-with-no-leap-year-timeseries
        v0 = v;
        v0(:,2:3) = 1;
        outputday = datenum(v) - datenum(v0) + 1;
%          fprintf('5-day output day = %d\n', outputday) % note output is averaged over five days that include and precede this day
        
        if (yrday > (outputday-5)) && (yrday <= outputday)
%             fprintf('textfile yrday %d\n', yrday)
%             fprintf('Using output day %d\n', outputday)
%             fprintf('Using output file %s\n', fnamesP(f).name)
                     

            % (i) Depth-resolved Chl down to 100 m
            fnameP = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/%s', yr, fnamesP(f).name);
            fprintf('Reading %s\n', fnameP);

            % Non-diatom
            t1 = ncread(fnameP, 'CHN', [lon lat 1 1], [1 1 depth_res 1]); % reading single point at [lon lat]
            nav_lat = ncread(fnameP, 'nav_lat', [lon lat], [1 1]);
            nav_lon = ncread(fnameP, 'nav_lon', [lon lat], [1 1]);
            fprintf('nav_lat = %4.2f, nav_lon = %4.2f\n', nav_lat, nav_lon);  
            t2 = squeeze(t1);

            %chn(:,pos) = t2;
            chn = t2;


            % Diatom
            t1 = ncread(fnameP, 'CHD', [lon lat 1 1], [1 1 depth_res 1]);
            t2 = squeeze(t1);

            % chd(:,pos) = t2;
            chd = t2;

            % (ii) Primary production

            fnameD = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/%s', yr, fnamesD(f).name);
            fprintf('Reading %s\n', fnameD);

            % Non-diatom
            t1 = ncread(fnameD,'PRN3',[lon lat 1 1], [1 1 depth_res 1]); % PRN3 is 3D model field
            t1 = squeeze(t1);

            % ndpp(:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)
            ndpp = t1;
            
            % Diatom
            t1 = ncread(fnameD,'PRD3',[lon lat 1 1], [1 1 depth_res 1]);
            t1 = squeeze(t1);

            % dpp(:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr) 
            dpp = t1;

            % (iii) Sea-ice cover
            fnameT = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/%s', yr, fnamesT(f).name);
            fprintf('Reading %s\n', fnameT);

            t1 = ncread(fnameT,'soicecov',[lon lat 1], [1 1 1]); % sea ice cover is a 2D model field, has no depth dimension
            t1 = squeeze(t1);

            icecover(pos) = t1; % lattotal, lon, no of 5-day means (73/yr)
            

        end
        
        
    end
    
    modchl = chn + chd;
    modppC = (ndpp + dpp) * 12.011 * 6.625; % including Redfieldian conversion
    
    
    % (1) depth-integrated chl and pp
    
    dep_nemo; % get model depth points into workspace
    thk = nemo_thk(1:14);
    
    int_modchl = nansum((chn + chd) .* thk); % integrating over depth (100m)
    int_modpp = nansum((ndpp + dpp) .* thk); % adding diatom and non-diatom primary production components
    int_modppC(pos) = int_modpp * 12.011 * 6.625; % conversion from mmolN (model) to mgC (satellites) is * 14 (g/mol) * 6.625 (assuming Redfieldian)
    
    
    
    % (2) interpolated (1m) chl and pp
    
    % model values, although given at box interfaces, really represent
    % values in the middle of boxes so here they're placed in the middle of
    % boxes and a seasurface (chl at z = 0) is added, then interp from 0-100m as per PPARR5 request. 
    nemo_mid = (nemo_dep(1:64) + nemo_dep(2:65)) / 2; 
    modz(1) = 0; modz(2:15) = nemo_mid(1:14);
    interpz = 0:1:100;
    t1 = modchl(1); t1(2:15) = modchl;
    t2 = modppC(1); t2(2:15) = modppC;

    
    interp_chl = interp1(modz,t1,interpz); % vq = interp1(x,v,xq) where vector x = sample points, v = corresponding values, xq = query points
    interp_ppC = interp1(modz,t2,interpz);
    
        % checking interpolation with a plot
    
%         %chl
%         figure(1);
%         plot(nemo_dep(1:14),modchl,'k'); hold on
%         plot(interp_chl,'y'); view([90 -90]); set(gca,'XDir','reverse')
%         legend('raw grid data','interpolated values'); xlabel('depth (m)'); ylabel('model variable chl')
%         hold off;
% 
%         % pp
%         figure(2);
%         plot(nemo_dep(1:14),modppC,'k'); hold on
%         plot(interp_ppC,'y'); view([90 -90]); set(gca,'XDir','reverse')
%         legend('raw grid data','interpolated values'); xlabel('depth (m)'); ylabel('model variable pp')
%         hold off;


    % AT MOMENT YEARS IN WHICH MODEL ISN'T RUN ARE STILL RETURNING VALUES
        % could be tidier by putting condition in line95 (before 'reading
        % in'): if size(fnamesP.name) = 0 in_modppC etc. = NaN; else
        % execute reading in section.
        
        
        % EDIT (27/5/14), following conversation with axy removed years 1980,
        % 1981 from submitted model data.
        
        % Note this also makes two 1/4 run submissions consistent ('MSA003' run
        % only has output for 1988-2006)
        
        if yr < 1988 | yr > 2006 % '|' = or
            int_modppC(pos) = NaN;
            column14(pos) = NaN;
            distance(pos) = NaN;
            icecover(pos) = NaN;

            interp_chl(:) = NaN;
            interp_ppC(:) = NaN;
        end
        
        
        % original
%     if yr < 1980 | yr > 2006 % '|' = or
%         int_modppC(pos) = NaN;
%         column14(pos) = NaN;
%         % dist_n(pos) = NaN; % for when you've calculated distances
%         icecover(pos) = NaN;
%         
%         interp_chl(:) = NaN;
%         interp_ppC(:) = NaN;
%         
%         
%     elseif (yr > 1981) & (yr < 1988)
%         int_modppC(pos) = NaN;
%         column14(pos) = NaN;
%         % dist_n(pos) = NaN; % for when you've calculated distances
%         icecover(pos) = NaN;
%         
%         interp_chl(:) = NaN;
%         interp_ppC(:) = NaN;      
%     end
%         
        
    
    
    % storing interpolated chlvals for each location in a single matrix 
    chlvals(:,pos)= interp_chl;
    ppvals(:,pos) = interp_ppC;

    pos = pos + 1;
    
    fprintf('iteration %d complete', i);
    
end


% creating chl array to read to text files - adds station and case# columns
chlarray = horzcat(double(station(1:i)),double(case_no(1:i)),chlvals');

% same for pp
pparray = horzcat(double(station(1:i)),double(case_no(1:i)),ppvals');




% write text files

% fileID = fopen('/noc/users/jl2g13/MATLAB/PPARR5/chl.txt','w');
% fprintf(fileID,'%6s %12s\n','station','case#'); % needs 100 more columns
% fprintf(fileID,'%6.2f %12.8f\n',chlarray;
% fclose(fileID);


% creating headers


    % for chl
    [rows cols] = size(chlarray);
    
    for q = 1:(cols-2)
        if q < 10, header = sprintf('Chl   %dm', q-1); % fixes formatting problem when write to text file (fills additional spaces with '<num>' otherwise; because matlab being forced to grow array in loop (not preallocated)?
        elseif q < 100, header = sprintf('Chl  %dm', q-1);
        else, header = sprintf('Chl %dm', q-1); end
        j = length(sprintf(header)); % to satisfy assignment in next line, need to know how long sprintf is.
        chl_headers(q,1:j) = sprintf(header);
    end


    % for pp
    [rows cols] = size(pparray);
    
    for q = 1:(cols-2)
        if q < 10, header = sprintf('PP   %dm', q-1); % fixes formatting problem when write to text file (fills additional spaces with '<num>' otherwise; because matlab being forced to grow array in loop (not preallocated)?
        elseif q < 100, header = sprintf('PP  %dm', q-1);
        else, header = sprintf('PP %dm', q-1); end
        j = length(sprintf(header)); % to satisfy assignment in next line, need to know how long sprintf is.   
        pp_headers(q,1:j) = sprintf(header); 
    end
    

% for chl
fid = fopen('chl.txt','wt');
fprintf(fid, 'Station\tCase#\t');
for i = 1:1:(cols-2), fprintf(fid, '%s\t', chl_headers(i,:)); end; fprintf(fid, '\n'); % reading headers to file
x = repmat('%8.4f\t',1,(cols-3)); % data output format (cols-3): station, case# and \n
x2 = repmat('%d\t',1,2); % for station and case# columns
fprintf(fid,[x2, x,'%8.4f\n'],chlarray'); % reading data to file
fclose(fid);

% % for all chl (i = 1-1005) - TEMP
% fid = fopen('arse_all', 'wt');
% fprintf(fid, 'Station\tCase#\t');
% for i = 1:1:(cols-2), fprintf(fid, '%s\t', chl_headers(i,1:8)); end; fprintf(fid, '\n'); % reading headers to file
% x = repmat('%-.4g\t',1,(cols-1)); % data output format
% fprintf(fid,[x,'%.4g\n'],chlarray'); % reading data to file
% fclose(fid);



% for pp
fid = fopen('pp.txt','wt');
fprintf(fid, 'Station\tCase#\t');
for i = 1:1:(cols-2), fprintf(fid, '%8s\t ', pp_headers(i,:)); end; fprintf(fid, '\n'); % reading headers to file
x = repmat('%8.4f\t',1,(cols-3));
x2 = repmat('%d\t',1,2);
fprintf(fid,[x2,x,'%8.4f\n'],pparray');
fclose(fid);

% % for all pp (i= 1-1005) - TEMP
% fid = fopen('arse2_all','wt');
% fprintf(fid, 'Station\tCase#\t');
% for i = 1:1:(cols-2), fprintf(fid, '%s\t', pp_headers(i,1:7)); end; fprintf(fid, '\n'); % reading headers to file
% x = repmat('%.4g\t',1,(cols-1));
% fprintf(fid,[x,'%.4g\n'],pparray');
% fclose(fid);

    

% % (i) old 'for chl'
% 
% 
% fid = fopen('/noc/users/jl2g13/MATLAB/PPARR5/chl.txt','wt');
% x = repmat('%-.4g\t',1,(cols-1)); % data output format
% x2 = repmat('%s\t',1,(cols-3));  % header output format
% 
% fprintf(fid, 'Station\tCase#\t');
% for i = 1:1:98, fprintf(fid, '%s\t', chl_headers(i,1:7)); end; fprintf(fid, '\n');
% 
% %        fprintf(fid,['%s\t %s\t %s\n'],'Station', 'Case#', 'Interp Chl'); % print station and case# headers only
% fprintf(fid,['%s\t %s\t', x2, '%15s\n'],'Station', 'Case#', chl_headers'); % print headers to file
% %       fprintf(fid,[x_fmt2,'%s\n'],headers')
% fprintf(fid,[x,'%.4g\n'],chlarray'); % print data to file
% fclose(fid);
% 
% % (ii) old 'for pp'
% fid = fopen('/noc/users/jl2g13/MATLAB/PPARR5/pp.txt','wt');
% x = repmat('%.4g\t',1,(cols-1));
% x2 = repmat('%7s\t',1,(cols-3));
% 
% fprintf(fid,['%7s\t %s\t', x2, '%7s\n'],'Station', 'Case#', pp_headers'); % print headers to file
% fprintf(fid,[x,'%.4g\n'],pparray');
% fclose(fid);


% (iii) appending int_ppC and icecover to GCM_cases.txt

% table for appended text file - in array form
missing = (ones(1005,4))*-999;
% dist_NaNs = NaN(1005,1); % TEMP filler until calculated mod-obs distance

i = 1; % temp whilst building and not running for all locations in GCM text file
stitch = [double(station) double(case_no) loc double(year) double(month_day) double(year_day) missing int_modppC' column14' distance' icecover']; % remove indexing once finished build

% for when  i =1. ie when building script
% stitch = [double(station(1:i)) double(case_no(1:i)) loc(1:i,:) double(year(1:i)) double(month_day(1:i,:)) double(year_day(1:i)) missing(1:i,:) int_modppC' column14 dist_NaNs(1:i) icecover']; % remove indexing once finished build



% patching back together original GCM cases with model PP and ice added

    % full table - bad formatting
    fid = fopen('/noc/users/jl2g13/MATLAB/PPARR5/PPARR5_GCM_CASES_edited.txt','wt');
    % titles_fmt = '%s\t %s\t %s\t %15s\t %s\t %s\t %s\t %-8s\t %-9s\t %-8s\t %-11s\t %10s\t %9s\t %8s\t %13s\t %s\n'; % earlier version
%     titles_fmt = '%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\n';
    titles_fmt2 = repmat('%s\t ',1,size(stitch,2)-1);
    data_fmt = '%d\t %d\t %f\t %10f\t %12d\t %3d\t %d\t %6d\t %13d\t %13d\t %13d\t %13d\t %24.4f\t %14d\t %18.3f\t %8.2f\n';
 
%     fprintf(fid,titles_fmt, 'Station','Case#','Latitude[Deg]','Longitude[East]','Year','Month','Day','Year_Day', 'NPP(I=1%)','z_eu(1%)','NPP(I=0.1%)','z_eu(0.1%)','Integrated_NPP','Grid_Info','Distance','Ice_Cover');
    fprintf(fid,[titles_fmt2,'%s\n'], 'Station','Case#','Latitude[Deg]','Longitude[East]','Year','Month','Day','Year_Day', 'NPP(I=1%)','z_eu(1%)','NPP(I=0.1%)','z_eu(0.1%)','Integrated_NPP[mgC/m2/d]','Grid_Info','Distance[km]','Ice_Cover[0-1]');
    
    % append_titles = char('pp','ice cover'); % these two lines an early attempt to redraw titles from input text file 
    % fprintf(fid, titles_fmt, [titles append_titles]);
    
    fprintf(fid, data_fmt, stitch');

    fclose(fid);
    
    
    


% WRITE TABULAR DATA TO A TEXT FILE: SEE EG.3 http://www.mathworks.co.uk/help/matlab/ref/fprintf.html