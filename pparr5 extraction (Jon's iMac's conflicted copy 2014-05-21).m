% fileID = fopen('/noc/users/jl2g13/MATLAB/PPARR5/PPARR5_GCM_CASES.txt')
% C = textscan(fileID,'%s %s %s %s %s %s %s %s', 'CollectOutput',1) % collect output returns single array rather than array-per-column
% lat_col = C{3}{(2:end)};
% lon_col = C{4}



fileID = fopen('/noc/users/jl2g13/MATLAB/PPARR5/PPARR5_GCM_CASES.txt');
N = 8; % number of columns in text file
C_titles = textscan(fileID,'%s',N); % uses formatspec N times so reads top row only (1x8) - the subsequent call to textscan makes it read from the place it last stopped
C_data = textscan(fileID,'%d %d %f %f %d %d %d %d','CollectOutput',1); % produces 3 cell arrays (collect output tries to produce 1 cell array, but gets split where  type (integer, floating-point) changes.
loc = C_data{2}; % text-file station lattitudes and longditudes
year = C_data{3}(:,1); % corresponding year station was sampled
year_day = C_data{3}(:,4); % and year day station was sampled

station = C_data{1}(:,1);
case_no = C_data{1}(:,2);
fclose(fileID); % finished reading in text file to Matlab


pos = 1; % iterator for each location

for i = 1:size(loc,1) % for every station in text file
    current_loc = loc(i,:);
    
    % converting eastings to regular lon
    if current_loc(2) > 180 % coord2index requires lon in range (-180,180), text file gives eastings (0-360)
        current_loc(2) = (360-current_loc(2)) * -1;
        fprintf('regular longditude is %g\n', current_loc(2))
    else
        fprintf('unchanged-format longditude is %g\n', current_loc(2))
    end
    
    
    % converting true locations to model grid cell locations
    [lon lat] = coord2index(current_loc(1),current_loc(2)); % why does coord2index require comma-separated or cell-array (why can't it just take 1x2 matrix)?.
    
    
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
        
        % extract year day from output file name
        str = fnamesP(f).name;
        str = str (14:end-7); % hard code
        dash = '/';
        d = strcat(str(1:4),dash,str(5:6),dash,str(7:8));
        v = datevec(d);
        v0 = v;
        v0(:,2:3) = 1;
        outputday = datenum(v) - datenum(v0) + 1;
        % fprintf('5-day output day = %d\n', outputday) % note output is averaged over five days that include and precede this day
        
        if (yrday > (outputday-5)) && (yrday <= outputday)
            fprintf('textfile yrday %d\n', yrday)
            fprintf('Using output day %d\n', outputday)
            fprintf('Using output file %s\n', fnamesP(f).name)
                        
        %if   
          % - take date from file name
          % - convert date to year day
          % -   if yrday (text file) is on day or within 5 days before filename
          % then read in variables, else skip this output file
          


            % (i) Depth-resolved Chl down to 100 m
            fnameP = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/%s', yr, fnamesP(f).name);
            fprintf('Reading %s\n', fnameP);

            % Non-diatom
            t1 = ncread(fnameP, 'CHN', [lon lat 1 1], [1 1 depth_res 1]); % reading single point at [lon lat]
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
    modppC = (ndpp + dpp) * 14 * 6.625; % including Redfieldian conversion
    
    
    % (1) depth-integrated chl and pp
    
    int_modchl = nansum(chn + chd); % integrating over depth (100m)
    int_modpp = nansum(ndpp + dpp); % adding diatom and non-diatom primary production components
    int_modppC(pos) = int_modpp * 14 * 6.625; % conversion from mmolN (model) to mgC (satellites) is * 14 (g/mol) * 6.625 (assuming Redfieldian)
    
    fprintf( 'model depth-integrated (100m) \n chl = %g (mgChl/m3) \n pp = %g (mgC/m2/d) \n icecover = %g', int_modchl, int_modppC, icecover)
    
    
    % (2) interpolated (1m) chl and pp
    
    dep_nemo; % get model depth points into workspace
    modz = nemo_dep(1:14);
    interpz = 1:1:100;
    
    interp_chl = interp1(modz,modchl,interpz); % vq = interp1(x,v,xq) where vector x = sample points, v = corresponding values, xq = query points
    interp_ppC = interp1(modz,modppC,interpz);
    
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

    
    % storing interpolated chlvals for each location in a single matrix 
    chlvals(:,pos)= interp_chl;
    ppvals(:,pos) = interp_ppC;

    pos = pos +1;
    
end

% I THINK AT THE MOMENT PPARRAY AND CHLARRAY ARE GIVING ME THE NUMBER OF
% NON-ZERO ELEMENTS IN THAT COLUMN (E.G. SEE PP.TXT FILE), NOT THE RAW
% VALUES OF EACH ELEMENT (20/5/14)

% creating chl array to read to text files
chlvals = permute(chlvals, [2 1]);

% chlvals = chlvals(2,:) % temporary
% chlarray = horzcat(station(2), case_no(2), chlvals(1,:))
chlarray = double(horzcat(station(1:i),case_no(1:i), chlvals)); % int32 -> double. 'save' function only takes double arrays

% same for pp
ppvals = permute(ppvals, [2 1]);
pparray = horzcat(station(1:i),case_no(1:i), ppvals);



% write text files
        % fout = fopen('/noc/users/jl2g13/MATLAB/PPARR5/chl.txt','wt'); % t = tab deliniated i think 
        % N = 100;
        % fprintf(fout,'%7s %5s \n','station','case#', N); % needs 100 more columns - note 5 and 7 would be set by default
        % fprintf(fout,'%6.2f %12.8f\n',chlarray;
        % fclose(fout);
        
        % worked example:
%          A = [1 2 4; 5 6 7];
%          
%          fid = fopen('example.txt','w');
%          [rows cols] = size(A);
%          x = repmat('%d\t',1,(cols-1));
%          fprintf(fid,[x,'%d\n'],A'); % note transpose because fprintf reads down columns not across rows
         
        % now for chlarray
         fid = fopen('/noc/users/jl2g13/MATLAB/PPARR5/chl.txt','w');
         [rows cols] = size(chlarray);
         x_fmt = repmat('%g\t',1,(cols-1)); % replicate format string to be size 1x(cols-1). cols-1 because last col wants \n not \t (see next code line) 
         x2_fmt = repmat('%s\t',1,(cols-1));
         
         pos = 1; % want array of 
         for i = 1:cols-2
            inom = char('Chl at %dm', i)
            t1 = sprintf(inom;
            headers(:,:,pos) = t1
            pos = pos + 1;
         end
         
         tsrt = sprintf('%s', mnom((q-4),:));
         
         headers = ['station' 'case#' headers]
            
         fprintf(fid,[x2_fmt,'%s\n'],headers); % titles
         
         fprintf(fid,[x_fmt,'%g\n'],chlarray');  % data
         fclose(fid);
         
         
         
        % and pparray 
         fid = fopen('/noc/users/jl2g13/MATLAB/PPARR5/pp.txt','w');
         [rows cols] = size(pparray);
         x_fmt = repmat('%g\t',1,(cols-1)); % replicate format string to be size 1x(cols-1). cols-1 because last col wants \n not \t (see next code line) 
         fprintf(fid,[x_fmt,'%g\n'],pparray');
         fclose(fid);
        
        
 
        
        
% to add integrated-PP, ice cover etc to existing text file use fopen with 'a' option (append).        
        


% WRITE TABULAR DATA TO A TEXT FILE: SEE EG.3 http://www.mathworks.co.uk/help/matlab/ref/fprintf.html

% for appending to existing cases file, see:
% http://www.mathworks.co.uk/help/matlab/import_export/writing-to-text-data-files-with-low-level-io.html#br5_kad-1
% (section 2)
