%  satsurfchl plot (Chl in mg/m3) (July 2005)

% to extract tar file: tar -xvf filename.tar
% to unzip extracted files: gunzip *.hdf.gz


pos = 1;
yr = 2005;
% for y = 1:1           % this loop and next line for using multiple years of satellite data
    % yr = 1996 + y     % seawifs data for end1997-mid2010
    
    
    % dnomP = sprintf('/noc/users/jl2g13/MATLAB/seawifsCHL/%d/*d05P.nc', yr); % when you download sat data en-masse organise it in subdirectories by year (i.e. ../MATLAB/seawifsCHL/1990/'8-day_means' etc.)
    dnomP = sprintf('/noc/users/jl2g13/MATLAB/seawifsCHL/%d/*.hdf', yr);
    fnamesP = dir (dnomP);
    
    [num, x] = size(fnamesP);
    
    
 
    % Preallocating
    satsurfchl = zeros(2160, 1080, num);  % lon x lat x num of 8-day averages reading in (note lon and lat both at 1/6°)
    
    
    for f = 1:num
        
      fnameP = sprintf('/noc/users/jl2g13/MATLAB/seawifsCHL/%d/%s', yr, fnamesP(f).name); % for when you start using multiple yrs of satellite data
      % fnameP = fnamesP(f).name;
      fprintf('Reading %s\n', fnameP);
        
        
        % For each file: extracting Chl from hdf file into matlab array
        % (for verbose documentation see HDF_to_array.m)
        fileID = hdfsd('start',fnameP, 'read'); % nb. no quotes around file name if it's already assigned as a string
                                                  % nb will only pick up file if you are in dir MATLAB/seawifsCHL/2005 (being in a parent dir (MATLAB) isn't enough
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
% end

ochl = flipdim(permute(satsurfchl, [2 1 3]), 1);
figure; pcolor(log(ochl(:,:,20))); shading flat; colorbar; title('ochl t = 20')



% Plotting in M_Map

% Calculating Chl over July (slices 24-27 in dim = 3; equals yDay =
% 185?209)
t1 = log10(ochl(:,:,24:27));
t2 = nansum(t1,3);


% creating grid for M_Map
x1 = (-180 + (1/12)):(1/6):(180 - (1/12));
y1 = (-90 + (1/12)):(1/6):(90 - (1/12));
[x2,y2] = meshgrid(x1, y1);


    
%t1 = log10(ochl(:,:,20));
% setting colorbar spacing
yt = [0.01 0.1 1 10];
lyt = log10(yt);


% plotting
figure; clf; odvpal(10);
m_proj('stereographic','lat',90,'long',30,'radius',30)
m_elev('contour',[-500 -500],'edgecolor','r'); hold on
m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
m_pcolor(x2,y2,t2); shading flat; hold on;
m_coast('patch',[.7 .7 .7],'edgecolor','k'); 
caxis ([lxt(1) lxt(end)]); cb = colorbar; set(cb, 'YTick', lyt, 'YTickLabel', yt); 
title('obs chl')    


