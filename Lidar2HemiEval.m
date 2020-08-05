function Lidar2HemiEval(sfile)

  %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %  % GENERAL DESCRIPTION                                                  
  %  %    This scripts creates synthethic hemispherical images from high     
  %  %   resolution LiDAR point cloud data of forest canopy and computes    
  %  %   sky-view fraction and direct/diffuse/total incoming shortwave      
  %  %   radiation for a given temporal resolution and time period.         
  %  %                                                                      
  %  %    Depending on usage, preparations may be required using L2HE_Prep.m 
  %  %   and L2HEPrep_Settings.m                                          
  %  %                                                                      
  %  % VERSION                                                              
  %  %    This is a first complete version of L2HEval to be referenced as    
  %  %   vs 1.0 | created 2018/08 | last modified on 2020/07/06             
  %  %                                                                      
  %  % AUTHORS                                                              
  %  %   C. Webster (1,2), T. Jonas(1)                                      
  %  %    (1) WSL Institute for Snow and Avalanche Research SLF, Davos, CH  
  %  %   Shortwave radiation function written by T. Jonas for HPEval and    
  %  %     edited here for synthetic hemispheric images.                    
  %  %    (2) University of Edinburgh, School of GeoSciences, Edinburgh, UK 
  %  %                                                                      
  %  % CONTRIBUTING MATERIAL                                                
  %  %   - solar position calculations from NOAA, found online at           
  %  %     www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html              
  %  %     (accessed 2018/06)                                               
  %  %   - function to split solar radiation into direct/diffuse component  
  %  %     as in the factorial snow model FSM 1.0 by R. Essery              
  %  %     available at www.geosci-model-dev.net/8/3867/2015                
  %  %     (accessed 2018/06)                                               
  %  %   - function to convert UTM into WGS84 coordinates, found online at  
  %  %     www.mathworks.com/matlabcentral/fileexchange/44242-utm2lonlat    
  %  %     (accessed 2018/08)                                               
  %  %
  %  % REQUIREMENTS                                                         
  %  %   - file structure and input files are formatted and initiated in    
  %  %     Lidar2HemiEval_Prep.m                                            
  %  %   - lasdata from file exchange required in file path                 
  %  %     (mathworks.com/matlabcentral/fileexchange/48073-lasdata)         
  %  %                                                                      
  %  % SETTINGS                                                             
  %  %    All settings such as paths to I/O data and model parameters are    
  %  %   handled through an external settings file the path/filename of     
  %  %   which is the only input argument to this script                    
  %  %                                                                      
  %  % USAGE                                                                
  %  %   > L2HEval('L2HEval_Settings.m') for single process                  
  %  %  
  %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tic
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% LOAD SETTINGS
  %% > read settings file
  if nargin == 0
      error('No settings file specified')
  elseif nargin == 1
      process = 1;  % to run single process or calibration
  elseif nargin == 2
      process = 3;  % to run calibration process
  else
      error('Incorrect number of input arguments. Possibilities are 1 or 2. Cannot proceed.');
  end
  [setp,setf,~]    = fileparts(sfile);
  addpath(setp);
  if exist(sfile,'file')
    try
        if process == 1
            l2heset       = feval(setf);
        elseif process == 2
            l2heset       = feval(setf,strname);
        end
    catch
      error('Error while reading L2HEval settings file');
    end
  else
    error('L2HEval settings file does not exist');
  end
%   try
%     version        = l2heset.version;
%   catch
%     error('Error while reading L2HEval settings file');
%   end

  %% > update settings if necessary
  % specify earliest version under which settings file does not require update
  if versioncompare('1.0',version) == 1
    % nothing to do / this is just an example line how to update hpeset if older than 1.0
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% CALCULATIONS
  %% > preparatory steps

  % check for compatability in user settings
  if l2heset.sw.calc_type == 1 && l2heset.sw.swr > 0
      error('Error: Incompatible parameter choice between file type and swr calculation option')
  end

  if l2heset.sw.calc_type == 3 && (isempty(l2heset.par.stime) || isempty(l2heset.par.etime))
      error('Error: Start and end time required for l2heset.sw.calc_type = 3');
  end

  if l2heset.sw.tershad == 1 && isempty(l2heset.in.dem)
      error('Error: Cannot calculate terrain mask without dem input')
  end

  if l2heset.sw.tershad < 2 && isempty(l2heset.in.dtm)
      error('Error: Cannot calculate local terrain mask without dtm input')
  end

  if l2heset.sw.svf == 0 && l2heset.sw.swr == 1
      warning('SVf calculation disabled with SWR calculation enabled. Calculating SVf but not saving result')
  end

  if (isfield(l2heset,'calibration') && l2heset.calibration == 1) && process ~= 3 || ...
          (isfield(l2heset,'calibration') && l2heset.calibration == 0) && process == 3
      error('Error: Inconsistency with input arguments and calibration switch in settings file')
  end

  if (isfield(l2heset,'batch') && l2heset.batch == 1) && process ~= 2 || ...
          (isfield(l2heset,'batch') && l2heset.batch == 0) && process == 2
      error('Error: Inconsistency with input arguments and batch processing switch in settings file')
  end

  % create output folder if not existent
  if process == 1 % folder structure for single runs
      if ~exist(l2heset.out.file,'dir')
        mkdir(l2heset.out.file)
      end
      if l2heset.out.save == 1
          shidir = strcat(l2heset.out.file,'/SHIs');
          if ~exist(shidir,'dir')
              mkdir(shidir)
          end
      end
      if l2heset.sw.swr > 0 && l2heset.sw.calc_type < 3
        swrdir = strcat(l2heset.out.file,'/SWR');
        if ~exist(swrdir,'dir')
          mkdir(swrdir)
        end
      end
    outfolder = l2heset.out.file;

  elseif process == 3 % running point size calibration procedure
    outfolder = fullfile(l2heset.out.file,strname);
      if ~exist(outfolder,'dir')
        mkdir(outfolder)
      end
      if l2heset.out.save == 1
        shidir = strcat(outfolder,'/SHIs');
          if ~exist(shidir,'dir')
              mkdir(shidir)
          end
      end
    ratiodir = strcat(outfolder,'/RingRatios');
      if ~exist(ratiodir,'dir'); mkdir(ratiodir); end
  end

  % gather information for "image"
  radius              = l2heset.par.radius;

  % define switch functions based on specified files
  if isempty(l2heset.in.dbh); include_trunks = 0;
  else; include_trunks = 1; end

  if isempty(l2heset.in.ltc); include_branches = 0;
  else; include_branches = 1; end

  % additional changes for calibration procedure to override user inpupt
  if isfield(l2heset,'calibration') && l2heset.calibration == 1
      l2heset.sw.swr     = 0; % turn off calculation of swr if not done already
      l2heset.sw.terrain = 1; % each image needs its own mask
      l2heset.sw.svf     = 1; % must calculate svf
  end

  % create folder for progress
  if ~exist(fullfile(outfolder,'ProgressLastPoint'),'dir')
      mkdir(fullfile(outfolder,'ProgressLastPoint'))
  elseif exist(fullfile(outfolder,'ProgressLastPoint'),'dir')
      rmdir(fullfile(outfolder,'ProgressLastPoint'),'s')
      mkdir(fullfile(outfolder,'ProgressLastPoint'))
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% > read DTM/DSM/DEM/point data
  
      % read specifications of points for which synthethic hemispherical images are to be calculated
      fid = fopen(l2heset.in.pts,'r');
      if fid == -1
        error('Pts file inaccessible');
      end

      try
        if l2heset.sw.calc_type == 1
            formatstr         = '%s%f%f%f';
            fgetl(fid);
            data              = textscan(fid,formatstr,'delimiter','\t');
            fclose(fid);
            pts.id            = data{1};
            pts.x             = data{2};
            pts.y             = data{3};
            pts.camera_height = data{4};
        elseif l2heset.sw.calc_type == 2
            formatstr         = '%s%f%f%f%s%s%f';
            fgetl(fid);
            data              = textscan(fid,formatstr,'delimiter','\t');
            fclose(fid);
            pts.id            = data{1};
            pts.x             = data{2};
            pts.y             = data{3};
            pts.camera_height = data{4};
            for tix = 1:length(pts.id)
                pts.t1(tix)   = datenum(data{5}(tix),'dd.mm.yyyy HH:MM:SS');
                pts.t2(tix)   = datenum(data{6}(tix),'dd.mm.yyyy HH:MM:SS');
            end
            pts.dt            = data{7}/60/24;
        elseif l2heset.sw.calc_type == 3
            formatstr         = '%f%f';
            data              = textscan(fid,formatstr,'delimiter','\t');
            fclose(fid);
            if isempty(data{1})
                error('Error while reading point specification file.\n%s',...
                  'Check l2heset.sw.calc_type and l2heset.in.pts file format are consistent.');
            end
            pts.x             = data{1};
            pts.y             = data{2};
            pts.camera_height = ones(size(pts.x)) .* l2heset.par.camera_height;
            pts.dt            = l2heset.par.time_int/60/24 .* ones(size(pts.x));
            pts.t1            = datenum(l2heset.par.stime,'dd.mm.yyyy HH:MM:SS') .* ones(size(pts.x));
            pts.t2            = datenum(l2heset.par.etime,'dd.mm.yyyy HH:MM:SS') .* ones(size(pts.x));
            clear data
        end
      catch
        try
          fclose(fid);
        catch
        end
        error('Error while reading point specification file.\n%s',...
          'Check l2heset.sw.calc_type and l2heset.in.pts file format are consistent.');
      end

      % arbritrary catches to check data imported properly
      if isempty(pts.x) || range(pts.x) > 5000
        error('Error while reading points specification file: l2heset.sw.calc_type and pts file format not consistent')
      end
      
    % additional changes needed for point size calibration procedure
    if process == 3
      % make a variable with all possible tolerance values
      msmin     = l2heset.par.tol_range(2);
      msmax     = l2heset.par.tol_range(1);
      tolerance = msmin:0.00025:msmax;

      % save the marker size
      dlmwrite(fullfile(outfolder,'Tolerance.txt'),tolerance,'delimiter','\n')

      % now make the pts files
      ptdx  = find(strcmp(pts.id,strname));
      xc = pts.x(ptdx); yc = pts.y(ptdx); camera_height = pts.camera_height(ptdx);
      pts = [];
      pts.id  = string(1:1:length(tolerance));
      pts.x   = xc .* ones(length(tolerance),1)';
      pts.y   = yc .* ones(length(tolerance),1)';
      pts.camera_height = camera_height .* ones(length(tolerance),1)';
    end
    
    % load and extract swr data
    if l2heset.sw.swr > 0 && ~isempty(l2heset.in.swr)
      fid = fopen(l2heset.in.swr,'r');
      if fid == -1
        error('SWR file inaccessible');
      end
      try %reading grid according to ASCII GIS format
        formatstr = '%s%f';
        fgetl(fid);
        data = textscan(fid,formatstr,'delimiter','\t');
        swr.time = datenum(data{1},'dd.mm.yyyy HH:MM:SS');
        swr.data = data{2};
      catch
        try
          fclose(fid);
        catch
        end
        error('Error while reading reference swr data file');
      end
    else
      swr = [];
    end
      
    % read DSM/nDSM data
    try
        [~,~,ext]   = fileparts(l2heset.in.dsm);
        if strcmpi(ext,'.las')                                                   % read DSM/nDSM data assuming las format
          ldat      = lasdata(l2heset.in.dsm,'loadall');
          dsm.x     = ldat.x;                                                    % change lasdata to structure
          dsm.y     = ldat.y;
          dsm.z     = ldat.z;
          if ~isempty(ldat.classification)
              dsm.class = ldat.classification;
              clear ldat
              dsm.x(dsm.class == 2) = [];
              dsm.y(dsm.class == 2) = [];
              dsm.z(dsm.class == 2) = [];
              dsm = rmfield(dsm,'class');
          end
        else
          error('Unknown input file format for dsm data')
        end
    catch
        error('Error reading DSM data')
    end
    
    % read DTM data
    try
      [~,~,ext]     = fileparts(l2heset.in.dtm);
      if strcmpi(ext,'.mat')
        load(l2heset.in.dtm);
      elseif strcmpi(ext,'.txt')                                                                     % read DTM data assuming text format
      dtm                 = load_ascii_grid(l2heset.in.dtm);         % read dtm data assuming ascii/text format
      xdtm                = [dtm.xllcorner:dtm.cellsize:dtm.xllcorner+dtm.cellsize*(dtm.ncols-1)]+dtm.cellsize/2 ; % convert reference from lower left corner of each cell to its center of cell
      ydtm                = [dtm.yllcorner:dtm.cellsize:dtm.yllcorner+dtm.cellsize*(dtm.nrows-1)]+dtm.cellsize/2'; % convert reference from lower left corner of each cell to its center of cell
      [dtm.x,dtm.y]       = meshgrid(xdtm,ydtm);
      dtm.z               = dtm.data(:); dtm.x = dtm.x(:); dtm.y = dtm.y(:);
      dtm.z(abs(dtm.z-dtm.NODATA_value)<eps) = NaN;
      dtm = rmfield(dtm,'data');
      clear xdtm ydtm
          clear xyz
      else
          error('Unknown input file format for dtm data')
      end
    catch
        error('Error reading DTM data')
    end
    
    % read DEM data
    if l2heset.sw.tershad < 2
      try
          if l2heset.sw.tershad == 1
              dem                 = load_ascii_grid(l2heset.in.dem);         % read DEM data assuming ascii/text format
              xdem                = [dem.xllcorner:dem.cellsize:dem.xllcorner+dem.cellsize*(dem.ncols-1)]+dem.cellsize/2 ; % convert reference from lower left corner of each cell to its center of cell
              ydem                = [dem.yllcorner:dem.cellsize:dem.yllcorner+dem.cellsize*(dem.nrows-1)]+dem.cellsize/2'; % convert reference from lower left corner of each cell to its center of cell
              [dem.x,dem.y]       = meshgrid(xdem,ydem);
              dem.z               = dem.data;
              dem.z(abs(dem.z-dem.NODATA_value)<eps) = NaN;
              dem = rmfield(dem,'data');
              clear xdem ydem
          end
      catch
          error('Error reading DEM data')
      end
    end
    
    % read DBH data
    if include_trunks
        try
            [~,~,ext]     = fileparts(l2heset.in.dbh);
            if strcmpi(ext,'.txt')
              dbhdat = importdata(l2heset.in.dbh);
              dbh.x  = dbhdat.data(:,1);
              dbh.y  = dbhdat.data(:,2);
              dbh.h  = dbhdat.data(:,3);
              dbh.r  = (dbhdat.data(:,4)./2); % diameter to radius
              clear dbhdat
            else
              error('Error reading trunk DBH data. Unknown file format.')
            end
        catch
            error('Error reading trunk DBH data')
        end
    end

    % read LTC data
    if include_branches
        try
            [~,~,ext]     = fileparts(l2heset.in.ltc);
            if strcmpi(ext,'.txt')
              ltcdat = importdata(l2heset.in.ltc);
              ltcdat.data(ltcdat.data == -9999) = NaN;
              ltc.px    = ltcdat.data(:,1);
              ltc.py    = ltcdat.data(:,2);
              ltc.pz    = ltcdat.data(:,3);
              ltc.tx    = ltcdat.data(:,4);
              ltc.ty    = ltcdat.data(:,5);
              ltc.hd    = ltcdat.data(:,6);
              ltc.tc    = ltcdat.data(:,7);
              ltc.ba    = ltcdat.data(:,8);
              clear ltcdat
            else
              error('Error reading las tree classification (LTC) data. Unknown file format.')
            end
        catch
            error('Error reading las tree classification (LTC) data.')
        end
    end
    
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% > further compatibility checks

  % 1. check pts coordinates are within the bounds of dsm and dtm
  if min(pts.x) < nanmin(dsm.x) || max(pts.x) > nanmax(dsm.x) || ...
          min(pts.y) < nanmin(dsm.y) || max(pts.y) > nanmax(dsm.y)
      error('Specified input PTS are outside the bounds of the input DSM');
  end
  if min(pts.x) < nanmin(dtm.x) || max(pts.x) > nanmax(dtm.x) || ...
          min(pts.y) < nanmin(dtm.y) || max(pts.y) > nanmax(dtm.y)
      error('Specified input PTS are outside the bounds of the input DTM');
  end

  % 2. check requested modelling period is within the bounds of the input
  %     swr data
  if l2heset.sw.calc_type == 2 && ~isempty(swr) && l2heset.sw.swr == 1
      for ptix = 1:length(pts.x)
          if pts.t1(ptix) < nanmin(swr.time) || pts.t2(ptix) > nanmax(swr.time)
              error('Requested time interval for modelling is outside the bounds of input SWR data');
          end
      end
  elseif l2heset.sw.calc_type == 3  && ~isempty(swr)
      if pts.t1(1) < nanmin(swr.time) || pts.t2(1) > nanmax(swr.time)
          error('Requested time interval for modelling is outside the bounds of input SWR data');
      end
  end

  % 3. check whether different modelling periods are requested for each
  %     point for l2heset.sw.calc_type = 2
  if l2heset.sw.calc_type == 2 && l2heset.sw.swr == 2
      if range(pts.t1) ~= 0 || range(pts.t2) ~= 0 || range(pts.dt) ~= 0
          l2heset.sw.swr = 1; % calculate different solar tracks
          warntxt = sprintf(['Different time intervals given for modelling at each location. \n',...
            'Solar track must be calculated for each individual location.\n',...
            'Over-riding user input for l2heset.sw.calc_type']);
          warning(warntxt); %#ok<SPWRN>
          clear warntxt
      end
  end

    
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% > preprocess DSM/DTM data

  %%% clip dsm to within eval-peri of min/max points
  eval_peri = l2heset.par.eval_peri;
  kpidx = find(dsm.x > (nanmin(pts.x)-eval_peri) & dsm.x < (nanmax(pts.x)+eval_peri) & ...
       dsm.y > (nanmin(pts.y)-eval_peri) & dsm.y < (nanmax(pts.y)+eval_peri));

  dsm.x = dsm.x(kpidx); dsm.y = dsm.y(kpidx); dsm.z = dsm.z(kpidx);
  kpidx = [];

  % clip dtm/dbh/ltc to dsm dimensions for faster processing (dtm go out
  % further)
  dtmbf = eval_peri*3;
  kpidx = find(dtm.x > (nanmin(dsm.x)-dtmbf) & dtm.x < (nanmax(dsm.x)+dtmbf) & ...
      dtm.y > (nanmin(dsm.y)-dtmbf) & dtm.y < (nanmax(dsm.y)+dtmbf));
  dtm.x = dtm.x(kpidx); dtm.y = dtm.y(kpidx); dtm.z = dtm.z(kpidx);
  kpidx = [];
  dtm.e       = zeros(size(dtm.z));

  if include_trunks
      kpidx = find(dbh.x > nanmin(dsm.x) & dbh.x < nanmax(dsm.x) & ...
          dbh.y > nanmin(dsm.y) & dbh.y < nanmax(dsm.y));
      dbh.x = dbh.x(kpidx); dbh.y = dbh.y(kpidx); dbh.h = dbh.h(kpidx); dbh.r = dbh.r(kpidx);
      kpidx = [];
  end

  if include_branches
      kpidx = find(ltc.px > nanmin(dsm.x) & ltc.px < nanmax(dsm.x) & ...
          ltc.py > nanmin(dsm.y) & ltc.py < nanmax(dsm.y));
      ltc.px = ltc.px(kpidx); ltc.py = ltc.py(kpidx); ltc.pz = ltc.pz(kpidx);
      ltc.tx = ltc.tx(kpidx); ltc.ty = ltc.ty(kpidx); ltc.ba = ltc.ba(kpidx);
      ltc.hd = ltc.hd(kpidx);
      kpidx = [];
  end

  %%% interpolate DTM to DSM coordinates (point to point)
  tsi = scatteredInterpolant(dtm.x,dtm.y,dtm.z);
  dsm.e = tsi(dsm.x,dsm.y);

  %%% convert DSM to nDSM (height above terrain) if not done already
  if mode(dsm.z) > 100
      try
          dsm.z = dsm.z - dsm.e;
      catch
          error('Cannot calculate height of DSM points above DTM.\n%s',...
              'Check PTS coordinates are within boundaries of DTM and/or DSM.');
      end
  end

  %%% ignore values below height cutoff hth
  dsm.z(dsm.z <= l2heset.par.hg_cutoff) = NaN;

  %%% determine dtm cellsize if it doesn't exist
  try
    if ~isfield(dtm,'cellsize')
      tempdtm = sort(unique(dtm.x));
      dtm.cellsize = tempdtm(2) - tempdtm(1);
      clear tempdtm
    end
  catch
    error('Error calculating cellsize of dtm')
  end

  % interpolate points to dtm
  if l2heset.sw.calc_type == 3
      pts.z = tsi(pts.x,pts.y);
  end
    
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% > generate extra canopy elements

  if include_trunks
      for tix = 1:length(dbh.x(:,1))
        x1 = dbh.x(tix);
        y1 = dbh.y(tix);
        r  = dbh.r(tix);
        h  = dbh.h(tix);

        [xnew,ynew,znew] = calculate_trunks(x1,y1,r,h,30,0.05);

        n = length(xnew);
        if tix == 1
            tsm.x(1:n,:) = xnew; tsm.y(1:n,:) = ynew; tsm.z(1:n,:) = znew;
        else
            tsm.x(end+1:(end+n),:) = xnew; tsm.y(end+1:(end+n),:) = ynew; tsm.z(end+1:(end+n),:) = znew;
        end

        xnew = []; ynew = []; znew = [];

      end

      clear xnew ynew znew

      % interpolate tree coordinates to dtm (tree heights currently normalised)
      tsm.e = tsi(tsm.x,tsm.y);

      % convert tsm to ntsm if not done already
      if mode(tsm.z) > 100
        tsm.z = tsm.z - tsm.e;
      end

      %%% ignore values below height cutoff
      tsm.z(tsm.z <= l2heset.par.hg_cutoff) = NaN;
  end

  if include_branches
      spacing = 0.1;

      % calclate angular and vertical (z) distance using hdist and ba
      ltc.ad = NaN(size(ltc.hd));
      dx1 = ltc.ba  > 90; ltc.ad(dx1) = ltc.hd(dx1) ./ cosd((90-ltc.ba(dx1)));
      dx2 = ltc.ba  < 90; ltc.ad(dx2) = ltc.hd(dx2) ./ cosd((ltc.ba(dx2)-90));
      dx3 = ltc.ba == 90; ltc.ad(dx3) = ltc.hd(dx3);

      ltc.zd = sqrt(ltc.ad.^2 - ltc.hd.^2);

      % calculate height at which branch intersects trunk
      ltc.tz = ltc.pz+ltc.zd;

      nwmatx = NaN((ceil(max(ltc.ad./spacing))),length(ltc.px));
      nwmaty = NaN((ceil(max(ltc.ad./spacing))),length(ltc.px));
      nwmatz = NaN((ceil(max(ltc.ad./spacing))),length(ltc.px));

      for bidx = 1:length(ltc.px)
          try % try catch in place because sometimes linspace won't work with the rounded npts.
              % if it fails after both floor and ceil are tried, the branch is ignored
              npts = floor(ltc.ad(bidx)/spacing);
              nwmatx(1:npts,bidx) = ltc.px(bidx) + (linspace(0,1,npts)).* (ltc.tx(bidx) - ltc.px(bidx));
              nwmaty(1:npts,bidx) = ltc.py(bidx) + (linspace(0,1,npts)).* (ltc.ty(bidx) - ltc.py(bidx));
              nwmatz(1:npts,bidx) = ltc.pz(bidx) + (linspace(0,1,npts)).* (ltc.tz(bidx) - ltc.pz(bidx));
          catch
               try
                   npts = ceil(ltc.ad(bidx)/spacing);
                   nwmatx(1:npts,bidx) = ltc.px(bidx) + (linspace(0,1,npts)).* (ltc.tx(bidx) - ltc.px(bidx));
                   nwmaty(1:npts,bidx) = ltc.py(bidx) + (linspace(0,1,npts)).* (ltc.ty(bidx) - ltc.py(bidx));
                   nwmatz(1:npts,bidx) = ltc.pz(bidx) + (linspace(0,1,npts)).* ((ltc.pz(bidx)+zd(bidx)) - ltc.pz(bidx));
               catch
                   continue % precautionary measure in case floor/ceil don't work
               end
          end

      end

      ltc.x = nwmatx(:); ltc.y = nwmaty(:); ltc.z = nwmatz(:);
      ltc.x(isnan(ltc.z)) = [];
      ltc.y(isnan(ltc.z)) = [];
      ltc.z(isnan(ltc.z)) = [];
      clear newmatx newmaty newmatz

      % interpolate ltc points to dem
      ltc.e  = tsi(ltc.x,ltc.y);

      % clear up memory
      ltc = rmfield(ltc,{'px','py','pz','tx','ty','ba','hd','ad','zd'});

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% create terrain mask

    %% Create blank image matrix
    [xgrid,ygrid]               = meshgrid(1:radius*2,1:radius*2);
    grid.rad                    = sqrt((xgrid - radius - 1).^2 + (ygrid - radius - 1).^2);
    grid.rad(grid.rad > radius) = NaN;
    grid.tht                    = (grid.rad/radius)*90; % zenith angle in degree, i.e. zenith is tht = 0; horizon is tht = 90
    grid.phi                    = flipud(atan2(ygrid - radius,xgrid - radius)); % azimuth angle

    grid.image = ones(size(grid.phi));
    grid.image(isnan(grid.rad)) = 1;
    grid.coorcrt = [xgrid(:),ygrid(:)];

    grid.coorcrt = ((grid.coorcrt - radius) ./ radius) .* 90;
    grid.coorcrt(isnan(grid.coorcrt)) = NaN;

    clear xgrid
    grid = rmfield(grid,{'tht','phi'});

    if l2heset.sw.tershad < 2 && l2heset.sw.terrain == 2 % calculate DEM data for a single centre point and generate the mask

        xdist = nanmax(pts.x) - nanmin(pts.x);
        ydist = nanmax(pts.y) - nanmin(pts.y);

        if xdist > 50 || ydist > 50
            formatspec = ['Warning: User input to calculate single DEM for all coordinates. \n'...
                'X coordinates are %3.0f m apart and Y coordinates are %3.0f m apart.\n'];
            fprintf(formatspec,xdist,ydist);
        end

        % specify points in centre of grid
        loc_x = mean([max(pts.x),min(pts.x)]);
        loc_y = mean([max(pts.y),min(pts.y)]);

        dem = dem2pol(dem,loc_x,loc_y,pts.camera_height(1),...
            l2heset.par.terrain_perim);

        % put in DEM
        tht_temp = dem.rtht + (linspace(0,1,size(min(dem.rtht):0.5:90,2))).* ...
                      ((ones(size(dem.rtht,1),1) .* 90) - dem.rtht);
        phi_temp = repmat(dem.rphi,1,size(min(dem.rtht):0.5:90,2));

        [tx,ty]          = pol2cart(phi_temp(:),tht_temp(:));

        demcart          = [tx,ty];
        idx              = ismembertol(grid.coorcrt,demcart,0.0085,'ByRows',true);
        ndx              = zeros(size(grid.image(:)));
        ndx(idx>0)       = 1;
        imdx             = logical(reshape(ndx,[radius*2,radius*2]));
        grid.image(imdx) = 0;

        clear dem tht_temp phi_temp

    end



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% calculate solar track -> if one track for all points

  if l2heset.sw.swr == 2 && l2heset.sw.calc_type > 1
    % define time axis for calculations
    if isempty(swr)
        loc_time = (pts.t1:pts.dt:pts.t2)';
    else
        loc_time = swr.time;
    end

    % specify points in centre of grid
    loc_x = mean([max(pts.x),min(pts.x)]);
    loc_y = mean([max(pts.y),min(pts.y)]);

    sol = calc_solar_track(loc_x,loc_y,loc_time,l2heset);

  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% > create matrix for data collection

  flat_SVF       = NaN(length(pts.x),1);
  hemi_SVF   = NaN(length(pts.x),1);
  flat_SVF_N     = NaN(length(pts.x),1);
  hemi_SVF_N = NaN(length(pts.x),1);
  flat_SVF_S     = NaN(length(pts.x),1);
  hemi_SVF_S = NaN(length(pts.x),1);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% > create synthetic hemispherical images

  t = toc;
  progtextinit = ['0. Pre-calculation took ' sprintf('%3.2f',t) ' seconds'];
  dlmwrite(fullfile(outfolder,'ProgressLastPoint',progtextinit),NaN);

  
   for coorix = 1:length(pts.x)
    tic

    %% transfer coordinates relative to pts
    dsm   = convert2pol(dsm,tsi,pts.x(coorix),pts.y(coorix),pts.camera_height(coorix),eval_peri);

    % trunks to polar coordinates
    tsmadd = [];
    if include_trunks
        % find any trunks within 4m of coordinates and make denser than the others
        tdst = 4;
        tidx = find(dbh.x <= pts.x(coorix)+tdst & dbh.x >= pts.x(coorix)-tdst ...
            & dbh.y <= pts.y(coorix)+tdst & dbh.y >= pts.y(coorix)-tdst);
        if size(tidx,1) > 0
            tsmadd = [];
             for tixt = 1:length(tidx)
                hdt         = sqrt(((dbh.x(tidx(tixt))-pts.x(coorix))).^2 + ((dbh.y(tidx(tixt))-pts.y(coorix)).^2));
                if hdt < 2
                    int1 = 150; int2 = 0.005; %int2 = hint
                else
                    int1 = 100; int2 = 0.01;
                end

                [xnew,ynew,znew] = calculate_trunks(dbh.x(tidx(tixt)),dbh.y(tidx(tixt)),...
                                                        dbh.r(tidx(tixt)),dbh.h(tidx(tixt)),int1,int2);

                % remove back side of the trunk
                horzdistt     = sqrt((xnew-pts.x(coorix)).^2 + ((ynew-pts.y(coorix)).^2));
                xnew(horzdistt > (min(horzdistt)+ (dbh.r(tidx(tixt))/100))) = [];
                ynew(horzdistt > (min(horzdistt)+ (dbh.r(tidx(tixt))/100))) = [];
                znew(horzdistt > (min(horzdistt)+ (dbh.r(tidx(tixt))/100))) = [];

                n = length(xnew);
                if tixt == 1
                    tsmadd.x(1:n,:) = xnew; tsmadd.y(1:n,:) = ynew; tsmadd.z(1:n,:) = znew;
                    tsmadd.e(1:n,:) = ones(n,1) .* tsi(dbh.x(tidx(tixt)),dbh.y(tidx(tixt)));
                else
                    tsmadd.x(end+1:(end+n),:) = xnew; tsmadd.y(end+1:(end+n),:) = ynew; tsmadd.z(end+1:(end+n),:) = znew;
                    tsmadd.e(end+1:(end+n),:) = ones(n,1) .* tsi(dbh.x(tidx(tixt)),dbh.y(tidx(tixt)));
                end
                xnew = []; ynew = []; znew = [];
             end

             clear xnew ynew znew

          %%% ignore values below height cutoff hth
          tsmadd.z(tsmadd.z <= l2heset.par.hg_cutoff) = NaN;

          % convert tsm to ntsm if not done already
          if mode(tsmadd.z) > 100
            tsmadd.z = tsmadd.z - tsmadd.e;                                                 % height about terrain
          end

          tsmadd   = convert2pol(tsmadd,tsi,pts.x(coorix),pts.y(coorix),pts.camera_height(coorix),4);
          % remove duplicate x/y coordinates based on radial distance
          [tsmadd.rad,sdx] = sort(tsmadd.rad);
          [tx,ty] = pol2cart(tsmadd.phi,tsmadd.tht);
          tsmaddcart = [tx,ty];
          [~,idc] = unique(roundn(tsmaddcart(sdx,:),-1),'rows');
          tsmadd.phi = tsmadd.phi(idc,:);
          tsmadd.tht = tsmadd.tht(idc,:);
          tsmadd.rad = tsmadd.rad(idc,:);
          clear tx ty tsmaddcart


        end

        % convert rest of tsm to polar coordinates
        tsm = convert2pol(tsm,tsi,pts.x(coorix),pts.y(coorix),pts.camera_height(coorix),0.5*eval_peri);

    end

    %%% do the same for the branch lines
    if include_branches
        ltc = convert2pol(ltc,tsi,pts.x(coorix),pts.y(coorix),pts.camera_height(coorix),...
            l2heset.par.branch_peri*eval_peri);
    end

    if l2heset.sw.tershad < 3 % to further include DTM data

      % transfer DTM to polar coordinates
      dtm = dtm2pol(dtm,tsi,pts.x(coorix),pts.y(coorix),pts.camera_height(coorix),dtmbf);

      if l2heset.sw.tershad < 2 && l2heset.sw.terrain == 1 % to further include DEM data for individual points

        %%% transfer DEM to polar coordinates
      dem = dem2pol(dem,pts.x(coorix),pts.y(coorix),pts.camera_height(coorix),...
          l2heset.par.terrain_perim);
      end
    end

    t=toc;
    if exist('progtext1','var')
           delete(fullfile(outfolder,'ProgressLastPoint',progtext1))
    end
    progtext1 = ['1. Transferring to polar coordinates took ' sprintf('%3.2f',t) ' seconds'];
    dlmwrite(fullfile(outfolder,'ProgressLastPoint',progtext1),NaN);



    %% create synthetic hemispherical image
    tic

    % put in surface features
    if include_trunks && include_branches
        dat.phi = [dsm.phi;tsm.phi;ltc.phi];
        dat.tht = [dsm.tht;tsm.tht;ltc.tht];
        dat.rad = [dsm.rad;tsm.rad;ltc.rad];
    elseif ~include_trunks && ~include_branches
        dat.phi = dsm.phi;
        dat.tht = dsm.tht;
        dat.rad = dsm.rad;
    elseif ~include_trunks && include_branches
        dat.phi = [dsm.phi;ltc.phi];
        dat.tht = [dsm.tht;ltc.tht];
        dat.rad = [dsm.rad;ltc.rad];
    elseif include_trunks && ~include_branches
        dat.phi = [dsm.phi;tsm.phi];
        dat.tht = [dsm.tht;tsm.tht];
        dat.rad = [dsm.rad;tsm.rad];
    end

    if ~isempty(tsmadd)
        dat.phi = [dat.phi;tsmadd.phi];
        dat.tht = [dat.tht;tsmadd.tht];
        dat.rad = [dat.rad;tsmadd.rad];
        tsmadd = [];
    end


    %% specify grid to classify
    image2ev = grid.image;
    
    gridcart   = grid.coorcrt;
    dat.tht(dat.tht > 90) = 90;
    slp = 0;

    % classify image with surface points

    rbins = 0:(eval_peri-0)/5:eval_peri; % get radial (distance) bins

    if process == 1 || process == 2
      told  = l2heset.par.tolerance.*(1:-0.025:0.9);
    elseif process == 3
      told  = tolerance(coorix).*(1:-0.025:0.9);
    end

    [tx,ty] = pol2cart(dat.phi,dat.tht);

    datcart = [tx,ty];

    % remove duplicate x/y coordinates based on radial distance
    [dat.rad,sdx] = sort(dat.rad);
    [datcart,idc] = unique(roundn(datcart(sdx,:),-1),'rows');
    dat.rad = dat.rad(idc,:);

    for rdx = 1:length(rbins)-1 % radial distance bins
     tdx = (dat.rad >= rbins(rdx)  & dat.rad < rbins(rdx+1));

     idx = ismembertol(gridcart,datcart(tdx,:),told(rdx),'ByRows',true);

     ndx              = zeros(size(image2ev(:)));
     ndx(idx>0)       = 1;
     imdx             = logical(reshape(ndx,[radius*2,radius*2]));
     image2ev(imdx)   = 0;
     gridcart(imdx,:) = NaN;
     datcart(tdx,:)   = NaN;
     tdx = [];
    end

    if l2heset.sw.tershad < 3 % to further include DTM data

        tht_temp = dtm.rtht + (linspace(0,1,size(min(dtm.rtht):0.5:90+slp,2))) .* ...
                      ((ones(size(dtm.rtht,1),1) .* 90+slp) - dtm.rtht);
        phi_temp = repmat(dtm.rphi,1,size(min(dtm.rtht):0.5:90+slp,2));

        [tx,ty]          = pol2cart(phi_temp(:),tht_temp(:));

        dtmcart          = [tx,ty];
        idx              = ismembertol(gridcart,dtmcart,0.0085,'ByRows',true);
        ndx              = zeros(size(image2ev(:)));
        ndx(idx>0)       = 1;
        imdx             = logical(reshape(ndx,[radius*2,radius*2]));
        image2ev(imdx)   = 0;

        dtm.phi = []; dtm.tht = []; dtm.rad = []; dtm.rphi = []; dtm.rtht = [];
        clear tht_temp phi_temp

        % to further include DEM data for individual points
        if l2heset.sw.tershad < 2 && l2heset.sw.terrain == 1
            tht_temp = dem.rtht + (linspace(0,1,size(min(dem.rtht):0.5:90+slp,2))).* ...
                          ((ones(size(dem.rtht,1),1) .* 90+slp) - dem.rtht);
            phi_temp = repmat(dem.rphi,1,size(min(dem.rtht):0.5:90+slp,2));

            [tx,ty]          = pol2cart(phi_temp(:),tht_temp(:));

            demcart          = [tx,ty];
            idx              = ismembertol(gridcart,demcart,0.0085,'ByRows',true);
            ndx              = zeros(size(image2ev(:)));
            ndx(idx>0)       = 1;
            imdx             = logical(reshape(ndx,[radius*2,radius*2]));
            image2ev(imdx)   = 0;
          dem.rphi = []; dem.rtht = []; dem.phi = []; dem.tht = []; dem.rad = [];
          clear tht_temp phi_temp

        end

    end

    image2ev(isnan(grid.rad)) = 1;
    
    t=toc;
    if exist('progtext2','var')
           delete(fullfile(outfolder,'ProgressLastPoint',progtext2))
    end
    progtext2 = ['2. Creating image matrix took ' sprintf('%3.2f',t) ' seconds'];
    dlmwrite(fullfile(outfolder,'ProgressLastPoint',progtext2),NaN);

    % make some memory space
    gridcoors = [];
    dat = [];
    dsm.phi = []; dsm.tht = []; dsm.rad = [];
    if include_trunks; tsm.phi = []; tsm.tht = []; tsm.rad = []; end
    if include_branches; ltc.phi = []; ltc.tht = []; ltc.rad = []; end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% save synthetic hemispherical image
    tic

    if l2heset.out.save == 1

        if l2heset.sw.calc_type == 3
            exname = ['SHI_PtID' num2str(coorix,'%04d') '_' sprintf('%06.0f',pts.x(coorix)) '_' ...
                sprintf('%06.0f',pts.y(coorix)) '.mat'];
        elseif l2heset.sw.calc_type < 3
            olabel  = [sprintf('%06.0f',pts.x(coorix)) '_' sprintf('%06.0f',pts.y(coorix))];
            exname  = char(strcat('SHI_',pts.id(coorix),'_',olabel,'.mat'));
        end

        save(fullfile(shidir,exname),'image2ev');

    end

    t=toc;
    if exist('progtext3','var')
          delete(fullfile(outfolder,'ProgressLastPoint',progtext3))
    end
    progtext3 = ['3. Saving image matrix took ' sprintf('%3.2f',t) ' seconds'];
    dlmwrite(fullfile(outfolder,'ProgressLastPoint',progtext3),NaN);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% evaluate synthetic hemispherical image
    tic
    %% calculate SVF
    if l2heset.sw.svf || l2heset.sw.swr > 0

      % prepare image for evaluation
      image2ev(image2ev>0) = 1;

      % create zenith rings
      n                                = 9;
      lens_profile_tht                 = 0:10:90;
      lens_profile_rpix                = 0:1/9:1; % linear profile per 10deg zenith angle
      ring_tht                         = 0:90/n:90;                                       % zenith angle in degree, i.e. zenith is tht = 0; horizon is tht = 90
      ring_radius                      = interp1(lens_profile_tht,lens_profile_rpix*radius,ring_tht);

      % loop through zenith angle rings
      white_to_all_ratio               = nan(1,length(ring_radius)-1);
      white_to_all_ratio_S             = nan(1,length(ring_radius)-1);
      white_to_all_ratio_N             = nan(1,length(ring_radius)-1);
      surface_area_ratio_hemi          = nan(1,length(ring_radius)-1);
      surface_area_ratio_flat          = nan(1,length(ring_radius)-1);
      for rix = 1:length(ring_radius)-1
        % identify pixels within zenith ring
        inner_radius                   = ring_radius(rix);                                     % inner radius of zenith ring in deg
        outer_radius                   = ring_radius(rix+1);                                   % outer radius of zenith ring in deg
        relevant_pix                   = find(grid.rad > inner_radius & grid.rad <= outer_radius);   % pixels within zenith ring
        relevant_pix_N                 = find(grid.rad > inner_radius & grid.rad <= outer_radius & ygrid > radius);
        relevant_pix_S                 = find(grid.rad > inner_radius & grid.rad <= outer_radius & ygrid <= radius);
        % calculate tranmissivity per weight and repsective weights
        white_to_all_ratio(rix)        = sum(image2ev(relevant_pix) == 1) / length(relevant_pix);              % raw transimissivity per ring
        white_to_all_ratio_N(rix)      = sum(image2ev(relevant_pix_N) == 1) / length(relevant_pix_N);              % raw transimissivity per ring
        white_to_all_ratio_S(rix)      = sum(image2ev(relevant_pix_S) == 1) / length(relevant_pix_S);              % raw transimissivity per ring
        surface_area_ratio_hemi(rix)   = cos(ring_tht(rix)/360*2*pi) - cos(ring_tht(rix+1)/360*2*pi);                % surface area*  of zenith ring per surface area*  of the entire hemisphere, * areas on hemisphere
        surface_area_ratio_flat(rix)   = sin(ring_tht(rix+1)/360*2*pi)^2 - sin(ring_tht(rix)/360*2*pi)^2;            % surface area** of zenith ring per surface area** of the entire hemisphere, ** areas projected on horizontal surface  
      end

      % save individual ring ratios for calibration
      if process == 3
        txtname = char(strcat((sprintf('%03.0f',pts.id(coorix))),'_',num2str(tolerance(coorix)),'_RingRatios.txt'));
        dlmwrite(fullfile(ratiodir,txtname),white_to_all_ratio,'delimiter','\t','precision',4);
      end

      % calculate SVF
      flat_SVF(coorix)   = sum(white_to_all_ratio.*surface_area_ratio_flat);  % SVF to reflect perspective of a horizonatal flat uplooking sensor surface (weights zenith rings according to their surface area projected onto a horizonatal flat surface)
      flat_SVF_N(coorix) = sum(white_to_all_ratio_N.*surface_area_ratio_flat);
      flat_SVF_S(coorix) = sum(white_to_all_ratio_S.*surface_area_ratio_flat);
      
      hemi_SVF(coorix)   = sum(white_to_all_ratio.*surface_area_ratio_hemi); % SVF to reflect perspective of hemipherically shaped sensor surface or plant (weights zenith rings according to their surface area on the hemisphere)
      hemi_SVF_N(coorix) = sum(white_to_all_ratio_N.*surface_area_ratio_hemi); % 
      hemi_SVF_S(coorix) = sum(white_to_all_ratio_S.*surface_area_ratio_hemi); % 

    end

    %% calculate SWR
    if l2heset.sw.swr > 0

          % define solar track for individual point
          if l2heset.sw.swr == 1 % creating solar track per point

              % define time axis for calculations
              if isempty(l2heset.in.swr)
                  loc_time = (pts.t1(coorix):pts.dt(coorix):pts.t2(coorix))';
              else
                  loc_time = swr.time;
              end

              sol = calc_solar_track(pts.x(coorix),pts.y(coorix),loc_time,l2heset);

          end


        if l2heset.sw.scatter == 1
            drad                  = l2heset.par.scatter_drad;
        else
            drad                  = hpeset.par.scatter_drad(1);
        end

        image_center          = fliplr(size(image2ev)/2);                   % [xcoor, ycoor]
        trans_for             = zeros(length(loc_time),length(drad));

        keepix                = find(sol.tht <= 90);      
        prad                  = interp1(lens_profile_tht,lens_profile_rpix*radius,sol.tht);

        x                     = image_center(1) + sin(deg2rad(sol.phi(keepix))) .* prad(keepix);
        y                     = image_center(2) + cos(deg2rad(sol.phi(keepix))) .* prad(keepix);
    
        for six = 1:length(x)
            for dix = 1:length(drad) % generates a square solar disc (now working in cartesian coordinates)
                dpix              = drad(dix)/90*radius*sqrt(pi)/2;    % half the side length of a square (in pixels) that has the same area than a circle of radius drad (in degree)
                xmm               = [max(1,round(x(six)-dpix)):min(size(image2ev,2),round(x(six)+dpix))];
                ymm               = [max(1,round(y(six)-dpix)):min(size(image2ev,1),round(y(six)+dpix))];
                switch dix
                    case 1
                        image2ev(ymm,xmm) = min(image2ev(ymm,xmm),0.6);
                        nolp_1          = length(find(image2ev(ymm,xmm) == 0.6));
                        noap_1          = numel(image2ev(ymm,xmm));
                        trans_for(keepix(six),1) = nolp_1/noap_1;
                    case 2
                        image2ev(ymm,xmm) = min(image2ev(ymm,xmm),0.61);
                        nolp_2          = length(find(image2ev(ymm,xmm) >= 0.60 & image2ev(ymm,xmm) <= 0.61));
                        noap_2          = numel(image2ev(ymm,xmm));
                        trans_for(keepix(six),2) = (nolp_2-nolp_1)/(noap_2-noap_1);
                    case 3
                        image2ev(ymm,xmm) = min(image2ev(ymm,xmm),0.62);
                        nolp_3          = length(find(image2ev(ymm,xmm) >= 0.60 & image2ev(ymm,xmm) <= 0.62));
                        noap_3          = numel(image2ev(ymm,xmm));
                        trans_for(keepix(six),3) = (nolp_3-nolp_2)/(noap_3-noap_2);
                end
            end
        end       
        
        trans_for(isnan(trans_for)) = 0;
    
        %% calculate transmissivity for direct shortwave radiation
        tweight = l2heset.par.scatter_dwgh;
        switch length(drad)
            case 3
                trans_for_wgt = (trans_for(:,1)*tweight(1) + trans_for(:,2)*tweight(2) + trans_for(:,3)*tweight(3))/sum(tweight(1:3));
            case 2
                trans_for_wgt = (trans_for(:,1)*tweight(1) + trans_for(:,2)*tweight(2))/sum(tweight(1:2));
            case 1
                trans_for_wgt = trans_for(:,1);
            otherwise
                error('Error while calculating transmissivity for direct swr, unknown option');
        end

        %% calculate shortwave radiation above canopy
        if isempty(l2heset.in.swr)
          % this is potential swr, assuming an atmospheric transmissivity of 1 (which also means that all swr is direct / none is diffuse)
          swr.opn  = max(1367*sol.sin_elev,0);
        else
          % this is swr as measured at the open site
          swr.opn = interp1(swr.time,max(swr.data,0),loc_time);
          % set to 0 when sun is below horizon
          swr.opn(sol.sin_elev<0)=0;
        end
    
        loc_slp = NaN; loc_asp = NaN; % option in the future to include slope angle
        
        %% calculate shortwave radiation below canopy
        % split swr into diffuse/direct radiation
        trans_atm      = swr.opn./max(1367*sol.sin_elev,0);               % atmospheric transmissivity
        dif_frac       = ones(size(loc_time));                             % fraction of diffuse radiation
        fix            = sol.sin_elev > 0;
        dif_frac(fix)  = 0.165;                                            % value from FSM implementation --> could look into accuracy/better options of mid-latitude alpine site
        fix            = sol.sin_elev > 0 & trans_atm < 0.80;
        dif_frac(fix)  = 0.9511 - 0.1604.*trans_atm(fix) + 4.388.*trans_atm(fix).^2 - 16.638.*trans_atm(fix).^3 + 12.336.*trans_atm(fix).^4;
        fix            = sol.sin_elev > 0 & trans_atm < 0.22;
        dif_frac(fix)  = 1 - 0.09.*trans_atm(fix);
        dir_frac       = 1-dif_frac;

        % determine maximum direct swr (if sensor is pointed towards sun)
        mdm            = 3; % based on aspect and slope
        switch mdm
        case 0
          max_dir      = 1367.*ones(trans_atm);                            % using solar constant
        case 1
          max_dir      = 1367.*trans_atm.*dir_frac;                        % using measured swr transmissivity and modelled split into direct and diffuse component
        case 2
          trans_dirmax = 1 - 0.165;                                        % value copied from from above fsm splitting function
          max_dir      = 1367.*exp(log(trans_dirmax)./max(sol.sin_elev,0)); % using de beer's law to acount for length of pathway in terms of zenith angle. The extinction coefficient is based on the maximum transmission of direct swr if sun is at zenith
        case 3
          trans_dirmax = 1 - 0.165;
          max_dir_1    = 1367.*exp(log(trans_dirmax)./max(sol.sin_elev,0));
          max_dir_2    = 1367.*trans_atm.*dir_frac;
          max_dir      = min(max_dir_1,max_dir_2);                         % applying both above criteria 1 and 2
        otherwise
          error('Error while calculating maximum direct swr, unknown option');
        end

        % calculate diffuse swr below canopy, neglecting the distinction of / possible differences between SVF_flat, SVF_tsun, and SVF_incl as the latter two are unavailable (unless you are using synthetic HPs)
        swr.for_dif    = dif_frac.*swr.opn.*flat_SVF(coorix);

        % calculate direct swr below canopy
        cang_1           = max(sin((90-sol.tht)/360*2*pi),0.001);              % cosine of angle between flat surface and position of sun
        cang_2           = max(sin((90-sol.tht)/360*2*pi).*cos(loc_slp/360*2*pi) + cos((90-sol.tht)/360*2*pi).*sin(loc_slp/360*2*pi).*cos((loc_asp-sol.phi)/360*2*pi),0); % cosine of angle between inclined surface and position of sun
        swr.for_dir_flat = min(dir_frac.*swr.opn        ,max_dir).*trans_for_wgt;
        swr.for_dir_incl = min(dir_frac.*swr.opn./cang_1,max_dir).*trans_for_wgt.*cang_2;

        % calculate total swr below canopy
        swr.for_all_flat = swr.for_dif + swr.for_dir_flat;                 % swr for a sensor pointing to the zenith
        swr.for_all_incl = swr.for_dif + swr.for_dir_incl;                 % swr for a sensor normal to local aspect/slope (= aspect/slope of nearest DEM grid cell)

        t=toc; % end of image evaluation and swr calculation timing
   
        %% save output: new file if ~run_grid otherwise append to file

        tic
        if l2heset.sw.calc_type == 3 % running the grids
           if coorix == 1 % create the file

               if ~exist('ilabel','var')
                   ilabel = [sprintf('%6.0f',pts.x(coorix)) '_' sprintf('%6.0f',pts.y(coorix))];
               end

               dlmwrite(fullfile(outfolder,strcat('TimeStamp_',ilabel,'.txt')),loc_time,'delimiter','\t','precision',21);
               dlmwrite(fullfile(outfolder,strcat('Coordinates_',ilabel,'.txt')),[pts.x,pts.y,pts.z],'delimiter','\t','precision',12);
               dlmwrite(fullfile(outfolder,strcat('SWR_direct_flat_calculated_',ilabel,'.txt')),swr.for_dir_flat','delimiter','\t','precision',3);
               dlmwrite(fullfile(outfolder,strcat('SWR_diffuse_calculated_',ilabel,'.txt')),swr.for_dif','delimiter','\t','precision',3);
               dlmwrite(fullfile(outfolder,strcat('SWR_total_flat_calculated_',ilabel,'.txt')),swr.for_all_flat','delimiter','\t','precision',3);
               dlmwrite(fullfile(outfolder,strcat('Forest_Trans_calculated_',ilabel,'.txt')),trans_for_wgt','delimiter','\t','precision',3);
           else % append to file from above
               dlmwrite(fullfile(outfolder,strcat('SWR_direct_flat_calculated_',ilabel,'.txt')),swr.for_dir_flat','-append','delimiter','\t','precision',3);
               dlmwrite(fullfile(outfolder,strcat('SWR_diffuse_calculated_',ilabel,'.txt')),swr.for_dif','-append','delimiter','\t','precision',3);
               dlmwrite(fullfile(outfolder,strcat('SWR_total_flat_calculated_',ilabel,'.txt')),swr.for_all_flat','-append','delimiter','\t','precision',3);
               dlmwrite(fullfile(outfolder,strcat('Forest_Trans_calculated_',ilabel,'.txt')),trans_for_wgt','-append','delimiter','\t','precision',3);
           end
        else % save individual files for different points
          olabel  = [sprintf('%6.0f',pts.x(coorix)) '_' sprintf('%6.0f',pts.y(coorix))];
          exname  = char(strcat(pts.id(coorix),'_SWR_',olabel));
        try
          output.timestamp        = loc_time;
          output.swr.open         = swr.opn;
          output.swr.open_dir     = dir_frac.*swr.opn;
          output.swr.for_all_flat = swr.for_all_flat;
          output.swr.for_dir_flat = swr.for_dir_flat;
          output.swr.for_dif      = swr.for_dif;
          output.swr.for_all_incl = swr.for_all_incl;
          output.swr.for_dir_incl = swr.for_all_incl;
          output.for_tau          = trans_for_wgt;

          save(fullfile(swrdir,[exname '.mat']),'output');

          clear output

        catch
          error('Error while saving swr data to file');
        end

        end

    end
    
    if exist('progtext4','var'); delete(fullfile(outfolder,'ProgressLastPoint',progtext4)); end
    progtext4 = ['4. Evaluating image took ' sprintf('%3.2f',t) ' seconds'];
    dlmwrite(fullfile(outfolder,'ProgressLastPoint',progtext4),NaN);

    if l2heset.sw.svf == 1
      % assemble different SVF calculations into matrix for saving.
      svfvals = [flat_SVF(coorix),flat_SVF_N(coorix),flat_SVF_S(coorix),...
                      hemi_SVF(coorix),hemi_SVF_N(coorix),hemi_SVF_S(coorix)];
      data = [coorix,pts.x(coorix),pts.y(coorix),svfvals];
      if process == 1 || process == 3
          fname = 'SVF_calculated';
      end
      if coorix == 1
          try
             fid = fopen(fullfile(outfolder,[fname '.txt']),'w');
             fprintf(fid,'%s\t','Point_ID');
             fprintf(fid,'%s\t','Easting');
             fprintf(fid,'%s\t','Northing');
             fprintf(fid,'%s\t','SVF_flat');
             fprintf(fid,'%s\t','SVF_flat_north');
             fprintf(fid,'%s\t','SVF_flat_south');
             fprintf(fid,'%s\t','SVF_hemi');
             fprintf(fid,'%s\t','SVF_hemi_north');
             fprintf(fid,'%s\n','SVF_hemi_south');
             fprintf(fid,'%0.0f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n',data);
          catch
             fclose(fid);
          end
          fclose(fid);
      else
          try
             fid = fopen(fullfile(outfolder,[fname '.txt']),'a');
             fprintf(fid,'%0.0f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n',data);
          catch
             fclose(fid);
          end
          fclose(fid);
      end
    end
    
    t=toc;
    if exist('progtext5','var')
          delete(fullfile(outfolder,'ProgressLastPoint',progtext5))
    end
    progtext5 = ['5. Exporting results took ' sprintf('%3.2f',t) ' seconds'];
    dlmwrite(fullfile(outfolder,'ProgressLastPoint',progtext5),NaN);

    % make some memory space
    clear image2ev

    percentdone = floor((coorix / length(pts.x)) * 100);
    if exist('outtext','var'); delete(fullfile(outfolder,outtext)); end
    outtext = ['Processing' sprintf('%3.0f',percentdone) '% ...' num2str(coorix) ' of ' num2str(length(pts.x)) '.txt'];
    dlmwrite(fullfile(outfolder,outtext),NaN);

   end % end loop through points
   
  % turn text files into .mat files
  if l2heset.sw.calc_type == 3 && l2heset.sw.swr > 0
    Coord    = importdata(fullfile(outfolder,strcat('Coordinates_',ilabel,'.txt')));
    TimeS    = importdata(fullfile(outfolder,strcat('TimeStamp_',ilabel,'.txt')));
    SWRdif   = importdata(fullfile(outfolder,strcat('SWR_diffuse_calculated_',ilabel,'.txt')));
    SWRdir   = importdata(fullfile(outfolder,strcat('SWR_direct_flat_calculated_',ilabel,'.txt')));
    SWRtot   = importdata(fullfile(outfolder,strcat('SWR_total_flat_calculated_',ilabel,'.txt')));
    ForTrans = importdata(fullfile(outfolder,strcat('Forest_Trans_calculated_',ilabel,'.txt')));

    fname    = fullfile(outfolder,strcat('SWROutput_',ilabel,'.mat'));

    % make mat
    Coord_row = Coord;
    Time_col  = TimeS;

    save(fname,'Coord_row');
    save(fname,'Time_col','-append');
    save(fname,'SWRdif','-append'); clear SWRdif %#ok<*NASGU>
    save(fname,'SWRdir','-append'); clear SWRdir
    save(fname,'SWRtot','-append'); clear SWRtot
    save(fname,'ForTrans','-append'); clear ForTrans

    % delete the .txt files
    files = dir(fullfile(outfolder,'*.txt'));

    for fidx = 1:length(files)

       if ~contains(files(fidx,1).name,'SVF')
           fname = fullfile(outfolder,files(fidx).name);
           delete(fname)
       end

    end

  end
   
   
end
%% > additional functions
function vm = versioncompare(v1,v2)
  % returns +1 if v1 > v2; returns -1 if v1 < v2; returns 0 if v1 = v2
  try
    vers1          = [0 0 0];
    vstr1          = v1;
    xix            = strfind(vstr1,'.');
    switch length(xix)
    case 0
      vers1(1)     = str2double(vstr1);
    case 1
      vers1(1)     = str2double(vstr1(1:xix(1)-1));
      vers1(2)     = str2double(vstr1(xix(1)+1:end));
    case 2
      vers1(1)     = str2double(vstr1(1:xix(1)-1));
      vers1(2)     = str2double(vstr1(xix(1)+1:xix(2)-1));
      vers1(3)     = str2double(vstr1(xix(2)+1:end));
    otherwise
      error(' ');
    end
    vers2          = [0 0 0];
    vstr2          = v2;
    xix            = strfind(vstr2,'.');
    switch length(xix)
    case 0
      vers2(1)     = str2double(vstr2);
    case 1
      vers2(1)     = str2double(vstr2(1:xix(1)-1));
      vers2(2)     = str2double(vstr2(xix(1)+1:end));
    case 2
      vers2(1)     = str2double(vstr2(1:xix(1)-1));
      vers2(2)     = str2double(vstr2(xix(1)+1:xix(2)-1));
      vers2(3)     = str2double(vstr2(xix(2)+1:end));
    otherwise
      error(' ');
    end
    if vers1(1) > vers2(1)
      vm = +1;
    elseif vers1(1) < vers2(1)
      vm = -1;
    elseif vers1(2) > vers2(2)
      vm = +1;
    elseif vers1(2) < vers2(2)
      vm = -1;
    elseif vers1(3) > vers2(3)
      vm = +1;
    elseif vers1(3) < vers2(3)
      vm = -1;
    else
      vm = 0;
    end
  catch
    vm   = NaN;
  end
end

function [xnew,ynew,znew] = calculate_trunks(x1,y1,r_cm,h,np,hint)

    bh = 1.5; % elevation of tree crown projected to ground + 1.5

    r = r_cm ./ 100; % radius must be in metres to match tree height units

    % lower trunk (below breast height, cylinder shaped)
    th = 0:pi/np:2*pi;
    xunit = (r * cos(th) + x1)';
    yunit = (r * sin(th) + y1)';
    zinc = [0.1:hint:bh];

    xnew1 = repmat(xunit,length(zinc),1);
    ynew1 = repmat(yunit,length(zinc),1);
    znew1 = repmat([0.1:hint:bh],length(yunit),1);
    znew1 = znew1(1:end)';

    % upper trunk - cone shaped
    int = (floor(h) - bh)/hint;
    rchange = fliplr([0.001:((r-0.001)/(int-1)):r]);

    xunit = (rchange(:) .* cos(th) + x1)';
    xnew2  = vertcat(xunit(:));
    yunit = (rchange(:) * sin(th) + y1)';
    ynew2  = vertcat(yunit(:));
    zunit = repmat((bh+hint):hint:floor(h),size(xunit,1),1);
    znew2 = vertcat(zunit(:));

    xnew = [xnew1;xnew2];
    ynew = [ynew1;ynew2];
    znew = [znew1;znew2];

end

function pcd = convert2pol(pcd,tsi,xcoor,ycoor,ch,peri)

  xpc           = pcd.x - xcoor;
  ypc           = pcd.y - ycoor;
  zpc           = pcd.z + pcd.e  - tsi(xcoor,ycoor) - ch;
  horzdistt     = sqrt(((pcd.x-xcoor)).^2 + ((pcd.y-ycoor).^2));
  keepix        = find(horzdistt < peri);
  xpc           = xpc(keepix);
  ypc           = ypc(keepix);
  zpc           = zpc(keepix);

  [pcd.phi,pcd.tht,pcd.rad] = cart2sph(xpc(:),ypc(:),zpc(:));
  pcd.tht           = ((pi/2) - pcd.tht) * (180/pi);
  pcd.tht           = pcd.tht ;
%   pcd.tht(pcd.tht >= 90) = NaN;
  pcd.rad(pcd.rad > peri) = NaN;

  pcd.phi(isnan(pcd.rad)) = []; pcd.tht(isnan(pcd.rad)) = []; pcd.rad(isnan(pcd.rad)) = [];
  pcd.phi(isnan(pcd.tht)) = []; pcd.rad(isnan(pcd.tht)) = []; pcd.tht(isnan(pcd.tht)) = [];


end

function [dtm,tbins,mintht] = dtm2pol(dtm,tsi,xcoor,ycoor,ch,peri)

  xpc           = dtm.x - xcoor;
  ypc           = dtm.y - ycoor;
  zpc           = dtm.z + dtm.e - tsi(xcoor,ycoor) - ch;
  horzdistt     = sqrt(((dtm.x-xcoor)).^2 + ((dtm.y-ycoor).^2));
  keepix        = find(horzdistt < peri);
  xpc           = xpc(keepix);
  ypc           = ypc(keepix);
  zpc           = zpc(keepix);

  [dtm.phi,dtm.tht,dtm.rad] = cart2sph(xpc(:),ypc(:),zpc(:));
  dtm.tht           = (pi/2) - dtm.tht;
  dtm.tht           = dtm.tht * (180/pi);
  dtm.tht(dtm.tht >= 90) = NaN;

  dtm.phi(isnan(dtm.tht)) = []; dtm.rad(isnan(dtm.tht)) = []; dtm.tht(isnan(dtm.tht)) = [];

  %%% calculate horizon lines
  rbins             = [2*dtm.cellsize:sqrt(2)*dtm.cellsize:peri];
  tbins             = [-pi+((pi/360)*3):pi/720:pi-((pi/360)*3)]';
  mintht            = ones(length(tbins),1).*90;                     % initilization of horizon line
  for rbix = length(rbins)-1:-1:1
    try
      fix1 = find(dtm.rad >= rbins(rbix) & dtm.rad < rbins(rbix+1));
      [~,fix2] = sort(dtm.phi(fix1));
      mintht = min(mintht,interp1q(dtm.phi(fix1(fix2)),dtm.tht(fix1(fix2)),tbins));
    catch
      continue
    end
  end

  dtm.rphi       = [-pi:pi/180:pi]';
  dtm.rtht       = interp1([tbins(end)-(2*pi); tbins; tbins(1)+(2*pi)],...
                    mintht([end,1:end,1]),dtm.rphi);

  dtm.rad(dtm.rad > peri) = NaN;
  dtm.phi(isnan(dtm.rad)) = []; dtm.tht(isnan(dtm.rad)) = []; dtm.rad(isnan(dtm.rad)) = [];

end

function [dem,tbins,mintht] = dem2pol(dem,xcoor,ycoor,ch,topo_peri)

  xpc           = dem.x - xcoor;
  ypc           = dem.y - ycoor;
  zpc           = dem.z - interp2(dem.x,dem.y,dem.z,xcoor,ycoor) - ch;
  horzdistdem   = sqrt(((dem.x-xcoor)).^2 + ((dem.y-ycoor).^2));
  keepixdem     = find(horzdistdem < topo_peri);
  xpc           = xpc(keepixdem);
  ypc           = ypc(keepixdem);
  zpc           = zpc(keepixdem);

  [dem.phi,dem.tht,dem.rad] = cart2sph(xpc(:),ypc(:),zpc(:));
  dem.tht                   = (pi/2) - dem.tht;
  dem.tht                   = dem.tht * (180/pi);
  dem.tht(dem.tht >= 90,1)  = 90; % do NOT delete those values but set them to tht = 90 to allow consistent interpolation of horizont line in line 342

  %%% calculate horizon lines
  rbins             = [2*dem.cellsize:sqrt(2)*dem.cellsize:topo_peri];
  tbins             = [-pi+((pi/360)*3):pi/720:pi-((pi/360)*3)]';
  mintht            = ones(length(tbins),1).*90;                     % initilization of horizon line
  for rbix = length(rbins)-1:-1:1
    fix1 = find(dem.rad >= rbins(rbix) & dem.rad < rbins(rbix+1));
    [~,fix2] = sort(dem.phi(fix1));
    mintht = min(mintht,interp1q(dem.phi(fix1(fix2)),dem.tht(fix1(fix2)),tbins));
  end

  for rbix = length(rbins)-1:-1:1
    try
      fix1 = find(dtm.rad >= rbins(rbix) & dtm.rad < rbins(rbix+1));
      [~,fix2] = sort(dtm.tht(fix1));
      mintht = min(mintht,interp1q(dtm.phi(fix1(fix2)),dtm.tht(fix1(fix2)),tbins));
    catch
      continue
    end
  end
  
  dem.rphi = [-pi:pi/360:pi]';
  dem.rtht = interp1([tbins(end)-(2*pi); tbins; tbins(1)+(2*pi)],...
              mintht([end,1:end,1]),dem.rphi);


end

function sol = calc_solar_track(loc_x,loc_y,loc_time,l2heset)
% the below calcuations are based on formulas found at https://www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html (accessed June 2018)

  time_zone = l2heset.set.time_zone;

  % conversion of x/y  to lat/lon
  switch l2heset.par.coor_system
    case 'CH1903'
      xd               = (loc_x - 600000)/1000000;
      yd               = (loc_y - 200000)/1000000;
      lon              = (2.6779094 + 4.728982 .* xd + 0.791484 .* xd .* yd + 0.1306 .* xd .* yd.^2 - 0.0436 .* xd.^3) .* 100 ./ 36;
      lat              = (16.9023892 + 3.238272 .* yd - 0.270978 .* xd.^2 - 0.002528 .* yd.^2 - 0.0447 .* xd.^2 .* yd - 0.0140 .* yd.^3) .* 100 ./ 36;
    case 'CH1903+'
      xd               = (loc_x - 2600000)/1000000;
      yd               = (loc_y - 1200000)/1000000;
      lon              = (2.6779094 + 4.728982 .* xd + 0.791484 .* xd .* yd + 0.1306 .* xd .* yd.^2 - 0.0436 .* xd.^3) .* 100 ./ 36;
      lat              = (16.9023892 + 3.238272 .* yd - 0.270978 .* xd.^2 - 0.002528 .* yd.^2 - 0.0447 .* xd.^2 .* yd - 0.0140 .* yd.^3) .* 100 ./ 36;
    case 'TM35FIN'
      xd               = (loc_x - 2600000)/1000000;
      yd               = (loc_y - 1200000)/1000000;
      lon              = (2.6779094 + 4.728982 .* xd + 0.791484 .* xd .* yd + 0.1306 .* xd .* yd.^2 - 0.0436 .* xd.^3) .* 100 ./ 36;
      lat              = (16.9023892 + 3.238272 .* yd - 0.270978 .* xd.^2 - 0.002528 .* yd.^2 - 0.0447 .* xd.^2 .* yd - 0.0140 .* yd.^3) .* 100 ./ 36;
   case 'UTM'
      zone             = sscanf(l2heset.par.utm_zone,'%f%c');
      hemi             = char(zone(2));
      [lat,lon]        = utm2deg(loc_x,loc_y,zone(1),hemi);
  otherwise
    error('Error while converting x/y to lat/lon, the coordinate system is not yet supported');
  end

  % conversion of time
  time_vec                    = datevec(loc_time);                     % if time stamps represent begin/end of interval, but should represent center of interval -> apply correction here
  day_of_year                 = (datenum(time_vec(:,1),time_vec(:,2),time_vec(:,3)));
  time_of_day                 = (datenum(00,00,00,time_vec(:,4),time_vec(:,5),time_vec(:,6)));
  julian_day                  = day_of_year + 1721058.5 + time_of_day - time_zone/24;
  julian_century              = (julian_day - 2451545) ./ 36525;

  % calculate solar elevation angle
  geom_mean_long_sun_deg      = mod(280.46646 + julian_century .* (36000.76983 + julian_century .* 0.0003032),360);
  geom_mean_anom_sun_deg      = 357.52911 + julian_century .* (35999.05029 - 0.0001537 .* julian_century);
  eccent_earth_orbit          = 0.016708634 - julian_century .* (0.000042037 + 0.0000001267 .* julian_century);
  sun_eq_of_ctr               = sin(deg2rad(geom_mean_anom_sun_deg)) .* (1.914602 - julian_century .* (0.004817 + 0.000014 .* julian_century)) + sin(deg2rad(2 .* geom_mean_anom_sun_deg)) .* ( 0.019993 - 0.000101 .* julian_century) + sin(deg2rad(3 .* geom_mean_anom_sun_deg)) .* 0.000289;
  sun_true_long_deg           = sun_eq_of_ctr + geom_mean_long_sun_deg;
  sun_app_long_deg            = sun_true_long_deg - 0.00569 - 0.00478 .* sin(deg2rad(125.04 - 1934.136 * julian_century));
  mean_obliq_ecliptic_deg     = 23 + (26 + ((21.448 - julian_century .* (46.815 + julian_century .* (0.00059 - julian_century .* 0.001813)))) / 60) / 60;
  obliq_corr_deg              = mean_obliq_ecliptic_deg + 0.00256 .* cos(deg2rad(125.04 - 1934.136 .* julian_century));
  sun_declin_deg              = rad2deg(asin(sin(deg2rad(obliq_corr_deg)) .* sin(deg2rad(sun_app_long_deg))));
  var_y                       = tan(deg2rad(obliq_corr_deg ./ 2)) .* tan(deg2rad(obliq_corr_deg ./ 2));
  eq_of_time_minutes          = 4 * rad2deg(var_y .* sin(2 .* deg2rad(geom_mean_long_sun_deg)) - 2 .* eccent_earth_orbit .* sin(deg2rad(geom_mean_anom_sun_deg)) + 4 .* eccent_earth_orbit .* var_y .* sin(deg2rad(geom_mean_anom_sun_deg)) .* cos(2 .* deg2rad(geom_mean_long_sun_deg)) - 0.5 .* var_y .* var_y .* sin(4 .* deg2rad(geom_mean_long_sun_deg)) - 1.25 .* eccent_earth_orbit .* eccent_earth_orbit .* sin(2 .* deg2rad(geom_mean_anom_sun_deg)));
  true_solar_time_min         = mod(time_of_day .* 1440 + eq_of_time_minutes + 4 .* lon - 60 .* time_zone,1440);
  hour_angle_deg              = NaN(size(true_solar_time_min));
  idx1                        = true_solar_time_min ./4 < 0;
  idx2                        = true_solar_time_min ./4 >= 0;
  hour_angle_deg(idx1)        = true_solar_time_min(idx1) ./4 + 180;
  hour_angle_deg(idx2)        = true_solar_time_min(idx2) ./4 - 180;
  solar_zenith_angle_deg      = rad2deg(acos(sin(deg2rad(lat)) .* sin(deg2rad(sun_declin_deg)) + cos(deg2rad(lat)) .* cos(deg2rad(sun_declin_deg)) .* cos(deg2rad(hour_angle_deg))));
  solar_elev_angle_deg        = 90 - solar_zenith_angle_deg;

  % calculate atmospheric diffraction dependent on solar elevation angle
  approx_atm_refrac_deg       = NaN(size(solar_elev_angle_deg));
  approx_atm_refrac_deg(solar_elev_angle_deg > 85) = 0;
  idx1                        = solar_elev_angle_deg > 5 & solar_elev_angle_deg <= 85;
  approx_atm_refrac_deg(idx1) = (58.1 ./ tan(deg2rad(solar_elev_angle_deg(idx1))) - 0.07 ./ (tan(deg2rad(solar_elev_angle_deg(idx1)))).^3 + 0.000086 ./ (tan(deg2rad(solar_elev_angle_deg(idx1)))).^5) ./ 3600;
  idx2                        = solar_elev_angle_deg > -0.757 & solar_elev_angle_deg <= 5;
  approx_atm_refrac_deg(idx2) = (1735 + solar_elev_angle_deg(idx2) .* (-518.2 + solar_elev_angle_deg(idx2) .* (103.4 + solar_elev_angle_deg(idx2) .* (-12.79 + solar_elev_angle_deg(idx2) .* 0.711)))) ./ 3600;
  idx3                        = solar_elev_angle_deg <= -0.757;
  approx_atm_refrac_deg(idx3) = (-20.772 ./ tan(deg2rad(solar_elev_angle_deg(idx3)))) ./ 3600;
  solar_elev_corr_atm_ref_deg = solar_elev_angle_deg + approx_atm_refrac_deg;

  % calculate solar azimuth angle depending on hour angle
  solar_azimuth_angle         = NaN(size(hour_angle_deg));
  idx1                        = hour_angle_deg > 0;
  solar_azimuth_angle(idx1)   = mod(round((rad2deg(acos(((sin(deg2rad(lat)) .* cos(deg2rad(solar_zenith_angle_deg(idx1)))) - sin(deg2rad(sun_declin_deg(idx1)))) ./ (cos(deg2rad(lat)) .* sin(deg2rad(solar_zenith_angle_deg(idx1)))))) + 180)*100000)/100000,360);
  idx2                        = hour_angle_deg <= 0;
  solar_azimuth_angle(idx2)   = mod(round((540 - rad2deg(acos(((sin(deg2rad(lat)) .* cos(deg2rad(solar_zenith_angle_deg(idx2)))) - sin(deg2rad(sun_declin_deg(idx2)))) ./ (cos(deg2rad(lat)) .* sin(deg2rad(solar_zenith_angle_deg(idx2)))))))*100000)/100000,360);

  % convert solar position to tht/phi
  sol.tht                     = 90 - solar_elev_corr_atm_ref_deg;          % convert to zenith angle in degree, i.e. sun @ zenith is tht = 0; sun @ horizon is tht = 90
  sol.az                      = solar_azimuth_angle;
  sol.phi                     = solar_azimuth_angle;
  sol.phi(sol.phi>=360)       = sol.phi(sol.phi>=360) - 360;
  sol.phi(sol.phi<0)          = sol.phi(sol.phi<0)+360;                            % convert to azimuth angle in degree, so that 0� is North, 90� is East, 180� is South, and 270� is East
  sol.sin_elev                = sin(solar_elev_corr_atm_ref_deg/360*2*pi); % sin of elevation angle

end

function answer = load_ascii_grid(varargin)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % DESCRIPTION:
  %   copy of function loadgrid within MetDataWizard for use outside of MDW
  %
  % IMPLEMENTATION:
  %   by TJ in November-2014 @ SLF Switzerland
  %   last changes 04.11.2014
  %
  % INPUT:
  %   prompts the user to select input file in case of nargin = 0
  %   else varargin{1} = path/file to grid to open
  %   supported formats: .mat / ASCII GIS / Massimiliano's binary format
  %
  % OUTPUT: grid structure (internal MDW format) or error message (char)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  answer = [];
  if nargin == 0
    [file,path] = uigetfile('*.*','Load MetDataWizard grid file');
    if isequal(file,0) || isequal(path,0)
      return;
    else
      filepath = fullfile(path,file);
    end
  else
    filepath = varargin{1};
  end

  fid = fopen(filepath,'r');
  if fid == -1
    answer = 'File inaccessible';
    return;
  end
  try %reading grid according to ASCII GIS format
    answer = 'Error reading grid.ncols';
    grid.ncols = fgets(fid);
    fix = strfind(lower(grid.ncols),'ncols');
    if isempty(fix)
      error(' ');
    end
    hlpstr = [grid.ncols(1:fix-1) grid.ncols(fix+5:end)];
    grid.ncols = str2double(hlpstr);
    answer = 'Error reading grid.nrows';
    grid.nrows = fgets(fid);
    fix = strfind(lower(grid.nrows),'nrows');
    if isempty(fix)
      error(' ');
    end
    hlpstr = [grid.nrows(1:fix-1) grid.nrows(fix+5:end)];
    grid.nrows = str2double(hlpstr);
    answer = 'Error reading grid.xllcorner';
    grid.xllcorner = fgets(fid);
    fix = strfind(lower(grid.xllcorner),'xllcorner');
    if isempty(fix)
      error(' ');
    end
    hlpstr = [grid.xllcorner(1:fix-1) grid.xllcorner(fix+9:end)];
    grid.xllcorner = str2double(hlpstr);
    answer = 'Error reading grid.yllcorner';
    grid.yllcorner = fgets(fid);
    fix = strfind(lower(grid.yllcorner),'yllcorner');
    if isempty(fix)
      error(' ');
    end
    hlpstr = [grid.yllcorner(1:fix-1) grid.yllcorner(fix+9:end)];
    grid.yllcorner = str2double(hlpstr);
    answer = 'Error reading grid.cellsize';
    grid.cellsize = fgets(fid);
    fix = strfind(lower(grid.cellsize),'cellsize');
    if isempty(fix)
      error(' ');
    end
    hlpstr = [grid.cellsize(1:fix-1) grid.cellsize(fix+8:end)];
    grid.cellsize = str2double(hlpstr);
    answer = 'Error reading grid.NODATA_value';   %#ok<*NASGU>
    grid.NODATA_value = fgets(fid);
    fix = strfind(lower(grid.NODATA_value),'nodata_value');
    if isempty(fix)
      error(' ');
    end
    hlpstr = [grid.NODATA_value(1:fix-1) grid.NODATA_value(fix+12:end)];
    grid.NODATA_value = str2double(hlpstr);
    answer = 'Error reading grid.data';
    formatstr = '';
    for cix = 1:grid.ncols
      formatstr = [formatstr '%f']; %#ok<AGROW>
    end
    data = textscan(fid,formatstr,grid.nrows);
    for cix = 1:grid.ncols
      grid.data(:,cix) = data{cix};
    end
    grid.data = flipud(grid.data); %conversion necessary in order to comply with own asciigrid standards
    clear data;
    fclose(fid);
    answer = grid;
  catch
    try
      fclose(fid);
    end
    answer = 'File is not a standard ASCII grid';
  end
end

function  [Lat,Lon] = utm2deg(xx,yy,zone,hemi)
  % Credit: based on utm2lonlat found on Matlab file exchange; accessed on 2018/08/25
  % Copyright (c) 2013, Erwin N. All rights reserved.
  x      = xx(:);
  y      = yy(:);
  sa     = 6378137.000000; sb = 6356752.314245;
  e2     = ( ( ( sa .^ 2 ) - ( sb .^ 2 ) ) .^ 0.5 ) ./ sb;
  e2squared = e2 .^ 2;
  c      = ( sa .^ 2 ) ./ sb;
  switch upper(hemi)
  case 'N'
    X    = x - 500000;
    Y    = y;
  case 'S'
    X    = x - 500000;
    Y    = y - 10000000;
  end
  S      = ( ( zone .* 6 ) - 183 );
  lat    =  Y ./ ( 6366197.724 .* 0.9996 );
  v      = ( c ./ ( ( 1 + ( e2squared .* ( cos(lat) ) .^ 2 ) ) ) .^ 0.5 ) .* 0.9996;
  a      = X ./ v;
  a1     = sin( 2 .* lat );
  a2     = a1 .* ( cos(lat) ) .^ 2;
  j2     = lat + ( a1 ./ 2 );
  j4     = ( ( 3 .* j2 ) + a2 ) ./ 4;
  j6     = ( ( 5 .* j4 ) + ( a2 .* ( cos(lat) ) .^ 2) ) ./ 3;
  alpha  = ( 3 ./ 4 ) .* e2squared;
  beta   = ( 5 ./ 3 ) .* alpha .^ 2;
  gamma  = ( 35 ./ 27 ) .* alpha .^ 3;
  Bm     = 0.9996 .* c .* ( lat - alpha .* j2 + beta .* j4 - gamma .* j6 );
  b      = ( Y - Bm ) ./ v;
  Epsi   = ( ( e2squared .* a.^2 ) ./ 2 ) .* ( cos(lat) ).^ 2;
  Eps    = a .* ( 1 - ( Epsi ./ 3 ) );
  nab    = ( b .* ( 1 - Epsi ) ) + lat;
  senoheps = ( exp(Eps) - exp(-Eps) ) ./ 2;
  Delt   = atan(senoheps ./ (cos(nab) ) );
  TaO    = atan(cos(Delt) .* tan(nab));
  longitude = (Delt .* (180/pi) ) + S;
  latitude = ( lat + ( 1 + e2squared .* (cos(lat).^2) - ( 3/2 ) ...
    .* e2squared .* sin(lat) .* cos(lat) .* ( TaO - lat ) ) ...
    .* ( TaO - lat ) ) .* (180/pi);
  Lat    = latitude;
  Lon    = longitude;
end
