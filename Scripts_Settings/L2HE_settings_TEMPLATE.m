function l2heset = L2HE_settings_TEMPLATE

%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %  % GENERAL DESCRIPTION                                                  
  %  %    L2HE_Settings defines setting for L2HEval, which is a 
  %  %   tool to generate synthetic hemispherical images from point cloud
  %  %   data of a forest canopy. Sky-view fraction and below-canopy incoming
  %  %   shortwave radiation can also be calculated from each of these
  %  %   synthetic images. 
  %  %  
  %  %   Created by C.Webster @ WSL/SLF
  %  %                                                                      
  %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% version of settings file
  % current version of settings file, any changes of the output structure 
  % array 'l2heset' (adding / removing a field, or changing the size / type
  % of an existing field) require creating a new version number so that 
  % necessary updates can be carried out upon loading l2heset within L2HEval.
  % The versioning system of HPEval supports up to 3 levels, i.e. 1.4.2
  l2heset.version           = '1.0';
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Process settings 
 
  l2heset.calibration       = 0;
  % set to 1 to enter marker size calibration procedure
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Access to input data
  % Provide all paths/files in absolute terms 

  basefolder                = 'D:\L2HEval_master';
  % Use a basefolder to simplify the below path definitions if required
  
  site                      = 'Laret_Edge';
  % Option to use the site name to simplify the below path definitions if required
  %  or change individually below

  l2heset.in.dsm            = fullfile(basefolder,'Data_Surface/DSM/',strcat(site,'.las'));
  % Path/file to data file that contains DSM/nDSM data.
  % Current functionality supports '.las', 'lasm'. Data in .laz format must
  %   be unzipped externally. 
  
  l2heset.in.dtm            = fullfile(basefolder,'/Data_Terrain/DTM/DTM_Laret_Edge.txt');
  % Path/file to data file that contains DTM (local terrain) data.
  % Current functionality supports '.las', 'lasm', '.txt' amd '.mat'.
  % If wanting to calculate radiation to a tilted surface (see
  %   l2heset.sw.tilt), dtm should contain slope and aspect information.
  % Current functionality for this option supports '.mat' and '.txt'.
  % DTM data is used to calculate a points relative height above the ground
  %   when converting to polar coordinates (in conjunction with the
  %   l2heset.par.camera_height setting). 
  % This field can be left empty if DTM data is unavailable (i.e. 
  %   l2heset.in.dtm = '';) and/or not wishing to include terrain shading
  %   (see and adjust settings below). 
  % If left empty, the ground is assumed flat and
  %   ground height will be taken as the lowest point in the input DSM
  %   data. This is suitable for flat areas but l2heset.par.camera_height
  %   will have to be specified with this in mind i.e. camera height will
  %   be metres above lowest point in lidar cloud. If working in
  %   mountainous terrain a DTM is required, particularly to represent
  %   blocking of local terrain. It is suggested a DTM is created using
  %   ground-classified points from a larger lidar tile and matlab's
  %   meshgrid algorithm. 
  
  l2heset.in.dem            = fullfile(basefolder,'Data_Terrain/DEM/dem_50m_clip.txt');
  % Path/file to data file that contains DEM data

  l2heset.in.pts            = fullfile(basefolder,'Data_Points/Laret_Edge_gridpts.txt');
  % Path\file to text file that contains coordinates for which
  %   synthethic hemispherical images and swr are to be created/calculated.
  % The type of file to be read is to be specified in the switch options 
  %   below as l2heset.sw.calc_type. Files can be in one of three formats: 
  % 1 > User defined .txt file for creating individual images at specified point
  %     locations. This .txt file should have four columns (with headers) 
  %     and include (1) PtID, (2) easting, (3) northing, in same coordinates 
  %     as files specified above, (4) camera height. 
  %     This file type is only applicable for calculating sky-view
  %     fraction. 
  % 2 > If swr is to be calculated, text file is as above but additionally 
  %     contains (with headers) (5) start time, (6) end time (both in the format 
  %     'dd.mm.yyyy HH:MM:SS') and (7) time interval in minutes. Starttime 
  %     and endtime may be different for each individual point. 
  % 3 > a 2-column .txt file containing a list of coordinates (1) X and 
  %     (2) Y (without headers)
  %     These points can either be generated separately or generated in the
  %     L2HE_Prep.m script which has the option to remove all points in the 
  %     grid that lie within the boundaries of tree trunks. 
  %     Note: if this type of text file is used, start and end time and 
  %     time interval settings for swr calculations must be included below. 

  l2heset.in.dbh            = fullfile(basefolder,'Data_Surface\DBH\',strcat(site,'_treeinfo.txt'));
  % Input data for creating tree trunks. 
  % This file should be created using Prep4L2HE.R. 
  % Path/file to data  (.txt format) that contains diameter at breast height
  %   and tree location information calculated from R script as columns 
  %   (1) X, (2) Y, (3) normalised height/m, 
  %   (4) diameter at breast height/cm, (5) crown area in cubic metres. 
  % Leave empty if not including tree trunks in synthetic images (i.e.
  %   l2heset.in.dbh = '';)
  
  l2heset.in.ltc            = fullfile(basefolder,'Data_Surface\LTC\',strcat(site,'_lastreeclass.txt'));
  % Input data for creating branches.
  % Path/file to data that contains classified las data from the precalc.R
  %   script but has (1) X coordinates, (2) Y coordinates, (3) Z coordinates
  %   of points that are within crown diameter of a tree, (4) X  and (5) Y 
  %   coordinates of associated tree centre.
  % Leave empty if not including branches in synthetic images (i.e.
  %   l2heset.in.ltc = '';)
 
  l2heset.in.swr            = '';%'D:\L2HEval_master\Data_SWR\Open_SWR_2019_2min_raw.txt';
  % Path/file to swr data from a reference open site or tower. 
  % Leave empty to calculate potential sub-canopy shortwave radiation 
  %   (i.e. l2heset.in.swr = '';) This setting calculates swr under clear sky
  %   conditions with an atmospheric transmissivity of 1.
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% access to output data
  
  l2heset.out.file          = fullfile(basefolder,'Output_Data',site);
  % Path/file to folder into which output results will be saved. 
  
  l2heset.out.save          = 1;
  % This switch enables/disables the saving of each synthetic hemispheric
  %  image as a .mat matrix. The following values are supported:
  % 0 > Will disable saving. Using this option when calculating large
  %     numbers of images across a grid will save CPU time and disk space.
  % 1 > Saves individual hemispheric images following generation. If this
  %     option is selected, a folder named 'SHIs' will be created in the
  %     l2heset.out.file directory. 
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% calculation switches
  
  l2heset.sw.calc_type      = 3;
  % This switch changes the initial processing of the input points file. 
  % The following values are supported  
  % 1 > To be used when only synthetic images are to be created at each 
  %     individual user-specified point. This option can be used with 
  %     l2heset.sw.svf either enabled or disabled but is not compatable when
  %     l2heset.sw.swr is enabled. 
  % 2 > To be used when generating synthetic images at user-specified
  %     locations and calculating sub-canopy shortwave radiation. 
  % 3 > To be used when generating synthetic images at gridded spatial 
  %     scales. This option works with either l2heset.sw.svf or 
  %     l2heset.sw.swr enabled or disabled. 
  %     Note: if this option is selected, start and end time and time 
  %     interval settings for swr calculations must be included below. 
    
  l2heset.sw.terrain        = 1;
  % This switch enables/disables masking of surrounding terrain using the dem.
  % The following values are supported
  % 1 > To enable masking for every individual image. 
  % 2 > To enable a single mask for all images. This option calculates a
  %     terrain mask for the centre point of all input coordinates (grid or
  %     individual points). This will speed up calculations.
  %     It is recommended to use this option when calculating at gridded 
  %     spatial scales (i.e. l2heset.sw.calc_type = 3). 
  % To disable masking change l2heset.sw.image_type below. 
  
  l2heset.sw.tershad    = 1;                                         
  % 1 > output image includes DSM/DTM/DEM
  % 2 > output image includes DSM/DTM
  % 3 > output image includes DSM
  % Note: if DEM data is unavailable or calculation of terrain mask is
  %   disabled, this parameter must be greater than 1. 
  
  l2heset.sw.svf            = 1;
  % This switch enables/disables calculation of sky view fraction. 
  % The following values are supported
  % 0 > will disable calculation of svf. Note that sky view fraction is 
  %     needed for the calculating of sub-canopy swr. If hpeset.sw.swr = 1
  %     svf will be calculated irrespective of this option, but in this 
  %     case no output file will be generated for svf data
  % 1 > will enable calculation of svf  

  l2heset.sw.swr            = 1;
  % This switch enables/disables calculation of below-canopy shortwave
  %   radation. The following values are supported
  % 0 > will disable calculation of swr
  % 1 > will enable calculation of swr, this is the recommended choice.
  % 2 > will enable calculation of swr, but will calculate the position of
  %     the sun once for the centre of all input coordinates. This solar track
  %     will be applied to all images in the calculation. This should be
  %     implemented if calculating across a grid (i.e. l2heset.sw.calc_type
  %     = 3). 
  % Note: if l2heset.sw.calc_type = 2 and l2heset.sw.swr = 2, the swr switch
  %   will only be implemented if the start and end times for the modelling
  %   period are the same for all points. If input modelling periods are
  %   different, the script will automatically calculate different solar
  %   tracks for each point. 
  
  l2heset.sw.scatter        = 1;
  % This switch enables/disables simulating the effect of scattered direct 
  %   swr due to haze, multiple reflections or refration within the canopy. 
  % The following values are supported
  % 0 > will disable simulating scattered direct swr   
  % 1 > will enable simulating scattered direct swr. Note that this option
  %     will entail increased calucation times  
  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% calculation parameters - location 
  
  l2heset.set.time_zone     = 1;
  % Zone of any time data in offset hours relative to GMT, e.g. select +1
  %   if providng swr data from a data logger with CET.
  
  l2heset.par.coor_system   = 'CH1903';
  % Coordinate system used in data / info files. While it is primarily 
  %   important that all coordinates specified (DEM, camera position) feature
  %   the same coordinate system, these need to be converted into lat/lon 
  %   (WGS84) if calculation of swr is enabled. Only for this purpose the 
  %   coordinate system needs to be specified. 
  % Supported option currently include:
  % 'CH1903'  : Swissgrid LV03
  % 'CH1903+' : Swissgrid LV95
  % 'UTM'     : please note that in this case the utm zone must be specified
  
  l2heset.par.utm_zone      = '32N';
  % Zone of utm coordinate system, using either the letter N or S as a 
  %   specifier, whether or not the coordinates are given for the northern
  %   or southern hemisphere, i.e. the letter shall not denote latitude bands
  %   also sometimes used for specifying UTM coordinates
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% calculation parameters - lidar/point cloud
 
  l2heset.par.terrain_perim = 10000;
  % This parameter defines a horizontal perimeter around the camera
  %   position within which the DEM is considered when calculating the
  %   terrain mask. The perimeter is to be given in the same units as
  %   the DEM (e.g meters in the case of UTM coordinates)
  
  l2heset.par.eval_peri     = 100;
  % This parameter defines a perimeter around every evaluation point within
  %   which DTM/DSM is considered when calculating a synthethic hemispherical
  %   image. The evaluation perimeter is defined as the radial direction from 
  %   a point, to be given in the same units as the DTM/DSM data
  %   (m in case of the demonstration dataset). Small evaluation perimeter
  %   will allow for quick calculation, but decrease the accuracy of the 
  %   hemispherical images at high zenith angles (close to to the horizon). 
  % Large evaluation perimeter on the other hand will allow more accurate 
  %   represenation of far distance canopy elements in the hemispherical 
  %   images but require more CPU time.
  % By default, new trunk points are included up to half this distance and
  %   DTM points are included up to twice this distance. 
  %
  
  l2heset.par.branch_peri   = 0.5;
  % This parameter defines the percentage distance of the above evalution
  %   perimeter within which branches will be plotted. Beyond this
  %   distance, the trunk and original lidar points are deemed dense enough
  %   to represent distant canopy. 
  % Value should be between 0 and 1.
  
  l2heset.par.camera_height = 0.5;
  % This parameter sets the theoretical height above the ground at which
  %   the synthetic hemispheric image will be created. 
  % If l2heset.sw.calc_type < 3, this parameter will already be in the pts
  %   file, allowing for different heights per analysis point. 
  % If l2heset.sw.calc_type = 3, this parameter should be specified here. 
  
  l2heset.par.tolerance     = 0.006;
  % This parameter defines the scaling of the tolerance around each point
  %   in the lidar data for populating the synthetic image matrix. 
  % The number is the "plotting" size of a single LiDAR point at zero
  %   distance to the "camera". The tolerance is then scaled linearly to
  %   half at maximum distance from camera (i.e. if l2heset.par.eval_peri
  %   is set to 100m, at 100m from the camera the tolerance will be half of
  %   that at the camera). Between these two values the plotting size is 
  %   a linear function of the distance to the camera with 10 increments.
  % This parameter can also be used to reduce or increase the tolerance in
  %   the algorithm to suit the point cloud density. 
  % A tolerance = 0.006  is suggested for point cloud densities
  %   between 20-30 points/m-2 (not including added trunk or branche
  %   points). 
  % It is suggested that this value is calibrated using a real 
  %   hemispherical image from a location within your the cloud using the
  %   built in calibration procedure (see documentation). 
  
  l2heset.par.hg_cutoff     = 0;
  % This parameter defines a height of canopy elements above the forest
  %   floor underneath which all points are neglected. A height cutoff
  %   is necessary to a) simulate a hemispherical image taken from a certain
  %   height above the forst floor, and b) to avoid unrealiable output images
  %   due to insufficient point cloud densities as you approach the forest
  %   floor (assuming LiDAR data to be taken from above the canopy). 
  
  l2heset.par.radius        = 500;
  % This parameter relates to the output size of the resulting synthethic
  %   hemispherical images. The value defines the radius of the image 
  %   circle in pixels. Make sure your display has a height of at least 2.5 
  %   times the radius.


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% calculation parameters - shortwave radiation 
  
  l2heset.par.scatter_drad  = [0.533 1.066 2.132]./2;
  % This parameter defines the apparent size of the sun disk in degrees 
  %   of the radius. While the true apparent diameter of the sun is about 
  %   0.533deg, this number can be modified to conceptually account for 
  %   scattered but yet directional swr (see hpeset.sw.scatter above). Up to 
  %   three evaluation disks are supported by HPEval, so that the parameter 
  %   can be a 1x1, 1x2, or 1x3 vector.
  
  l2heset.par.scatter_dwgh  = [60 30 10];
  % This parameter defines the relative weight each of the above sun disk
  %   evaluation will be attributed to. The transmissivity of direct
  %   radiation throught the canopy evaluated based on the first sun disk
  %   scatter_drad(1) will be weighted with scatter_dwgh(1), and so on, to 
  %   arrive at a weighted mean of all evaluation. Note that the the vector
  %   scatter_dwgh must have at least the size of the vector scatter_drad.
  
  l2heset.par.swr_maxdir    = 3;
  % This parameter selects a particular way of constraining direct
  %   radiation from the sun. A limiter is necessary to avoid unrealistic
  %   radiation values in certain situations, such as calculating swr per
  %   inclined surface area when the elevation angle of the sun is very low
  %   and the terrain is facing towards the sun. 
  % Available choices are
  % 0 > using the solar constant
  % 1 > using measured swr transmissivity and modelled split into direct 
  %     and diffuse component
  % 2 > using de beer's law to acount for length of pathway in terms of 
  %     zenith angle. The extinction coefficient is based on the maximum 
  %     transmission of direct swr if sun is at zenith
  % 3 > applying constratins as in 1 and 2, which ever is stricter. This 
  %     is the default option
  
  l2heset.par.time_int      = 10;
  % Time interval for calculating SWR if calculating on a grid. Units are
  %   in minutes i.e. 1 = 1 minute, 0.25 = 15 seconds. 
  
  l2heset.par.stime         = '21.01.2019 00:00:00';
  l2heset.par.etime         = '22.01.2019 00:00:00';
  % Start and end times for calculating sub-canopy incoming shortwave
  %   radiaiton. Should be in format 'dd.mm.yyyy HH:MM:SS'
  % To be used when generating synthetic images at gridded spatial scales
  %   (i.e. if l2heset.sw.calc_type = 3).
  % Leave empty if l2heset.sw.calc_type = 1 or 2.
  
end