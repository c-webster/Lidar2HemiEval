function prepset = L2HEPrep_settings_TEMPLATE

  %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %  % GENERAL DESCRIPTION                                                  %
  %  %   Preparatory script that generates canopy height models and calls
  %  % associated R script that computes canopy maxima and assoicated
  %  % information for input into L2HEval. 
  %  %                                                                      %  
  %  % REQUIREMENTS
  %  %   Functionality is currently written for Windows 10, but should work
  %  % on Windows 7. LASTools is written for Windows using .exe commands,
  %  % which will not run through Linux or OSX systems. Installation of
  %  % 'wine' would enable running of these commands on non-Windows
  %  % systems, but changes must be made below. Contact author for more
  %  % info if required. 
  %  %                                                                      %
  %  % VERSION   
  %  %
  %  % AUTHORS                                                              %
  %  %   C. Webster 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% USER INPUT

prepset.in.lastoolpath   = 'C:/Workspace/LASTools/bin/';
% Full file path to LASTools depository

prepset.in.basefolder    = 'D:/L2HEval_master/';
% Full path to where output data should be kept.
% If running multiple CHMS, the output file structure will be created in
% this folder.
% If only generating one CHM, the output files will be saved into this
% folder.

prepset.in.rpath         = 'C:/"Program Files"/R/R-3.6.1/bin/';
% Full file path to R depository where Rscript.exe is contained
% If file path contains a space (e.g. 'Program Files'), enclose the folder
%   name in quotes i.e. "Program Files"

prepset.in.preppath      = fullfile(prepset.in.basefolder,'Scripts_Associated/L2HE_Prep_R.R');
% Full path and file to location of prep R script. 

prepset.in.lasfile       = 'D:\L2HEval_master\Data_Surface\DSM\Laret_Edge.las';
% prepset.in.lasfile       = fullfile(prepset.in.basefolder,'Data_Surface','DSM','Laret_Edge.las');
% Full path and file to location of large lidar from which to clip smallar
%   point cloud
% If the desired analysis area covers more than one point cloud, suggest
%   using lasmerge.exe before running L2HEPrep.m. The author cautions that
%   the output from lasmerge.exe might not work with the canopy
%   segmentation algorithm due to errors experienced during code
%   development.

prepset.in.cliplas = 0;
% 1 > clips lasfile to within limits and buffer boundaries specified below
% 0 > uses the full extent of the lasfile for CHM calculation and
%     segmentation

prepset.out.tempfolder    = fullfile(prepset.in.basefolder,'TEMP');
% Full path to where output data should be kept. This folder is cleared at
%   the end of the process. 
% If only generating one CHM, the output files will be saved into this
%   folder.

prepset.in.chmstep   = 0.5;
% Grid-size for creating the CHM (in same units as LiDAR data). 
% It is recommended to use 0.5m for the CHM segmentation process in
%   identifying treetops to work. 

prepset.out.saveCHM   = 1;
% Option to keep the CLM within the file structure, or delete after
%   segmentation process

prepset.out.chmfmt    = '.tif';
% Format for saving output chm. Can be '.tif' or '.asc'

prepset.in.xlimits       = [785665,785674]; 
prepset.in.ylimits       = [190995,191004];
% These variables are the minimum and maximum x and y coordinates
%   representing the area you wish to analyse. These values are used to
%   generate the grid points where synthetic images will be calculated. 

prepset.in.buffer        = 100;
% Buffer around analysis area for including in lidar data if clipping a
%   larger .las file. This is
%   important is it will include points outside the core analysis area that
%   are needed to generate the hemispherical images. This value is very
%   closely realted to the eval_peri value in the L2HE_settings file. If you
%   select an evaluation perimeter of 100m, then buffer here should be set to
%   the same value. 
% Units are the same as the lidar point cloud units. 

prepset.in.trunks        = 1;
% setting to calculate trunk points 
% 1 > enable
% 0 > disable
% If enabled, prepset.in.biome must be specified below.

prepset.in.biome         = 19;
% Calculation of diameter at breast height for generating the trunks use
%   the R function itcSegment. A different allometric equation is used
%   depending on the input biome, which must be specified here. This number
%   must be between 0-24. Further information can be found at 
%   https://cran.r-project.org/web/packages/itcSegment/itcSegment.pdf and
%   Jucker et al. 2017: https://doi.org/10.1111/gcb.13388
% Example values are:
% 7  : Nearctic-Boreal forests-Angiosperm
% 8  : Nearctic-Boreal forests-Gymnosperm
% 9  : Nearctic-Temperate coniferous forests-Angiosperm
% 10 : Nearctic-Temperate coniferous forests-Gymnosperm
% 11 : Nearctic-Temperate mixed forests-Angiosperm
% 12 : Nearctic-Temperate mixed forests-Gymnosperm
% 16 : Palearctic-Boreal forests-Angiosperm
% 17 : Palearctic-Boreal forests-Gymnosperm
% 18 : Palearctic-Temperate coniferous forests-Angiosperm
% 19 : Palearctic-Temperate coniferous forests-Gymnosperm
% 20 : Palearctic-Temperate mixed forests-Angiosperm
% 21 : Palearctic-Temperate mixed forests-Gymnosperm

prepset.in.branches      = 1;
% setting to create branch points
% 1 to enable
% 0 to disable
% if enabled, prepset.in.species must be specified below

prepset.in.species          = 'NS';
% Species of the forest is required for selecting branch angle when
%   classifying branch points. Various studies have presented measured and
%   modelled branch angle data for certain species. So far the species
%   supported are: 
% 'NS' : Norway Spruce
% 'DF' : Douglas Fir
% 'SP' : Scot Pine
% 'SB' : Silver Birch
% 'NA' : Uses horizontal branches
% If a forest is mixed, it is recommended the most common species is
%   selected. 

prepset.in.gengrid       = 1;
% Option to generate an analysis grid of evenly spaced points based on x 
%   y limits.
% 1 to enable
% 0 to disable
% Alternative is to specify own points in 'Data_Points' folder. 

prepset.in.spacing       = 1;
% Spacing between points for analysis grid generation. If left
%   empty, the default is 1. Units should be the same as the lidar point
%   cloud. 
% Note: this variable is only used if prepset.in.gengrid = 1, but it should
%   not be left empty.

prepset.in.epsgstr        = '2056';
% When generating the analysis grid, points are removed if they lie within
%   the bounds of the calculated trunks. For faster calculation, trunk
%   coordinates are put into a spatial points dataframe, which requires the
%   EPSG number of the coordinate reference system of the lidar dataset. 
% e.g sample dataset is Swiss grid 1903 (epsg = '2506'). If coordinates are
%   in UTM, the epsg code of the UTM zone must be specified e.g. UTM zone
%   32N (epsg = '32632'). 
% Set to '0000' to avoid requirement for epsg number - this will not remove
%   points that lie within tree trunk bounds. 


end