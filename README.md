
Note: A more efficient and flexible version of this code is available as L2R/las2rad in the [CanRad](https://github.com/c-webster/CanRad.jl) julia package


## Lidar2HemiEval

Three companion scripts that together:
	1. Enhances airborne LiDAR data of forests to include opague trunks and increase tree crown density
	2. Creates synthetic hemispherical images from enhanced LiDAR data
	3. Calculates sky-view fraction and time-varying canopy transmissivity for estimating total incoming sub-canopy shortwave radiation


 ### CITATION (to be updated)
Webster, C., Mazzotti, G., Essery, R., and Jonas T. (in review) Enhancing airborne LiDAR data for improved forest structure representation in shortwave transmission models. Remote Sensing of Environment. 


 ### REQUIREMENTS
Scripts were written and tested for Matlab 2019a and R version 3.6.1 (2019-07-05) 
   - LAStools by Martin Isenburg (https://rapidlasso.com/lastools/) - a license is required for larger areas/denser Lidar point clouds. 
   - lasdata.m package from Matlab File Exchange by Teemu Kumpumäki(https://mathworks.com/matlabcentral/fileexchange/48073-lasdata) (first accessed 2018/10)
   - Usage of the Matlab 'system' command for running LAStools is currently implemented for Windows only. If running on OSX/Linux, installation of wine and relevant adjustments to Lidar2HemiEval_Prep.m are both required.
  
   
 ### CONTRIBUTING MATERIAL
   - Calculation of canopy height models follows 'generating_a_pit_free_chm.bat' from LAStools, based on http://www.riegl.com/uploads/tx_pxpriegldownloads/khosravipour_SilviLaser2013.pdf
   - R packages rLiDAR, raster, rlas, itcSegment, sp and pracma. These are installed/added whenever L2HE_Prep_R.R is run. 
   - Solar radiation algorithm from HPEval (https://github.com/Tobias-Jonas-SLF/HPEval)
   - Solar position calculations from NOAA, found online at www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html (accessed 2018/06)
   - Function to split solar radiation into direct/diffuse component as in the factorial snow model FSM 1.0 by R. Essery available at www.geosci-model-dev.net/8/3867/2015 (accessed 2018/06)
   - Function utm2lonlat.m from Matlab File Exchange to convert UTM into WGS84 coordinates (https://mathworks.com/matlabcentral/fileexchange/44242-utm2lonlat) (accessed 2018/08)


 ### EXECUTION 
From the Matlab command line:
Lidar2HemiEval_Prep('L2HEPrep_settings.m') [Note L2HE_Prep.R is run from inside this script]
Lidar2HemiEval('L2HEval_Settings.m')

All settings such as paths to I/O data, input parameters and switch settings are handled through the external settings files, which is the only input argument to the scripts. 
Settings for each script are explained in detail in the *_settings.m files. 

Specific comments:
- lasdata.m must be downloaded from the file exchange and placed in the Matlab filepath (see contributing material)


 ### SAMPLE DATASET
   - Lidar data over the Laret field area in Davos, Switzerland
   - Local digital terrain model at 1 m resolution ('DTM') and regional digital terrain model at 50 m resolution ('DEM')
   - Evaluation points in three formats (see Lidar2HemiEval for further details)
   - Shortwave radiation data measured at open site 200m from centre of Lidar area. 
   All spatial datasets are in the SwissGrid coordinate system CH1903/LV03 EPSG:2056. 


 ### DESCRIPTION
	1. Lidar2HemiEval_Prep.m
		i. Calculates a canopy height model over the field area
		ii. Gathers input and runs L2HE_Prep_R.R
	
	1a. L2HE_Prep_R.R
		i. Identifies canopy maxima and segments the canopy into individual tree crowns following http://quantitativeecology.org/using-rlidar-and-fusion-to-delineate-individual-trees-through-canopy-height-model-segmentation/.
		ii. Calculates diameter at breast height (dbh) using equations from Jucker et al. 2017 (doi: 10.1111/gcb.13388) for each tree. Results are saved with x/y coordinates of each tree maxima and tree height. 
		iii. Each point in the lidar data is classified to the tree crown to which is lies within. Points outside of segmented tree crowns are excluded from the output. Each tree crown is random assigned a branch angle.
		iv. Analysis points evenly spaced along a grid (as per user settings) are generated, and points that lie within trunks are excluded. Data is saved as 

	2. Lidar2HemiEval.m
		i. Trunk points are calculated using data in .*_treeinfo.txt'. Cylinders are created up to dbh, then the trunk linearly tapers up to tree height. Data is saved as 'tsm' structure. 
		ii. Branch points are calculated using data in .*_lastreeclass.txt'. Additional points at 0.15 m spacing are calculated along a Euclidean trajectory between the classified lidar point and the tree trunk at the crown-specific branch angle. Data is saved as 'bsm' structure.
		iii. The script loops over each input coordinate (e.g. in '*_gridpts.txt') and calculates:
			a. Synthetic hemispherical images at each point
			b. Sky-view fraction for both hemispherical and flat sensor perspective
			c. Time-varying canopy direct shortwave transmissivity
			Progress of loop, and timing of each step is recorded in the output folder.


 ### OUTPUT
The scripts have several optional outputs (outlined in more detail in the settings files):

Lidar2HemiEval_Prep
	1. Canopy height model over field area as .tif or .asc
	2. Locations and allometric statistics of trunks as input for Lidar2HemiEval ('Data_Surface/DBH/*_treeinfo.txt')
	3. Lidar points classified to each tree crown as input for Lidar2HemiEval ('Data_Surface/LTC/*_lastreeclass.txt')
	4. Grid point locations for calculating synthetic hemispherical images ('Data_Points/*_gridpts.txt')

Lidar2HemiEval
	1. Synthetic hemispherical images at each location as 'SHI_*.mat' files. See 'MakeSHIs.m' to convert to .png files. 
	2. Sky-view fraction for either a flat or hemispherical sensor as 'SVF_calculated.txt'  as '*.mat'
	3. Time-varying direct shortwave transmissivity, incoming total, diffuse and direct shortwave radiation in Wm/2 for the given time period as 'SWR_*.mat'

