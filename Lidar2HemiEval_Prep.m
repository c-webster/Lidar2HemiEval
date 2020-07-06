function Lidar2HemiEval_Prep(sfile)

%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  % GENERAL DESCRIPTION                                                  %
%  %   Preparatory script that generates canopy height models 
%  %                                                                      %
%  % VERSION                                                              %
%  %   v1.0 updated 06.07.2020                                            %
%  %   
%  % AUTHOR:                                                              %
%  %   Clare Webster  (1,2)
%  %    (1) WSL Institute for Snow and Avalanche Research SLF, Davos, CH  %
%  %    (2) University of Edinburgh, School of GeoSciences, Edinburgh, UK %

%  % USAGE
%  %   > Lidar2HemiEval_Prep('L2HEPrep_settings.m')

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %% LOAD SETTINGS
  if nargin == 0
    error('No settings file specified')
  end

  [setp,setf,~]    = fileparts(sfile);
  cd(setp);

  try
    if exist(sfile,'file')
        prepset  = feval(setf);
    else
        error('Specified settings folder is not in directory')
    end
  catch
        error('Error while reading L2HEPrep settings file')
  end

  % Check input compatability
  if prepset.in.gengrid ~= 0 && prepset.in.gengrid ~= 1
    error('Unknown setting for grid generation. Options are 0 or 1 to disable or enable functionality')
  end

  if ~isempty(prepset.out.chmfmt) && ~strcmp(prepset.out.chmfmt,'.tif') && ~strcmp(prepset.out.chmfmt,'.asc')
    error('Unknown CHM output format. Options are .tif or .asc or empty')
  end

  if (strcmp(prepset.out.chmfmt,'.asc') || isempty(prepset.out.chmfmt)) &&...
            (prepset.in.trunks || prepset.in.branches)
    disp('Warning: CHM output format changed from to .tif to enable calculation of trunks and branches')
    prepset.in.chmfmt = '.tif';
  end

  if prepset.in.biome < 0 || prepset.in.biome > 24
    error('Unknown input biome for calculating diameter at breast height')
  end

  % Species settings
  switch prepset.in.species
    case 'NS'
        species = 1;
    case 'DF'
        species = 2;
    case 'SP'
        species = 3;
    case 'SB'
        species = 4;
    case 'NA'
        species = 0;
    otherwise
        error('Input species not recognised');
  end

  % Check spacing input
  if isempty(prepset.in.spacing) && prepset.in.gengrid == 1
    disp('Grid spacing left empty. Using default of 1 metre')
    prepset.in.spacing = 1;
  elseif isempty(prepset.in.spacing) && prepset.in.gengrid == 0
    prepset.in.spacing = 1; % will be unused but variable needs value for R script to run
  end


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Clip large lidar to analysis area + buffer

  tmpdir = fullfile(prepset.in.basefolder,'temp');
  if ~exist(tmpdir,'dir')
    mkdir(tmpdir)
  end
  
  tx1 = prepset.in.xlimits(1); 
  tx2 = prepset.in.xlimits(2);
  ty1 = prepset.in.ylimits(1);
  ty2 = prepset.in.ylimits(2);  
  
   
  [~,lasname] = fileparts(prepset.in.lasfile);

  if prepset.in.cliplas 
      new(1,:) = [tx1-prepset.in.buffer,ty1-prepset.in.buffer];
      new(2,:) = [tx2+prepset.in.buffer,ty1-prepset.in.buffer];
      new(3,:) = [tx2+prepset.in.buffer,ty2+prepset.in.buffer];
      new(4,:) = [tx1-prepset.in.buffer,ty2+prepset.in.buffer];
      new(5,:) = [tx1-prepset.in.buffer,ty1-prepset.in.buffer];

      clipfilename = fullfile(tmpdir,strcat(lasname,'_clipdim.txt'));
      dlmwrite(clipfilename,new,'Precision',12,'Delimiter','\t');

      fprintf('\n Clipping Lidar.... \n \n')

      lasname = strcat(lasname,'_clipped');
      lasfilename = fullfile(prepset.in.basefolder,'Data_Surface','DSM',strcat(lasname,'.las'));
      system([prepset.in.lastoolpath 'lasclip.exe -i ' prepset.in.lasfile ' -poly '...
          clipfilename ' -o ' lasfilename ]);
        
  else
    lasfilename = prepset.in.lasfile;
  end
  
   if prepset.in.gengrid
    outdir = fullfile(prepset.in.basefolder,'Output_Prep',lasname);
    if ~exist(outdir,'dir'); mkdir(outdir); end
      dlmwrite(fullfile(outdir,strcat(lasname,'_analysisarea.txt')),[tx1,tx2,ty1,ty2],'Precision',12,'Delimiter','\t')
   end
  
  [chmfname, normlidar] = create_CHM(lasfilename,prepset,tmpdir);
    
  if prepset.in.trunks || prepset.in.branches
    
    if prepset.in.trunks; tdir = fullfile(prepset.in.basefolder,'Data_Surface','DBH'); if ~exist(tdir,'dir'); mkdir(tdir); end; end
    if prepset.in.branches; bdir = fullfile(prepset.in.basefolder,'Data_Surface','LTC'); if ~exist(bdir,'dir'); mkdir(bdir); end; end

      fprintf('\n Starting segmentation in R.... \n \n')
          try
            system([prepset.in.rpath 'Rscript.exe ' prepset.in.preppath ' ' prepset.in.basefolder ' '...
              chmfname ' ' normlidar ' ' lasname ' ' num2str(prepset.in.trunks) ' ' num2str(prepset.in.branches) ' ' ...
              num2str(prepset.in.gengrid) ' ' num2str(prepset.in.spacing) ' ' ...
              prepset.in.epsgstr ' ' num2str(prepset.in.biome) ' ' num2str(species)]); 
          catch
            error('Failed to execute R segmentation script. Check R file path / version number')
          end

      r_finfile = fullfile(tmpdir,'R_finish_confirmation.txt');
      if exist(r_finfile,'file')
        fprintf('\n Segmentation in R finished successfully.... \n')
        delete(r_finfile) % then delete it because it's not needed
      else 
        error('Failure in R segmentation script')
      end

      delete(normlidar) % deleted here and not in make CHM function because the segmentation algorithm needs it
      
  end

  
  % clear tmpdir   
   if exist(tmpdir,'dir')
     rmdir(tmpdir,'s')
   end
   
   if ~prepset.out.saveCHM
     delete(chmfname)
     if ~prepset.in.gengrid
       delete(fullfile(prepset.in.basefolder,'Output_Prep',lasname))
     end
   end
   
  fprintf('\n Finished. \n \n')



end

function [chmfname, normlidar] = create_CHM(lasfilename,prepset,tmpdir)

    if ~exist(tmpdir,'dir'); mkdir(tmpdir); end

    [~,lasname] = fileparts(lasfilename);
    outdir = fullfile(prepset.in.basefolder,'Output_Prep',lasname);
    if ~exist(outdir,'dir'); mkdir(outdir); end
    chmfname    = fullfile(prepset.in.basefolder,'Output_Prep',lasname,strcat(lasname,'_chm',prepset.out.chmfmt));
    
    grndlidar = fullfile(tmpdir,strcat(lasname,'_ground.laz'));
    normlidar = fullfile(tmpdir,strcat(lasname,'_norm.laz'));  

    fprintf('\n Calculating ground points... \n')
       
    system([prepset.in.lastoolpath 'lasground_new.exe -i ' lasfilename ' -wilderness -o ' grndlidar]);

    fprintf('\n Making CHM.... \n \n')
    
    %%% The following steps were translated into Matlab from
    %%% 'generating_a_pit_free_chm.bat' from LAStools:
    %%% :: a batch script for generating a pit-free CHM as outlined
    %%% :: in the Silvilaser 2013 poster by A. Khosravipour et al.

    % CHM settings
    step = prepset.in.chmstep;
    kill = 1;

    % create normalised point cloud
    system([prepset.in.lastoolpath 'lasheight.exe -i ' grndlidar ' -replace_z -o ' normlidar]);

    % make the chm   
    system([prepset.in.lastoolpath 'blast2dem.exe -i ' normlidar ' -keep_first -step ' num2str(step) ...
        ' -odir ' tmpdir ' -odix _00 -obil']);

    system([prepset.in.lastoolpath 'blast2dem.exe -i ' normlidar ' -keep_first -drop_z_below 2 -step ' ...
        num2str(step) ' -kill ' num2str(kill) ' -odir ' tmpdir ' -odix _02 -obil']);

    system([prepset.in.lastoolpath 'blast2dem.exe -i ' normlidar ' -keep_first -drop_z_below 5 -step ' ...
        num2str(step) ' -kill ' num2str(kill) ' -odir ' tmpdir ' -odix _05 -obil']);     

    system([prepset.in.lastoolpath 'blast2dem.exe -i ' normlidar ' -keep_first -drop_z_below 10 -step ' ...
        num2str(step) ' -kill ' num2str(kill) ' -odir ' tmpdir ' -odix _10 -obil']); 

    system([prepset.in.lastoolpath 'blast2dem.exe -i ' normlidar ' -keep_first -drop_z_below 15 -step ' ...
        num2str(step) ' -kill ' num2str(kill) ' -odir ' tmpdir ' -odix _15 -obil']);

    system([prepset.in.lastoolpath 'blast2dem.exe -i ' normlidar ' -keep_first -drop_z_below 20 -step ' ...
        num2str(step) ' -kill ' num2str(kill) ' -odir ' tmpdir ' -odix _20 -obil']);       

    system([prepset.in.lastoolpath 'blast2dem.exe -i ' normlidar ' -keep_first -drop_z_below 25 -step ' ...
        num2str(step) ' -kill ' num2str(kill) ' -odir ' tmpdir ' -odix _25 -obil']); 

    system([prepset.in.lastoolpath 'blast2dem.exe -i ' normlidar ' -keep_first -drop_z_below 30 -step ' ...
        num2str(step) ' -kill ' num2str(kill) ' -odir ' tmpdir ' -odix _30 -obil']);           

    system([prepset.in.lastoolpath 'lasgrid -i ' fullfile(tmpdir,'*.bil') ' -merged -step ' ...
        num2str(step) ' -highest -o ' chmfname]);

    delete(grndlidar)

    fprintf('\n ....done \n \n')


   
end