function MakeSHIs


datdir   = 'D:\L2HEval_master\Output_Data\Laret_Edge\SHIs';


    files = dir(fullfile(datdir,'*SHI*.mat'));

    if isempty(files)
        disp('No SHI.mat files found in given directory')
    end


    for fx = 1:length(files)

        load(fullfile(datdir,files(fx).name));

        imwrite(image2ev,fullfile(datdir,strcat(files(fx).name(1:end-4),'.png')))

    end
     
end
