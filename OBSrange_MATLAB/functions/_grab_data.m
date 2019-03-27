%% get the OBS survey data files

% directories
wd = '/Users/zeilon/Documents/MATLAB/OBS_range_dist/LOCATE_OBSv3/data/';
datadir = '~/Work/Fieldwork/PacArray18/OBS_survey/OrcaAcousticSurvey/';
odir    = 'OrcaAcousticSurvey';
odirsio = 'SIO_relocated';

%% do SIO correction
cd(datadir);
!./SIOlocate_all.bash
!for i in *; do mv $i `echo $i sed s/_corr/_SIOcorr/g`; done'
cd(wd)

%% get station files and copy across
files = dir(datadir);
for ii = 1:length(files)
    if ~strcmp(files(ii).name(1),'.') && ~any(regexp(files(ii).name,'_corr')) && ~any(regexp(files(ii).name,'bash'))
        copyfile([datadir,files(ii).name],odir)
    end
end

%% get SIO correction files
files = dir(datadir);
for ii = 1:length(files)
    if any(regexp(files(ii).name(1),'_SIOcorr'))
        copyfile([datadir,files(ii).name],odirsio)
    end
end

