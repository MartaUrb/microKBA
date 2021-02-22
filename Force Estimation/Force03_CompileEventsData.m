clear all

%% list files 
cd(fileparts(matlab.desktop.editor.getActiveFilename));
currentpath = cd();
 idcs   = strfind(currentpath,'/');
 motherpath = currentpath(1:idcs(end)-1);
 
sourcepath = fullfile(motherpath,  'Data', 'Events', 'Events_byExp_PlusForce');
dest = fullfile(motherpath, 'Data', 'Events', 'Events_compiled');

%list tracking result files
cd(sourcepath)
dirinfo = dir('*.txt');
dirinfo = dirinfo(~[dirinfo.isdir]);
filename = {dirinfo.name}.';

for ii = 1:length(dirinfo)
    dataimp = importdata(dirinfo(ii).name);
    dataESel{ii,1} = dataimp.data;
    exp_name{ii,1} = filename{ii}(1:end-7);
    headers = dataimp.textdata;
    headers = replace(headers,'angle (deg)','angle_deg');
    headers = split(headers).';

    clear dataimp
end

dataCol = vertcat(dataESel{:});


%% Save compiled event data
format = char();
for ii = 1:length(headers)
    format = strcat(format, '%s\t');
end
format = strcat(format, '\n');

fullname = fullfile(dest, strcat('KBA_allExperiment_Events_Sel3.txt'));
fid = fopen(fullname, 'wt');
fprintf(fid,  format, headers{:});  % header
fclose(fid);
dlmwrite(fullname, dataCol, 'delimiter', '\t','precision', 16, '-append')
