clear all

%%% list files 
cd(fileparts(matlab.desktop.editor.getActiveFilename));
currentpath = cd();
 idcs   = strfind(currentpath,'/');
 motherpath = currentpath(1:idcs(end)-1);
 
sourcepath = fullfile(motherpath, 'Data', 'Traces');
segmentspath = fullfile(motherpath, 'Data', 'Events','Events_byRegion_Sel');
dest = fullfile(motherpath,  'Data', 'Events', 'Events_byExp_PlusMSDy');

sourcepath = uigetdir(sourcepath, 'Select Experiment');
cd(sourcepath)

%list tracking result files
dirinfo = dir('*.txt');
dirinfo = dirinfo(~[dirinfo.isdir])
filename = {dirinfo.name}.';
clear dirinfo
section = {'x1y1'; 'x1y2'; 'x2y1'; 'x2y2'; 'x3y1'; 'x3y2'};

indT0 =  find(cell2mat(cellfun (@(x) contains(x, 'start', 'IgnoreCase',true), filename, 'UniformOutput',false)));
indR0 =  find(cell2mat(cellfun (@(x) contains(x, 'result', 'IgnoreCase',true), filename, 'UniformOutput',false)));


%%% Isolate hand-selected ypos section for each trajectory & compute msd %%%

for ii = 1:length(section)
    %find files corresponding to given FOV section
    indS = find(cell2mat(cellfun (@(x) contains(x, section(ii), 'IgnoreCase',true), filename, 'UniformOutput',false))); 
if ~isempty(indS)
    %find file with start/end time for msd evalutaion in given FOV
    indT = intersect(indT0, indS);
    %find file with tracking results for msd evalutaion in given FOV
    indR = intersect(indR0, indS);

    %open file with start/end time  for given FOV
    dataAllT=importdata(filename{indT});
    dataT{ii} = dataAllT.data(:,2:end);
    HeadersT = dataAllT.colheaders(2:end);
    clear dataAllT
    
    %open file with tracking results  for given FOV
    dataAllR=importdata(filename{indR});
    dataR{ii} = dataAllR.data(:,2:end);
    HeadersR = dataAllR.colheaders(2:end);
    clear dataAllR

    %for every segment in given FOV extract position in y between start and end time
    %specified
    for aa = 1:length(dataT{ii})
    yposall = dataR{ii}(find(dataR{ii}(:,9)==dataT{ii}(aa,1)), :);
    trajNo{aa,ii} = dataT{ii}(aa,1);
    yposP{aa,ii} = yposall(find(yposall(:,12)>=dataT{ii}(aa,2) & yposall(:,12)<=dataT{ii}(aa,3)),11);
    clear yposall
        if length(yposP{aa,ii})<=100;
        yposP{aa,ii}=[];
        end
    end
    clear indS indT indR
    end
end

msd0 = cellfun(@(x) (mean((x.*10^-6).^2)-mean(x.*10^-6).^2), yposP,'UniformOutput',false);
msd = cell2mat(msd0);


%%%%% Load segment files and add msd values 

%list segment files
cd(segmentspath)
dirinfo = dir('*.txt');
dirinfo = dirinfo(~[dirinfo.isdir])
filenameS = {dirinfo.name}.';
indSeg = find(contains(filenameS, filename{1,1}(1:10), 'IgnoreCase', true));
filenameSsel = filenameS(indSeg);
clear dirinfo

%load segment files
for ii = 1:length(indSeg)
dataAllS=importdata(filenameS{indSeg(ii)});
dataS{ii} = dataAllS.data;
HeadersS = dataAllS.colheaders;
clear dataAllS
end

%add msd to segments
HeadersS{1,end+1} = 'msd_y';
for ii = 1:length(section);
%find index of segments data corresponding to FOV section
aa = find(contains(filenameSsel, section(ii)));
if ~isempty(aa)
%for every trajectory number in segments,find index in msd matrix
for bb = 1:length(dataS{aa}(:,1))
    if ~isempty (find(cell2mat(trajNo(:,ii))==dataS{aa}(bb,1)))
    indMsd = find(cell2mat(trajNo(:,ii))==dataS{aa}(bb,1));
    msd_sel(bb,1) = msd(indMsd(1),ii);
    else
    msd_sel(bb,1) = NaN;
    end
    clear indMsd
end

ind_rem = find(isnan(msd_sel)); %remove NaNs (traj length <100)
dataS{aa} = horzcat(dataS{aa},msd_sel);
dataS{aa}(ind_rem,:) = [];
clear msd_sel ind_rem
end
end
dataScol= vertcat(dataS{:});

%%%%% Save segments with msd, 1 file per experiment
cd(dest)
ind_name = regexp(filenameS{indSeg(1),1},'_');
exp_name = strcat(filenameS{indSeg(1),1}(1:ind_name(2)-1), filenameS{indSeg(1),1}(ind_name(3):end-4), '_msd','.txt')
fullname = fullfile(dest, exp_name);

fid = fopen(fullname, 'wt');
fprintf(fid,  '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', HeadersS{:});  % header
fclose(fid);
dlmwrite(fullname, dataScol, 'delimiter', '\t','precision', 16, '-append')

    dataCheck=importdata(fullname);
    length(dataCheck.data)
    