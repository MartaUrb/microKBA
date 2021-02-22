clear all

% Localization precision: 
% The localization precison was calculated using background traces from exp 141112_fc2
% To determine the localization precision for bead tracking, 
% we estimated the standard deviation of x and y positions for 12 beads stuck to 
% the surface over the time of 5 min. The used traces were corrected for drift 
% by subtracting the average of all stuck bead traces. 

%% data opening and structuring

% create list of all txt files in folder
cd(fileparts(matlab.desktop.editor.getActiveFilename))
fileList = dir('*.txt')
fileList = fileList(~[fileList.isdir]) %clean some subdirectories
fileList(contains({fileList.name}, '._')) = [];

namebase = fileList(1).name(1:10);
% create cell array with data for one background bead in each cell 

for ii = 1:numel(fileList)
    dataall{ii} = dlmread(fileList(ii).name, '\t', 1, 0);
end
 
% make the trajectory field contain no for each bead
for ii = 1:numel(fileList)
    ind_b1 = regexp(fileList(ii).name,'bck\d*', 'start');
    ind_b2 = regexp(fileList(ii).name,'bck\d*', 'end');
    dataall{ii}(:,9) = str2num(fileList(ii).name(ind_b1+3:ind_b2));
    bead(ii,1) = str2num(fileList(ii).name(ind_b1+3:ind_b2));
end
[~,ind] = sort(bead);
bead = bead(ind);
dataall = dataall(ind);
clear ind

imp = importdata(fileList(1).name);
headers = imp.colheaders;
clear imp

pixelsize = 1.6;
frames = max(cellfun(@length, dataall));
frames_perbead = cellfun(@length, dataall);

data=dataall;

%%%%%%%%%%%%%%%%%%% Discard some beads %%%%%%%%%%%%%%%%%%
ind_incomplete = find(frames_perbead<frames); %beads not lasting whole stream
for ii =1:length(dataall)
    stdxall(ii,1) = std(dataall{ii}(:,4));
    stdyall(ii,1) = std(dataall{ii}(:,5));
end
ind_badquality = find(stdxall>=0.099); %with high stdy
ind_remove = union(ind_badquality, ind_incomplete);

data(ind_remove) = [];

%% Calculate Background trace
clear x_corr y_corr
x_corr = NaN(frames,length(data));
y_corr = NaN(frames,length(data));

for ii = 1:length(data)
    avx(ii,1) = mean(data{ii}(:,4));
    avy(ii,1) = mean(data{ii}(:,5));
    x_corr(1:frames_perbead(ii),ii) = (data{ii}(1:frames_perbead(ii),4)-avx(ii));
    y_corr(1:frames_perbead(ii),ii) = (data{ii}(1:frames_perbead(ii),5)-avy(ii));
end

clear data_bck
data_bck(:,1) = nanmean(x_corr,2);
data_bck(:,2) = nanmean(y_corr,2);

%% subtract background and calculate sd
frame_st = round(200./0.15);
frame_end = frame_st+300/0.15;
for aa = 1:length(data)
data{aa}=data{aa}(frame_st:frame_end,:);
end

for ii = 1:numel(data);
    data{ii}(:,10) = data{ii}(:,4)-data_bck(frame_st:frame_end,1);
    data{ii}(:,11) = data{ii}(:,5)-data_bck(frame_st:frame_end,2);
end

for ii = 1:length(data)
    stdx(ii,1) = nanstd(data{ii}(:,10));
    stdy(ii,1) = nanstd(data{ii}(:,11));
end

avstdx = mean(stdx).*pixelsize
avstdy = mean(stdy).*pixelsize
