clear all

%% open a file

cd(fileparts(matlab.desktop.editor.getActiveFilename));
currentpath = cd();
 idcs   = strfind(currentpath,'/');
 motherpath = currentpath(1:idcs(end)-1);
 
sourcepath = fullfile(motherpath,  'Data', 'Traces', '141111_fc2_20ulmin');
dest = fullfile(motherpath, 'Figures');
cd(sourcepath)

%%for the trajectory presented in figure 3a
sourcename = '141111_fc2_x3y2_SDy_results.txt';

%%% for selecting any experiemnt, chose any of the file ending with "_SDyresults":
%%[sourcename,~]=uigetfile('*.txt');


dataAll=importdata(sourcename);
dataNum = dataAll.data(:,1:end);
% for going from SDy short
%dataNum = horzcat([1:1:length(dataNum)].',dataNum);
Headers = dataAll.colheaders;

% create cell array with each cell being single trajectory
changeIndices = dataNum(1:end-1,10)-dataNum(2:end,10);
switchin = find(changeIndices<0);

beads{1} = dataNum(1:switchin(1),:);
for kk = 2:length(switchin)
    beads{kk} = dataNum(switchin(kk-1)+1:switchin(kk),:);
end
beads{kk+1} = dataNum(switchin(kk)+1:end,:);

for kk = 1:length(beads)
    beadno(kk) = beads{kk}(1,10);
end

%% selected trace

%%%%% give trajectory number %%%%%%%%
ii = 672;
namebase = strcat(sourcename(1:16), num2str(ii));

% find index of the trajectory
ind = find(beadno ==ii);
data = beads{ind};


%% plot a trajectory of choice - raw

 set(0,'DefaultAxesFontName', 'Helvetica')
 set(0,'defaultLineLineWidth',0.5); 

 f = figure
  set(gca,'FontSize', 12)

 subplot (3,1,1)  
  plot(data(:,11),data(:,12),'-k','linewidth', 0.5)
    xlabel('x position (um)', 'fontsize', 18)
    ylabel('y position (um)', 'fontsize', 18)
    set(gca, 'xtick', [-50:2:50])
   
  subplot(3,1,2)
   plot(data(1:end,13),data(1:end,11), '-k', 'linewidth', 0.5)
    xlabel('time (s)', 'fontsize', 18)
    ylabel('x position (um)', 'fontsize', 18)

  subplot (3,1,3)  
   plot(data(1:end,13),data(1:end,12),'-k','linewidth', 0.5)
   xlabel('time (s)', 'fontsize', 18)
   ylabel('y position (um)', 'fontsize', 18)
    
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 -2 22 29]);


%% plot a trajectory of choice - center and flip
t_tot = 40;
frame_tot = round(40/0.15);
if frame_tot < length(data)
data1 = data(length(data)-frame_tot:end,:);
else
data1 = data;
end


time1 = data1(:,13)-data1(1,13);
%center data
x1=data1(:,11)-mean(data1(1:end-20,11));
y1=data1(:,12)-mean(data1(1:end-20,12));
%flip to match direction convention in figure
x1=-(data1(:,11)-mean(data1(1:70,11)));
y1=data1(:,12)-mean(data1(1:70,12));
%set plot boundaries for x
xgr = [min(time1) max(time1)];

 set(0,'DefaultAxesFontName', 'Helvetica')
 set(0,'defaultLineLineWidth',0.5); 

 f = figure
  %figure('Visible','off')
  set(gca,'FontSize', 12)
    
  subplot (3,1,1)  
  plot(x1,y1,'-k','linewidth', 0.5)
    xlabel('x (um)', 'fontsize', 18)
    ylabel('y (um)', 'fontsize', 18)
 
    subplot(3,1,2)
   plot(time1,x1, '-k', 'linewidth', 0.5)
    xlabel('time (s)', 'fontsize', 18)
    ylabel('x (um)', 'fontsize', 18)
    xlim(xgr)
    
  subplot (3,1,3)  
   plot(time1,y1,'-k','linewidth', 0.5)
   xlabel('time (s)', 'fontsize', 18)
   ylabel('y (um)', 'fontsize', 18)
    xlim(xgr)
   

cd(dest)
print(fullfile(dest,strcat('Fig3b', '_trajectory_',namebase)), '-dpng', '-r300')


%% Save Source Data
 headersT = {'time_ms', 'x_um', 'y_um'};
dataT = [time1,x1,y1]; 

format = char();
for ii = 1:length(headersT)
    format = strcat(format, '%s\t');
end
format = strcat(format, '\n');

fullname = fullfile(dest, strcat('Fig3b_SourceData_', namebase,'.txt'));

fid = fopen(fullname, 'wt');
fprintf(fid,  format, headersT{:});  % header
fclose(fid);
dlmwrite(fullname, dataT, 'delimiter', '\t','precision', 16, '-append')
