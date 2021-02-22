
clear all

%% load data
cd(fileparts(matlab.desktop.editor.getActiveFilename));
currentpath = cd();
 idcs   = strfind(currentpath,'/');
 motherpath = currentpath(1:idcs(end)-1);
 
sourcepath = fullfile(motherpath,  'Data', 'Events', 'Events_compiled');
dest = fullfile(motherpath, 'Figures');
cd(sourcepath)

%%compiled event data
sourcename = 'KBA_allExperiment_Events_Sel3.txt';

dataimp = importdata(fullfile(sourcepath,sourcename));
data = dataimp.data;
headers = dataimp.textdata;
headers = split(headers).';
clear dataimp
indF = 15;

anglesel = 45;
ind_ang = [find(abs(data(:,11))<=anglesel); find(abs(data(:,11))>=(180-anglesel))];
data = data(ind_ang,:);

Fsel = [[-3.3:0.3:-0.6],[0.6:0.3:3.3]];
indW = (find(Fsel<0));
indA = (find(Fsel>0));

% prepare data binned for Force 
 for   ii = 1:length(Fsel);
   x = Fsel(ii);
   pm = 0.15;
   xL = x - pm;
   xH = x +pm;
   dataBin{ii} = data(data(:,indF)>=xL&data(:,indF)<=xH,:);
 end
 
set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 10)

%plot colors
    colorA = [239, 76, 67]./255;
    colorW = [85, 61, 151]./255;
 
%% Estimate run lengths and interaction times, bin by bin
dataFit(:,1) = Fsel;
headersFits = {'F', 'mu_RL', 'dmu_RL', 'x0_RL', 'mu_IT', 'dmu_IT', 'x0_IT'};

for ii = 1:length(dataBin)
    
runLengths = dataBin{1,ii}(:,9);
itimes = dataBin{1,ii}(:,8);

%estimate run lengths
obs = runLengths; %vector of measurements
s = zeros(1,100);
x0=s;
for m=1:100
   pk = randi(length(obs),1,length(obs));
   resample = obs(pk);
   [s(m),x0(m)] =mexpfit(resample,[]);
end
dataFit(ii,2) = mean(s);
dataFit(ii,3) = 2*std(s);
dataFit(ii,4) = mean(x0);

clear obs

%estimate interaction times
obs = itimes; %vector of measurements
s = zeros(1,100);
x0=s;
for m=1:100
   pk = randi(length(obs),1,length(obs));
   resample = obs(pk);
   [s(m),x0(m)] =mexpfit(resample,[]);
end
dataFit(ii,5) = mean(s);
dataFit(ii,6) = 2*std(s);
dataFit(ii,7) = mean(x0);

clear obs
clear runLengths itimes
end


%% plot RL vs Force plus patch

% rL_meanW = mean(dataFit(1:9,2));
% rL_stdW = std(dataFit(1:9,2));
% 
% rL_meanA = mean(dataFit(10:15,2));
% rL_stdA = std(dataFit(10:15,2));

rL_meanAll = mean(dataFit(2:18,2))
rL_stdAll = 2.*std(dataFit(2:18,2))

  set(0,'DefaultAxesFontName', 'Arial')
  set(0,'DefaultAxesFontSize', 10)

    f = figure
xline(0)
hold on
patch([-4 -4 4 4 ],[rL_meanAll-rL_stdAll, rL_meanAll+rL_stdAll, rL_meanAll+rL_stdAll, rL_meanAll-rL_stdAll],'k', 'facealpha', 0.1, 'edgecolor', 'none');
hold on
plot([-4,4], [rL_meanAll,rL_meanAll], '--','color','k', 'linewidth', 1.5)
hold on
    errorbar(dataFit(indW,1), dataFit(indW,2),dataFit(indW,3), '.', 'linewidth',1, 'color', colorW)
    hold on
    errorbar(dataFit(indA,1), dataFit(indA,2),dataFit(indA,3), '.', 'linewidth',1, 'color', colorA)
    hold on
    scatter(dataFit(indW,1), dataFit(indW,2), 15.8,'o', 'fill', 'markerfacecolor', colorW, 'markeredgecolor', 'none')
    hold on
    scatter(dataFit(indA,1), dataFit(indA,2), 15.8,'o', 'fill', 'markerfacecolor', colorA, 'markeredgecolor', 'none')

   xlim([-4,4])
   ylim([0, 1.0])
   ax = gca;
    ax.FontSize = 8;
    ax.XTick = [-6:1:6]
   xlabel ('force (pN)', 'FontSize',10)
   ylabel('run length (um)', 'FOntSize',10)
   box on 
   grid off
      
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 7.67 3.2]); 

Plotname = strcat('Fig4b_RunLengths');

print(fullfile(dest,strcat(Plotname,'.png')), '-dpng')


%% plot itime vs force plus patch

iT_meanW = mean(dataFit(1:10,5));
iT_stdW = 2.*std(dataFit(1:10,5));

iT_meanA = mean(dataFit(11:19,5));
iT_stdA = 2.*std(dataFit(11:19,5));

  set(0,'DefaultAxesFontName', 'Arial')
  set(0,'DefaultAxesFontSize', 10)
    
    
 f = figure
xline(0)
hold on
patch([-4 -4 0 0],[iT_meanW-iT_stdW, iT_meanW+iT_stdW, iT_meanW+iT_stdW, iT_meanW-iT_stdW],'k', 'facealpha', 0.1, 'edgecolor', 'none');
hold on
plot([-4,0], [iT_meanW,iT_meanW],  '--','color','k', 'linewidth', 1.5)
hold on
patch([0 0 4 4],[iT_meanA-iT_stdA, iT_meanA+iT_stdA, iT_meanA+iT_stdA, iT_meanA-iT_stdA],'k', 'facealpha', 0.1, 'edgecolor', 'none');
hold on
plot([0,4], [iT_meanA,iT_meanA],  '--','color','k', 'linewidth', 1.5)
hold on

    errorbar(dataFit(indW,1), dataFit(indW,5),dataFit(indW,6), '.', 'linewidth',1, 'color', colorW)
    hold on
    errorbar(dataFit(indA,1), dataFit(indA,5),dataFit(indA,6), '.', 'linewidth',1, 'color', colorA)
    hold on
    scatter(dataFit(indW,1), dataFit(indW,5), 20,'o', 'fill', 'markerfacecolor', colorW, 'markeredgecolor', 'none')
    hold on
    scatter(dataFit(indA,1), dataFit(indA,5), 20,'o', 'fill', 'markerfacecolor', colorA, 'markeredgecolor', 'none')

   xlim([-4,4])
   ylim([0, 2.2])
   ax = gca;
    ax.FontSize = 8;
    ax.XTick = [-6:1:6]
        ax.YTick = [0:0.5:2]

   xlabel ('force (pN)', 'FontSize',10)
   ylabel('interaction time (s)', 'FOntSize',10)
   box on 
   grid off
   
   cd(dest)
   
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 7.67 3.2]); 

Plotname = strcat('Fig4c_IntTimes');
print(fullfile(dest,strcat(Plotname,'.png')), '-dpng')




%% Save Bin Data
format = char();

for ii = 1:length(headersFits)
    format = strcat(format, '%s\t');
end
format = strcat(format, '\n');

fullname = fullfile(dest, strcat('Fig4bc_RunL_IntTime_BinData','.txt'));
fid = fopen(fullname, 'wt');
fprintf(fid,  format, headersFits{:});  % header
fclose(fid);
dlmwrite(fullname, dataFit, 'delimiter', '\t','precision', 16, '-append')

