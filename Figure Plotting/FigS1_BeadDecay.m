clear all

cd(fileparts(matlab.desktop.editor.getActiveFilename));
currentpath = cd();
 idcs   = strfind(currentpath,'/');
 motherpath = currentpath(1:idcs(end)-1);
 clear idcs
 
sourcepath = fullfile(motherpath, 'Data', 'PeakCount');

cd(sourcepath)

dirinfo = dir('*.txt');
dirinfo = dirinfo(~[dirinfo.isdir])
filename = {dirinfo.name}.';
experiemnt = cellfun(@(x) x(1:21), filename, 'UniformOutput', false);
experiemnt = cellfun(@(x) strrep(x,'_',' '), experiemnt,'UniformOutput', false);

for ii = 1:length(dirinfo)
peakcounts{ii} = dlmread(filename{ii}, '\t', 1, 1);
peakcounts{ii}(:,3) = peakcounts{ii} (:,1)*0.15;
peakcounts{ii}(:,4) = peakcounts{ii} (:,2)./max(peakcounts{ii}(:,2));
end

cmap = parula(length(peakcounts));
%% plotting peaks all experiemnts

h = figure
set(0,'DefaultAxesFontName', 'Myriad Pro')
set(gca,'FontSize', 16)

for ii = 1:length(peakcounts)
    plot (peakcounts{ii}(:,3),peakcounts{ii}(:,4), 'color', cmap(ii,:),'linewidth', 4)
    xlabel('time [s]', 'fontsize', 20)
    ylabel('Relative Counts', 'fontsize', 20)
     title('Decrease of Beads Number over Time', 'fontsize', 24)
hold on
end

grid on 


 ylim ([0, 1.19])   
 xlim ([0, 800])

desc =  legend(experiemnt)

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 20 15]);
  
%% plotting cnt vs selected

h = figure
set(0,'DefaultAxesFontName', 'Myriad Pro')
set(gca,'FontSize', 16)

    plot (peakcounts{1}(:,3),peakcounts{1}(:,4), 'color', [0.2 0.2 0.2],'linewidth', 4)
    xlabel('time [s]', 'fontsize', 20)
    ylabel('Relative Counts', 'fontsize', 20)
hold on
    plot (peakcounts{2}(:,3),peakcounts{2}(:,4), 'color', [0.9, 0.1, 0.1],'linewidth', 4)

     title('Decrease of Beads Number over Time', 'fontsize', 24)

grid on 


 ylim ([0, 1.19])   
 xlim ([0, 400])

desc =  legend(experiemnt)

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 20 15]);
  
