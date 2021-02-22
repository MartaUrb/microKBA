clear all

%% list files 
cd(fileparts(matlab.desktop.editor.getActiveFilename));
currentpath = cd();
 idcs   = strfind(currentpath,'/');
 motherpath = currentpath(1:idcs(end)-1);
 
sourcepath = fullfile(motherpath,  'Data', 'Events', 'Events_byExp_PlusForce');
dest = fullfile(motherpath, 'Figures');

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


%colors for plotting
cc = [57 177 206;
    123 39 127;
    156 197 69;
    232 74 30]./255;

%% prepare data with classes for selected experiment
ind_ex = [4, 1, 7, 3]; % selected experimnets
for aa = 1:length(ind_ex)
    ii = ind_ex(aa)
    force{aa} = dataESel{ii}(:,14);
    class{aa} = aa.*ones(size(force{aa},1),1);
    medianF (aa,1) = median(force{aa});
end
forceCol = vertcat(force{:});
classCol = vertcat(class{:});
forceT = [classCol, forceCol];

forceSource = nan(max(cellfun(@length, force)),length(force));
v = violinplot(forceT(:,2), forceT(:,1));

for aa = 1:length(v)
boxTable(1,aa) = v(aa).WhiskerPlot.YData(1);
boxTable(2,aa) = v(aa).BoxPlot.Vertices(1,2);
boxTable(3,aa) = v(aa).MedianPlot.YData;
boxTable(4,aa) = v(aa).BoxPlot.Vertices(3,2);
boxTable(5,aa) = v(aa).WhiskerPlot.YData(2);
forceSource(1:length(force{aa}),aa) = force{aa};
end

%% PLOT Force Violin
  set(0,'DefaultAxesFontName', 'Helvetica')
  set(0,'DefaultAxesFontSize', 12)

f = figure
names = {'10','20','30','40'};

f = figure
  f.Units = 'centimeters';
  f.Position = [0 0 7.23 5.75];
boxplot(forceT(:,2), forceT(:,1),'color', 'k','Widths', 0.6, 'Symbol', '');

hold on
for ii = 1:length(unique(forceT(:,1)))
scatter(v(ii).ScatterPlot.XData,v(ii).ScatterPlot.YData,5,cc(ii,:),'o','filled','MarkerFaceAlpha',0.2);
hold on 
end


set(gca, 'xticklabel',names)
set(gca, 'Ytick',[0:2:6])

ylabel ('force (pN)')
xlabel ('flow rate (ul min^{-1})')
box on 
ylim([0 6])

   
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 7.23 5.75]); 

Plotname = strcat('Fig2b_ForceDist');

 print(fullfile(dest, strcat(Plotname,'.png')), '-dpng')
 
 

%% Save Source Data

headersT = {'10ulmin', '20ulmin', '30ulmin', '40ulmin'};

format = char();
for ii = 1:length(headersT)
    format = strcat(format, '%s\t');
end
format = strcat(format, '\n');

fullname = fullfile(dest, strcat('Fig2b_SourceData_Scatter','.txt'));

fid = fopen(fullname, 'wt');
fprintf(fid,  format, headersT{:});  % header
fclose(fid);
dlmwrite(fullname, forceSource, 'delimiter', '\t','precision', 16, '-append')



headersBox = {'10ulmin', '20ulmin', '30ulmin', '40ulmin'};

format = char();
for ii = 1:length(headersBox)
    format = strcat(format, '%s\t');
end
format = strcat(format, '\n');

fullname = fullfile(dest, strcat('Fig2b_SourceData_BoxPlots','.txt'));

fid = fopen(fullname, 'wt');
fprintf(fid,  format, headersBox{:});  % header
fclose(fid);
dlmwrite(fullname, boxTable, 'delimiter', '\t','precision', 16, '-append')

