clear all

%% import data

cd(fileparts(matlab.desktop.editor.getActiveFilename));
  currentpath = cd();
  idcs   = strfind(currentpath,'/');
  motherpath = currentpath(1:idcs(end)-1);
 
% import MTs measurments data
sourcepathMT = fullfile(motherpath,  'Figures', 'FigS3_MTs', '141111_fc2_MT_distribution');
cd(sourcepathMT)
dirinfo = dir('*.xls')
dataAll = [];
for ii = 1:length(dirinfo)
    datatemp =dlmread(dirinfo(ii).name,'\t',1,0);
    dataAll = vertcat(dataAll, datatemp);
    clear datatemp
end
headers = {'no','BX','BY','Width','Height','Angle', 'Length'};
    
angleMT = dataAll(:,6);

indN = find(angleMT<-90);
angleMT(indN) = 180-abs(angleMT(indN));
indB = find(angleMT>90);
angleMT(indB) = abs(angleMT(indB))-180;

bars = 39;

% import events data from the corresponding experimnet
sourcepathEv = fullfile(motherpath,  'Data', 'Events', 'Events_byExp_PlusForce');
%sourcepathEv = fullfile(motherpath,  'Data', 'Events', 'Events_compiled');
dest = fullfile(motherpath, 'Figures');
%sourcenameEv = 'KBA_allExperiment_Events_Sel3.txt';
sourcenameEv ='141111_fc2_segments_sel3_t120_Force.txt';

dataimp = importdata(fullfile(sourcepathEv,sourcenameEv));
dataEv = dataimp.data;
headers = dataimp.textdata;
headers = split(headers).';
clear dataimp
dataSelAgainst = dataEv(abs(dataEv(:,11))<90,:);
dataSelWith = dataEv(abs(dataEv(:,11))>90,:);

%%%define colors for plots
    colorA = [239, 76, 67]./255;
    colorW = [85, 61, 151]./255;
    colorMT = [13 178 174]./255;

%% generate histogram data

%%%%% building length-weighted histogram for microtubules
set(0,'DefaultAxesFontName', 'Helvetica')

barN = bars;
vv = angleMT;
ww = dataAll(:, 7)
ee = linspace (-90,90,(barN+1));
      centers = linspace (-90,90, barN)
    
    if ~isvector(vv) || ~isvector(ww) || length(vv)~=length(ww)
        error('vals and weights must be vectors of the same size');
    end
    
    Nedge = length(ee);
    h = zeros(size(ee)-1);
    
  for n = 1:Nedge-2
      ind = find(vv >= ee(n) & vv < ee(n+1));
      h(n) = sum(ww(ind))/sum(ww);
  end
  
  for n = Nedge-1
      ind = find(vv >= ee(n) & vv <= ee(n+1));
      h(n) = sum(ww(ind))/sum(ww);
  end


%%% build histogram for events data

% building histogram from defined edges AGAINST
clear centersEA, clear hsEA, clear edgesEA, clear RangeA
barNA = bars;
valuesEA = dataSelAgainst(:,11);
edgesEA = linspace (-90, 90, (barNA+1));
edgesEA(end)=edgesEA(end)+1;  

    NEA = length(edgesEA);
    hsEA = zeros(size(edgesEA)-1);
    RangeA = edgesEA(2)-edgesEA(1);
    
    for n = 1:NEA-1
        indA{n} = find(valuesEA >= edgesEA(n) & valuesEA < edgesEA(n+1));
        hsEA(n) = length(indA{n})/length(dataSelAgainst);
        centersEA (n) = edgesEA(n)+RangeA./2;
    end
    
  % building histogram from defined edges WITH
clear centersEW, clear hsEW 

barNW = bars;
anglesW = dataSelWith(:,11);

indN = find(anglesW<-90);
anglesW(indN) = 180-abs(anglesW(indN));
indB = find(anglesW>90);
anglesW(indB) = abs(anglesW(indB))-180;
valuesEW=anglesW;
edgesEW = linspace (-90, 90, (barNW+1));
    
    NEW = length(edgesEW);
    hsEW = zeros(size(edgesEW)-1);
    RangeW = edgesEW(2)-edgesEW(1);
    
    for n = 1:NEW-1
        indW = find(valuesEW >= edgesEW(n) & valuesEW < edgesEW(n+1));
        hsEW(n) = length(indW)/length(dataSelWith);
        centersEW (n) = edgesEW(n)+RangeW./2;
    end

%% Plot all three histograms
    
set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 8)

   f = figure('Position',    [440   403   323   395])
   subplot(3,1,1)
   stairs(ee,[h,0], 'b-','LineWidth', 2, 'color', colorMT)
       set(gca,'XTick',[-180:45:180])
        xlim([-90,90])
        ylim([0, 0.4])
        xlabel('angle (deg)', 'fontsize', 10)
        ylabel('relative counts', 'fontsize', 10)
    
   subplot(3,1,2)
   stairs(edgesEA,[hsEA,0], 'r-','LineWidth', 2, 'color', colorA)
      set(gca,'XTick',[-180:45:180])
        xlim([-90,90])
        ylim([0, 0.25])
        xlabel('angle (deg)', 'fontsize', 10)
        ylabel('relative counts', 'fontsize', 10)
    
   subplot(3,1,3)
   stairs(edgesEW,[hsEW,0],'k-','LineWidth', 2, 'color', colorW)
     set(gca,'XTick',[-180:45:180])
        xlim([-90,90])
        ylim([0, 0.2])
        xlabel('angle (deg)', 'fontsize', 10)
        ylabel('relative counts', 'fontsize', 10)

        
Plotname = strcat('FigS3_MTs_Ev_orientation');
print(fullfile(dest,strcat(Plotname,'.pdf')), '-dpdf')

