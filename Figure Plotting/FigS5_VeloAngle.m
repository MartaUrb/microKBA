
clear all

%% Load data
cd(fileparts(matlab.desktop.editor.getActiveFilename));
currentpath = cd();
 idcs   = strfind(currentpath,'/');
 motherpath = currentpath(1:idcs(end)-1);
 
sourcepath = fullfile(motherpath,  'Data', 'Events', 'Events_compiled');
dest = fullfile(motherpath, 'Figures');
cd(sourcepath)

%%Load compiled event data
sourcename = 'KBA_allExperiment_Events_Sel3.txt';

dataimp = importdata(fullfile(sourcepath,sourcename));
data = dataimp.data;
headers = dataimp.textdata;
headers = split(headers).';
clear dataimp
 
 
 %% Plot angle resolved for selected force bin - 6 bins
   ind_F = 14;
   binC = [0.5:0.5:3.25];
   f = figure
for  ii = 1:length(binC);
    x = binC(ii);
   pm = 0.25;
   xL = x - pm;
   xH = x +pm;
   dataTemp = data(data(:,ind_F)>=xL&data(:,ind_F)<=xH,:);
   nel = length(dataTemp);

   clear valuesE
   
barN = 31;
ratesE = dataTemp(:,10);
valuesE = dataTemp(:,11);
edgesE = linspace (-180, 180, (barN+1));
edgesE(1)= -180.001;
edgesE(length(edgesE))= 180.001;

bin = edgesE(2)-edgesE(1);
binH = bin./2;
centersE = edgesE(1:end-1)+binH;

    NE = length(edgesE);
    hsE = zeros(length(edgesE)-1,1);
    
    for n = 1:NE-1
        ind1 = find(valuesE >= edgesE(n) & valuesE < edgesE(n+1));
        indices1{n}=ind1;
        rateValues{n} = dataTemp(indices1{n},10);
        meanRates(n,1) = mean(rateValues{n});
        medianRates(n,1) = median(rateValues{n});
        stdev(n,1) = std(rateValues{n});
        SEM(n,1) = stdev(n,1)/sqrt(length(rateValues{n}));
        hsE(n,1) = length(ind1)/length(valuesE);
        nB(n,1) = length(ind1);  
        perc1(n,1) = prctile(rateValues{n},25);
        perc3(n,1) = prctile(rateValues{n},75);
        meanAD(n,1) = mad(rateValues{n});
        medAD(n,1) = mad(rateValues{n},1);
    end
    
%hsEcorrN = hsEcorr./sum(hsEcorr);

% rate(angle) scatter DOUBLE
  set(0,'DefaultAxesFontName', 'Arial')
    set(0,'DefaultAxesFontSize', 10)

set(0,'DefaultFigureVisible','on');
patch1 = [-180 0; -135 0; -135 1.2; -180 1.2];
patch2 = [-45 0; 45 0; 45 1.2; -45 1.2];
patch3 = [180 0; 135 0; 135 1.2; 180 1.2];

subplot(3,2,ii)
%h = figure
cA = [239, 76, 67]./255;
cW = [85, 61, 151]./255;

patch(patch1(:,1), patch1(:,2), cW, 'facealpha', 0.3, 'edgecolor', 'none')
hold on
patch(patch2(:,1), patch2(:,2), cA, 'facealpha', 0.3, 'edgecolor', 'none')
hold on
patch(patch3(:,1), patch3(:,2), cW, 'facealpha', 0.3, 'edgecolor', 'none')
hold on
scatter(dataTemp(:,11), dataTemp(:,10), 15,'o', 'fill','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.1,'markerfacecolor', [0.4 0.4 0.4 ], 'markeredgecolor', 'none')
hold on
scatter(centersE,meanRates, 8, 'fill', 'markerfacecolor', 'k', 'markeredgecolor', 'none')
hold on
errorbar(centersE,meanRates,stdev,'.','CapSize',3, 'Marker', 'none', 'LineWidth', 0.5, 'color', 'k')

box on 
grid on
ylim([0 1.2])
xlim([-180 180])
dim = [.2 .8 .3 .05];
ax = gca;
ax.FontSize = 8;
ax.XTick = [-180:45:180]
ax.YTick = [0:0.25:1]

str1 = strcat(num2str(x), {' '},'±', {' '}, num2str(pm), {' '}, 'pN, {\it n} =',{' '}, num2str(nel));
str = str1{1,1};
%annotation('textbox',dim,'String',str1,'linestyle', 'none','FitBoxToText','on', 'Fontsize',14);

title(str1, 'FontWeight', 'normal')
xlabel ('angle (°)', 'FontSize', 10)
ylabel('velocity (um s^{-1})', 'FontSize', 10)

end


set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 1.144.*14.37 1.113.*13.84]); 

Plotname = 'FigS5_VeloAngle';

 print(fullfile(dest,strcat(Plotname,'.pdf')), '-dpdf')
 