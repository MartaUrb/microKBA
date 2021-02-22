
clear all

%% list files 
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

%% Building histogram of veocities from bin centers
clear setcolsel

anglesel = 45;
ind_ang = [find(abs(data(:,11))<=anglesel); find(abs(data(:,11))>=(180-anglesel))];
setcolsel = data(ind_ang,:);

ind_F = 15;
emax = 3.45;
emin = 0.3;
step = 0.3;
FminWE = max(setcolsel(setcolsel(:,ind_F)<0,15) );
FminAE = min(setcolsel(setcolsel(:,ind_F)>0,15) );

clear barN values rates edgesA edgesW N hsAE hsAE2 hsWE hsWE2 Range indAE indWE DmedAE DmedWE DavAE DavWE centersAE centersWE prc1WE prc2WE prc1AE prc2AE
barN = 23;
values = setcolsel(:,ind_F);
rates = setcolsel(:,10);
edgesW = [-emax:step: -emin];
edgesA = fliplr(edgesW*-1);
    
    N = length(edgesW);
    hsWE = zeros(size(edgesW)-1);
    hsAE = zeros(size(edgesA)-1);

    Range = abs(edgesW(2)-edgesW(1));
    
    for n = 1:N-1
        indWE{n} = find(values >= edgesW(n) & values < edgesW(n+1));        
        indAE{n} = find(values >= edgesA(n) & values < edgesA(n+1));

        hsWE(n) = nanmedian(rates(indWE{n}));
        prc1WE(n) = prctile(rates(indWE{n}),25);
        prc2WE(n) = prctile(rates(indWE{n}),75);
        hsWE2(n) = nanmean(rates(indWE{n}));
        DmedWE(n) = mad(rates(indWE{n}),1);
        DavWE(n) = nanstd(rates(indWE{n}));
        DavWE2(n) = nanstd(rates(indWE{n}))./sqrt(length(indWE{n}));
        centersWE(n) = edgesW(n)+Range./2;
        lengthWE (n) = length(indWE{n});
        
        hsAE(n) = nanmedian(rates(indAE{n}));
        prc1AE(n) = prctile(rates(indAE{n}),25);
        prc2AE(n) = prctile(rates(indAE{n}),75);
        hsAE2(n) = nanmean(rates(indAE{n}));
        DmedAE(n) = mad(rates(indAE{n}),1);
        DavAE(n) = nanstd(rates(indAE{n}));
        DavAE2(n) = nanstd(rates(indAE{n}))./sqrt(length(indAE{n}));
        centersAE(n) = edgesA(n)+Range./2;
        lengthAE (n) = length(indAE{n});

    end
    centersAE(1) = centersAE(2)-Range;

 %% Plot F-V Mean+/-SD
  set(0,'DefaultAxesFontName', 'Arial')
  set(0,'DefaultAxesFontSize', 10)
    
    colorR = [239, 76, 67]./255;
    colorB = [85, 61, 151]./255;
    
    f = figure
xline(0)
hold on
      scatter(setcolsel(:,15), setcolsel(:,10), 10,'o', 'fill','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1,'markerfacecolor', [0.4 0.4 0.4 ], 'markeredgecolor', 'none')
hold on
    errorbar(0, 0.662, 0.072,'ok', 'linewidth',1, 'MarkerSize',5,'markerfacecolor', [0 0 0], 'markeredgecolor', 'none')

hold on

    errorbar(centersWE, hsWE2,DavWE, '.', 'linewidth',1, 'color', colorB)
    hold on
    errorbar(centersAE, hsAE2,DavAE, '.', 'linewidth',1, 'color', colorR)
    hold on
    scatter(centersWE, hsWE2, 20,'o', 'fill', 'markerfacecolor', colorB, 'markeredgecolor', 'none')
    hold on
    scatter(centersAE, hsAE2, 20,'o', 'fill', 'markerfacecolor', colorR, 'markeredgecolor', 'none')

   xlim([-4,4])
   ylim([0, 1.3])
   ax = gca;
    ax.FontSize = 8;
    ax.XTick = [-6:1:6]
   xlabel ('force (pN)', 'FontSize',10)
   ylabel('velocity (um s^{-1})', 'FOntSize',10)
   box on 
   grid off
      
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 7.67 6.44]); 

Plotname = 'Fig4a_FVcurve';

 print(fullfile(dest,strcat(Plotname,'.png')), '-dpng')
%% 
 headersBin = {'Force', 'Velocity', 'SD', 'SEM','n'};
 dataBin = [[centersWE.';centersAE.'], [hsWE2.';hsAE2.'], [DavWE.';DavAE.'],[DavWE2.';DavAE2.'],[lengthWE.';lengthAE.']];
 
 % Save Bin Data
format = char();
for ii = 1:length(headersBin)
    format = strcat(format, '%s\t');
end
format = strcat(format, '\n');

fullname = fullfile(dest, strcat(Plotname,'.txt'));

fid = fopen(fullname, 'wt');
fprintf(fid,  format, headersBin{:});  % header
fclose(fid);
dlmwrite(fullname, dataBin, 'delimiter', '\t','precision', 16, '-append')
