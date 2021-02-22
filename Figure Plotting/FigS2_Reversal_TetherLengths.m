clear all

%% list files 
cd(fileparts(matlab.desktop.editor.getActiveFilename));
currentpath = cd();
 idcs   = strfind(currentpath,'/');
 motherpath = currentpath(1:idcs(end)-1);
 
sourcepath = fullfile(motherpath,  'Figures', 'FigS2_Tether_Reversal');
dest = fullfile(motherpath, 'Figures', 'FigS2_Tether_Reversal');
cd(sourcepath)

%%for the histograms in figure 3
sourcename = 'reversal_09122013_60_73_78_composite_length mes.txt';

%%% for selecting any experiemnt:
%%[sourcename,~]=uigetfile('*.txt');


%open data
alldata = importdata(sourcename);
data = alldata.data;
headers = alldata.colheaders;

ind_width = find(contains(headers,'Width'));
ind_len = find(contains(headers,'Length'));

lengthsum = (((data(:,ind_len)*1.6)-2.8)/2);


%% hist relative counts
binN = 50;
binCenters = linspace (2.5, 97.5, binN)
    [counts binCenters] = hist(lengthsum, binN, binCenters); 
       relcounts = counts / sum(counts);
       Gau = fit(binCenters',relcounts.','gauss1')
       xVal = [0:0.1:30];
       yVal = Gau.a1.*exp(-((xVal-Gau.b1)/Gau.c1).^2);
   
h = figure
    patch([Gau.b1-3.*Gau.c1, Gau.b1+3.*Gau.c1,Gau.b1+3.*Gau.c1,Gau.b1-3.*Gau.c1], [0 0 0.5 0.5], [0 0 0], 'facealpha', 0.3)
    hold on
         c = findobj(gca,'Type','patch')
        set(c,'FaceColor',[0.5 0.5 0.5],'EdgeColor','none')
    bar(binCenters, relcounts, 'BarWidth', 1, 'FaceColor',[224 184 220]./255,'EdgeColor','none')
    hold on
    plot(xVal, yVal, 'linewidth', 3, 'color', [170 63 155]./255);
    hold on
    xline(16.2, '--','linewidth',3, 'color', [0 151 20]./255);
    hold on
    xline(Gau.b1, 'linewidth',3, 'color', [170 63 155]./255);
    

ylim ([0, 0.3])
xlim([0, 25])
xlabel ('projected tether length (um)', 'fontsize', 12)
ylabel ('relative counts', 'fontsize', 12)
set(gca,'FontSize',10)

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 12 7]); 
box on
Plotname = strcat('FigS2d_Histogram_tetherLength', '_binN', num2str(binN));

 print(fullfile(dest, strcat(Plotname,'.png')), '-dpng', '-r300')
 
 