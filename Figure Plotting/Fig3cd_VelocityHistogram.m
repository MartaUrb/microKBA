 clear all
    
%% data opening and structuring

cd(fileparts(matlab.desktop.editor.getActiveFilename));
currentpath = cd();
 idcs   = strfind(currentpath,'/');
 motherpath = currentpath(1:idcs(end)-1);
 
sourcepath = fullfile(motherpath,  'Data', 'Events', 'Events_byExp_PlusForce');
dest = fullfile(motherpath, 'Figures');
cd(sourcepath)

%%for the histograms in figure 3
sourcename = '141112_fc2_segments_sel3_t92_Force.txt';

%%% for selecting any experiemnt:
%%[sourcename,~]=uigetfile('*.txt');


dataimp = importdata(sourcename);
headers = dataimp.textdata;
headers = split(headers);

dataAll = dataimp.data;
ind_ang = find(contains(headers,'ang'));
ind_vel = find(contains(headers,'rate'));

%find indeces of events with/against the flow +/-45Deg
dataSelAgainst = dataAll(abs(dataAll(:,ind_ang))<45,:);
dataSelWith = dataAll(abs(dataAll(:,ind_ang))>135,:);

meanVA = mean(dataSelAgainst(:,ind_vel));
meanVW = mean(dataSelWith(:,ind_vel));
stdVA = std(dataSelAgainst(:,ind_vel));
stdVW = std(dataSelWith(:,ind_vel));
sEMVA = std(dataSelAgainst(:,ind_vel))./sqrt(length(dataSelAgainst));
sEMVW = std(dataSelWith(:,ind_vel))./sqrt(length(dataSelWith));

%define colors 
ccA = [233 78 68]./255;
ccW = [85 63 146]./255;
ccFrame = [51 51 51]./255;

  %% building histogram from defined edges AGAINST/WITH the flow
clear edges  indA hsEA centersEA indW hsEW centersEW
edges = [0:0.05:2];
valuesEA = dataSelAgainst(:,ind_vel);
pdEA = fitdist(valuesEA, 'normal');

    NEA = length(edges);
    hsEA = zeros(size(edges)-1);
    RangeA = edges(2)-edges(1);
    
    for n = 1:NEA-1;
        indA{n} = find(valuesEA >= edges(n) & valuesEA < edges(n+1));
        hsEA(n) = length(indA{n})/length(dataSelAgainst);
        centersEA (n) = edges(n)+RangeA./2;
    end
Gau1A = fit(centersEA.',hsEA.','gauss1');

valuesEW = dataSelWith(:,ind_vel);
pdEW = fitdist(valuesEW, 'normal');
    NEW = length(edges);
    hsEW = zeros(size(edges)-1);
    RangeW = edges(2)-edges(1);
    
    for n = 1:NEW-1;
        indW = find(valuesEW >= edges(n) & valuesEW < edges(n+1));
        hsEW(n) = length(indW)/length(dataSelWith);
        centersEW (n) = edges(n)+RangeW./2;
    end
Gau1W = fit(centersEW.',hsEW.','gauss1');

    %% plot 
   set(0,'DefaultAxesFontSize', 8)
   set(0,'DefaultAxesFontName', 'Arial')

   f = figure
   f.Units = 'centimeters';
   f.PaperPosition = [0 0 20 5.4]; 
    
   xval = [0:0.001:1];
   yvalA = Gau1A.a1.*exp(-((xval-Gau1A.b1)/Gau1A.c1).^2) ;
   yvalW = Gau1W.a1.*exp(-((xval-Gau1W.b1)/Gau1W.c1).^2) ;
   
   subplot(1,2,1)
   g1 = bar(centersEA,hsEA, 'facealpha', 0.7)
   hold on
   plot(xval, yvalA, '-','linewidth', 2, 'color', ccA);
   hold on
   xline(mean(valuesEA), '--', 'linewidth',3, 'color', ccA)
    set(gca,'FontSize', 10)

    set(g1,'FaceColor',ccA,'EdgeColor',ccFrame)
    set(g1,'barwidth',1) 
    %set(gca,'XTick',linspace(-180,180,10))
    %set(gca, 'XTickLabel',linspace(-180,180,10))
    xlim([0, 1])
 %  ylim([0, max(hsEA)+0.1]) 
    ylim([0 0.2])
    xlabel('Velocity (um s^{–1})', 'fontsize', 10)
    ylabel('Relative Counts', 'fontsize', 10)
hold on


   subplot(1,2,2)
   g2 = bar(centersEW,hsEW,'facealpha', 0.7)    
    set(g2,'FaceColor',ccW,'EdgeColor',ccFrame)
    set(g2,'barwidth',1) 
     hold on
   plot(xval, yvalW, '-','linewidth', 2, 'color', ccW);
   hold on
   xline(mean(valuesEW), '--', 'linewidth',3, 'color', ccW)
   
    xlim([0, 1])
    ylim([0, 0.16]) 
    xlabel('Velocity (um s^{–1})', 'fontsize', 10)
    ylabel('Relative Counts', 'fontsize', 10)


   
    Plotname = strcat('Fig3cd_','velocity_histograms_',sourcename(1:10))
    print(fullfile(dest, strcat(Plotname,'.png')), '-dpng', '-r300')

    
    %% inspect plot data
   histoA =  [centersEA.',hsEA.']
   histoW = [centersEW.',hsEW.']

    Gau1A
    Gau1W
    
    
%% Save Bin Data
format = char();

headersBins = {'centers','rel_counts'}

for ii = 1:length(headersBins)
    format = strcat(format, '%s\t');
end
format = strcat(format, '\n');

fullnameA = fullfile(dest, strcat('Fig3c_Velo_BinData_Against','.txt'));
fidA = fopen(fullnameA, 'wt');
fprintf(fidA,  format, headersBins{:});  % header
fclose(fidA);
dlmwrite(fullnameA, histoA, 'delimiter', '\t','precision', 16, '-append')

fullnameW = fullfile(dest, strcat('Fig3d_Velo_BinData_With','.txt'));
fidW = fopen(fullnameW, 'wt');
fprintf(fidW,  format, headersBins{:});  % header
fclose(fidW);
dlmwrite(fullnameW, histoW, 'delimiter', '\t','precision', 16, '-append')

