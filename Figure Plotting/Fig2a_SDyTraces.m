clear all

%%% list files 
cd(fileparts(matlab.desktop.editor.getActiveFilename));
currentpath = cd();
 idcs   = strfind(currentpath,'/');
 motherpath = currentpath(1:idcs(end)-1);
 dest = fullfile(motherpath, 'Figures');

sourcepath = fullfile(motherpath, 'Data', 'Traces');
segmentspath =  fullfile(motherpath,  'Data', 'Events', 'Events_byExp_PlusForce');

cd(sourcepath)
dirinfo = dir();
dirinfo = dirinfo([dirinfo.isdir]);
dirinfo(find(contains({dirinfo.name}, '.'))) = [];

% colors
cc(1,:) = [57 177 206]./255; %blue
cc(2,:) = [123 33 127]./255; %purple
cc(3,:) = [156 197 69]./255; %lime green
cc(4,:) = [232 74 30]./255;  %koral red

%% open selected tracjestories
cd(sourcepath)
exp = {'150123_fc2_10ulmin', '141111_fc2_20ulmin', '150128_fc2_30ulmin', '150123_fc1_40ulmin'};
region = {'x2y1','x3y2','x1y1','x3y1'};
traj = [317, 727, 214, 656];
offset = [70, 70, 70, 70];
framest = [2427, 563, 643, 573];
for ii = 1:length(exp)
   tracedata = importdata(fullfile(sourcepath, exp{ii}, strcat(exp{ii}(1:end-8), '_', region{ii}, '_SDy_results.txt')));
%   startend =  importdata(fullfile(sourcepath, exp{ii}, strcat(exp{ii}(1:end-8), '_', region{ii}, '_SDy_startend.txt')));
%   times(ii,:)= startend.data(find(startend.data(:,2)==traj(ii)),:);
   indtrace = find(tracedata.data(:,10)==traj(ii) & tracedata.data(:,9)>=framest(ii));
  % indst = find(tracedata.data(indtrace,13)==times(ii,3));
  % indend = find(tracedata.data(indtrace,13)==times(ii,4));
   posY{ii} = tracedata.data(indtrace,12);
   timeY{ii} = tracedata.data(indtrace,13);
   timeYNew{ii} = timeY{ii}-timeY{ii}(offset(ii));
   msd_um2{ii} = mean(posY{ii}.^2)-mean(posY{ii}).^2;
   av{ii} = mean(posY{ii});
   posYNew{ii} = posY{ii}-av{ii};
   Sdy(ii,1) = sqrt(msd_um2{ii});
   clear tracedata startend indst indend
   trajname{ii} = strcat(exp{ii}, '_', region{ii}, '_traj',num2str(traj(ii)));
end

%% Prepare data for ploting
% Running average
step = 100; %step size for rolling average 

%calculate rolling st dev 
for ii = 1:length(posY)
        [fBeadsAvY{ii}, fBeadsStdY{ii}] = fctRunningAverage(posYNew{ii}, step);
end

 % prepare histpgram data
 clear barN values edges N hs Range centers
for ii = 1:4
barN = 20
values = posYNew{ii}(70:403);
maxV = max(values);
minV = min(values);
edges = linspace(-1,1,barN);
    
    N = length(edges);
    hs{ii} = zeros(size(edges)-1);
    Range = edges(2)-edges(1);
    
    for n = 1:N-1
        ind = find(values >= edges(n) & values < edges(n+1));
        hs{ii}(n) = length(ind)/length(values);
        centers{ii}(n) = edges(n)+Range./2;
    end
end

%% Plot
  
set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 6)

set(0,'defaultLineLineWidth',0.5); 

f = figure
f.Units = 'centimeters';
f.Position = [1.5 1.5 12 14];

for ii=1:4
subplot(4,2,2*ii-1, 'Position',[0.1300 0.1100+(4-ii).*0.2201 0.44 0.1577])
  clear x y yy y1 y2
  er = repmat(Sdy(ii),length(timeYNew{ii}),1);

    x = timeYNew{ii}.';
    y1=(fBeadsAvY{ii}-er).';                      %#create first curve
    y2=(fBeadsAvY{ii}+er).';                   %#create second curve
    y=[x,fliplr(x)];                %#create continuous x value array for plotting
    yy=[y1,fliplr(y2)];              %#create y values for out and then back

    h = fill(y,yy,'b');
    set(h, 'FaceColor', cc(ii,:),'FaceAlpha', 0.2, 'LineStyle', 'none')
    hold on
    er2 = fBeadsAvY{ii}(1,1);

    plot(timeYNew{ii}(:), posYNew{ii}(:), 'linewidth', 1, 'color', cc(ii,:))
    hold on
    plot(timeYNew{ii}(:), fBeadsAvY{ii}(:), 'linewidth', 2, 'color', cc(ii,:))
    hold on
    
    set(gca,'YTick',[-0.25, 0, 0.25])
    set(gca,'XTick',[0:10:50])

    ylim([-0.5 0.5])
    xlim ([0 50])
    ylabel('y (um)')
    if ii == 4
        xlabel('time (s)')
    end
    
s2(ii) = subplot(4,2,2*ii,'Position',[0.62 0.1100+(4-ii).*0.2201 0.12 0.1577])

   g = barh(centers{ii},hs{ii}, 'FaceColor',cc(ii,:),'EdgeColor','black')
   set(gca,'yticklabel', '')
   hold on
     
  c = findobj(gca,'Type','patch')
    set(c,'FaceColor',cc(ii,:),'EdgeColor','none')
    set(g,'barwidth',1) 
    set(gca,'XTick',[0:0.2:0.4])
    ylim([-0.5, 0.5])
    xlim([0, 0.5]) 
     if ii == 4
        xlabel('counts')
    end
  end
    
    print(fullfile(dest,strcat('Fig2a_traces.png')), '-dpng', '-r300')

    
%% Save Source Data

headersT = {'time_ms','y_um'};
indS = 70;
indE = 404;

for jj = 1:length(exp) 
dataT = [timeYNew{jj}(indS:indE),posYNew{jj}(indS:indE)]; 

format = char();
for ii = 1:length(headersT)
    format = strcat(format, '%s\t');
end
format = strcat(format, '\n');

fullname = fullfile(dest, strcat('Fig2a_SourceData_tr', num2str(jj),'_',trajname{jj},'.txt'));

fid = fopen(fullname, 'wt');
fprintf(fid,  format, headersT{:});  % header
fclose(fid);
dlmwrite(fullname, dataT, 'delimiter', '\t','precision', 16, '-append')
end

%% Save histogram data
headershsT = {'bin_centers','rel_counts_hs1','rel_counts_hs2','rel_counts_hs3','rel_counts_hs4'};

hsT = centers{1,1}(:);

for ii = 1:length(exp)
    hsT = [hsT, hs{ii}(:)];
end

format = char();
for ii = 1:length(headershsT)
    format = strcat(format, '%s\t');
end
format = strcat(format, '\n');

fullname = fullfile(dest, strcat('Fig2a_SourceData_histograms','.txt'));

fid = fopen(fullname, 'wt');
fprintf(fid,  format, headershsT{:});  % header
fclose(fid);
dlmwrite(fullname, hsT, 'delimiter', '\t','precision', 16, '-append')
