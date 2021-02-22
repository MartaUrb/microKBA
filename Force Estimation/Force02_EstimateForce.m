clear all

%% list files 
cd(fileparts(matlab.desktop.editor.getActiveFilename));
currentpath = cd();
 idcs   = strfind(currentpath,'/');
 motherpath = currentpath(1:idcs(end)-1);
 
sourcepath = fullfile(motherpath,  'Data', 'Events', 'Events_byExp_PlusMSDy');
dest = fullfile(motherpath, 'Data', 'Events', 'Events_byExp_PlusForce');

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
    clear dataimp
end

ind_msd = find(contains(headers,'msd'));

%flowrates used in the different experiemnts
fr = [20;20;40;10;30;10;30;50];

%%% constants for Force estimation
kB = 1.3806485279*10^-23; %Boltzmann constant
TC = 23; %23C
T = TC + 273.15; %temperature in Kelvin
kBT = kB*T;
P = 50*10^-9; %Persistance length of dsDNA [m]
b = kBT/P; %[N]
l0 = 16.2*10^-6; %tether contour length 
W = 0.15; %image integration time in seconds
R = 0.54*10^-6; %radius of the bead
mu = 0.9321.*10.^-3; %[Pa*s] dynamic viscosity at 23C from https://wiki.anton-paar.com/en/water/
gam = 6.*pi.*mu.*R; %friction factor
a = gam/W;

%%%%% further parameters for bead height estimation %%%%%
cf = 10.^-9./60; %conversion factor for flowrate from ul/min -> m3/s
% channel dimensions in m
h = 100.*10.^-6; % channel height[m]
w = 3.*10.^-3; % channel width [m]
%Radius of the bead - diameter 1.08 um
R = 0.54.*10.^-6; %[m]
%Magnetic force
Fmag = 0.1.*10.^-12; %[N]


%% Solve simple set of equations - plus Faxen ut flow velocity after AJP 2011
Flowrate_unq = unique(fr);

for ii = 1:length(Flowrate_unq)
clear z l
syms z l

%%%%% Stoke's drag for Poiseuille flow %%%%%
vmax = 3.*Flowrate_unq(ii)*cf./(2.*w.*h);
Fdrag = 6.*pi.*mu.*R.*4.*vmax.*(z./h).*(1-z./h);
% correction factor for surface proximity, for movement in the direction
% parallel to the surface, after Faxen 1923:
lmbd_par = (1- 9./16 .* R./z + 1./8 .* (R./z).^3 - 45./256.*(R./z).^4 - 1./16.*(R./z).^5 ).^(-1);
Fdrag_corr = lmbd_par.*Fdrag;

%%%%% Entropic restoring force of the DNA tether %%%%%
Ftether = b.*(0.25./(1-l./l0).^(2)-0.25+l./l0);

sin_th = z/(l+R);
cos_th = sqrt(1-sin_th^2);
% with tetha being the angle between the tether and the surface

f1 =  Fdrag_corr - cos_th*Ftether;
f2 =  Fmag - sin_th*Ftether;

initial_guess = [0.54 16.2; 0 16.2] *10^-6;

[z_sol, l_sol] = vpasolve([f1,f2], [z,l], initial_guess);

z_sol_um(ii,1) = z_sol*10^6;
l_sol_um(ii,1) = l_sol*10^6;
lmbd_perp_sol(ii,1) = [1- 9./16 .* R./z_sol + 1./8 .* (R./z_sol).^3 - 45./256.*(R./z_sol).^4 - 1./16.*(R./z_sol).^5 ].^(-1);

Fdrag_sol(ii,1) =  double(6.*pi.*mu.*R.*4.*vmax.*(z_sol./h).*(1-z_sol./h)*10^12);
Fdrag_sol_corr(ii,1) = double(lmbd_perp_sol(ii,1).*Fdrag_sol(ii,1));
Ftether_sol(ii,1) = double(b.*(0.25.*(1-l_sol./l0).^(-2)-0.25+l_sol./l0)*10^12);
end

sinus_th= double(z_sol_um)./(double(l_sol_um)+R*10^6);
cosinus_th = sqrt(1-sinus_th.^2);
thetha = rad2deg(asin(sinus_th));
posx = cosinus_th.*(double(l_sol_um)+R*10^6);

resultsFB1 = [Flowrate_unq , double(Ftether_sol), double(Fdrag_sol_corr), posx, double(z_sol_um), double(l_sol_um)+R*10^6, double(l_sol_um), thetha];
zposFB = double(z_sol_um).*10^-6;
lmbd_par = (1- 9./16 .* R./zposFB + 1./8 .* (R./zposFB).^3 - 45./256.*(R./zposFB).^4 - 1./16.*(R./zposFB).^5 ).^(-1);


 %% Calculate Force with msd corrected for blurring & gamma corrected for surface proximity
% blurring after Wong und Harvelson 2006
% gamma corrected using z position from force-balance estimation
syms F l
bounds = [0, 50*10^-12; 0.1*10^-12, 16.2*10^-6]; %upper and lower boundaries for F and l

for aa = 1:length(dataESel)
msdTemp = dataESel{aa}(:,ind_msd);
%msdTemp = msdTemp(~isnan(msdTemp));
zTemp = zposFB(find(Flowrate_unq ==fr(aa)));
lmbd_par_temp = (1- 9./16 .* R./zTemp + 1./8 .* (R./zTemp).^3 - 45./256.*(R./zTemp).^4 - 1./16.*(R./zTemp).^5 ).^(-1);
atemp = a.*lmbd_par_temp;
for ii = 1:length(msdTemp)
SB(ii) = vpasolve([F*P/kBT == 1/4*(1-l/l0 )^(-2)-1/4+l/l0, ...
                   2*atemp*l/F - 2*(atemp^2)*(l/F)^2*(1-exp(-F/(atemp*l))) == F*msdTemp(ii)/(kBT*l)],...
                   bounds);

LxB(ii,1) = double(SB(ii).l).*10^6.'; % the estimated tether length in um
FxB(ii,1) = double(SB(ii).F).*10^12.'; % the estimated force in pN
end

dataESel{aa}(:,13) = LxB;
dataESel{aa}(:,14) = FxB;
clear SB LxB FxB msdTemp atemp lmbd_par_temp
aa
end

headers{1,13} = 'L_um';
headers{1,14} = 'F_pN_abs';

%% F not absolute 
% Force is defined as negative for molecules stepping against the flow direction
% stepping angle > 90 deg

for aa = 1:length(dataESel)
      dataESel{aa}(:,15) =  dataESel{aa}(:,14);
      dataESel{aa}(find(abs(dataESel{aa}(:,11))>90),15) =    -1.* dataESel{aa}(find(abs(dataESel{aa}(:,11))>90),15) ;
end
headers{1,15} = 'F_pN';
dataESel_col = vertcat(dataESel{:});

%% Save segments with Force, 1 file per experiment
format = char();
for ii = 1:length(headers)
    format = strcat(format, '%s\t');
end
format = strcat(format, '\n');

for ii =1:length(dataESel)
fullname = fullfile(dest, strcat(exp_name{ii},'Force.txt'));
fid = fopen(fullname, 'wt');
fprintf(fid,  format, headers{:});  % header
fclose(fid);
dlmwrite(fullname, dataESel{ii}, 'delimiter', '\t','precision', 16, '-append')
clear fullname fid 
end