clear all; close all;

%% Define Parameters
Re              = logspace(6,12,30);
N_CRIT          = 6;
AoA             = -2.5:0.25:20;

%% Generate Airfoil List
afcode          = {'0012', '2412', '4415', '6409', '23012', '25112'};
airfoils        = cell(size(afcode));
for ii = 1:length(afcode)
   airfoils{ii} = generateNACAairfoil(afcode{ii}, 100);
end

%% Load All Airfoils
% [rootDir, ~, ~] = fileparts(mfilename('fullpath'));
% afdir   = [rootDir filesep 'airfoils'];
% affiles = dir(afdir);
% affiles = affiles(3:end);
% airfoils = cell(size(affiles));
% for ii = 1:length(affiles)
%     airfoils{ii} = scanAirfoil([afdir filesep affiles(ii).name]);
% end

%% Solve Airfoils
tic;
pool            = gcp;
polars          = XFOILbatch(airfoils, Re, N_CRIT, AoA, false, pool);
delete(pool);
toc;

%% Save Results
%save('exampleResults.mat', 'polars', 'airfoils', 'Re', 'N_CRIT', 'AoA');

%% Plot Results
cmap = hsv(size(polars,1));
h    = figure;
plot(AoA, 2*pi^2*AoA/180, 'k--');
xlabel('Angle of Attack (deg)');
ylabel('Coefficient of Lift');
hold on;
for ii = 1:size(polars,1)
    for jj = 1:size(polars,2)
        plot(polars{ii,jj}.alpha, polars{ii,jj}.CL, '-', 'Color', cmap(ii,:));
    end
end
hold off;