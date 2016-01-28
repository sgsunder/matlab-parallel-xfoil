% INPUTS:
% showPlot : (boolean)    Show XFOIL Plotting Window
% Airfoil  : (structure)  Airfoil Geometry, Normalized to Unit Chord
%            Airfoil.name : (string) Name of Airfoil
%            Airfoil.UX   : (vector) Upper X coordinates
%            Airfoil.LX   : (vector) Lower X coordinates
%            Airfoil.UY   : (vector) Upper Y coordinates
%            Airfoil.LY   : (vector) Lower Y coordinates
% 
function [polar, computationInfo] = XFOIL(showPlot, Airfoil, Re, N_CRIT, AoA)
%% Genereate Filepath Strings
[rootDir, ~, ~] = fileparts(mfilename('fullpath'));

% Initialize Temporary Directory
xfoilDir = tempname([rootDir filesep 'tmp']);
mkdir(xfoilDir);

% Copy XFOIL instance to Temporary Directory
if ispc
    InstallLoc = [rootDir filesep 'xfoil.exe'];
    XFOIL_Loc  = [xfoilDir filesep 'xfoil.exe'];
    [~,~] = dos(sprintf('mklink /H "%s" "%s"', XFOIL_Loc, InstallLoc));
elseif isunix
    % Use system local version of xfoil, if possible:
    [notInstalled, InstallLoc] = unix('which xfoil');
    if notInstalled ~= 0
        % If XFOIL doesn't exist locally, use the provided copy
        InstallLoc = [rootDir filesep 'xfoil'];
    end
    % Symlink the XFOIL executable to the appropriate directory
    XFOIL_Loc = [xfoilDir filesep 'xfoil'];
    unix(sprintf('ln -s %s %s', ...
        strrep(strtrim(InstallLoc), ' ', '\ '), ...
        strrep(strtrim(XFOIL_Loc) , ' ', '\ ') ...
    ));
end

%% Input and Process the Airfoil
% Ensure that All Vectors are Column Vectors
Airfoil.UX = Airfoil.UX(:);
Airfoil.UY = Airfoil.UY(:);
Airfoil.LX = Airfoil.LX(:);
Airfoil.LY = Airfoil.LY(:);

% Allocate Temporary Name
airfoil_nr   = 'airfoil.dat';
airfoil_n    = [xfoilDir filesep airfoil_nr]; % Absolute Path

% Write to Airfoil
fID = fopen(airfoil_n, 'w+');
fprintf(fID, '%s\n', Airfoil.name);
fprintf(fID, ' %f %f \n', [flipud(Airfoil.UX), flipud(Airfoil.UY)].');
fprintf(fID, ' %f %f \n', [Airfoil.LX(2:end),  Airfoil.LY(2:end) ].');
fclose(fID);

%% Reserve Polar File
polar_nr     = 'polar.out';
polar_n      = [xfoilDir filesep polar_nr];
% Mimics *nix touch command:
% fclose(fopen(polar_n, 'w+'));

%% Initialize Input File
infile_nr    = 'command.txt';
infile_n     = [xfoilDir filesep infile_nr];
fID          = fopen(infile_n, 'w+t');

% Show Plot if Asked
if ~showPlot
    fprintf(fID, 'PLOP\nG\n\n');
end

% Add Airfoil Name
fprintf(fID, 'LOAD %s\n', airfoil_nr);

% Add Filtering Step
fprintf(fID, '\nPANE\nMDES\nFILT 1.00\n');
n_filtering = 5;
for ii = 1:n_filtering
    fprintf(fID, 'EXEC\n');
end
fprintf(fID,'\nPANE\n');

% Add Processing Operation
MACH      = 0;
Vacc      = 0.01;
XTrTop    = 1;
XTrBottom = 1;

fprintf(fID, 'OPER\nVPAR\n');
fprintf(fID, 'N %1.2f\n', N_CRIT);
fprintf(fID, 'VACC %1.4f\n', Vacc);
fprintf(fID, 'XTR\n');
fprintf(fID, '%1.4f\n', XTrTop);
fprintf(fID, '%1.4f\n', XTrBottom);
fprintf(fID, '\n');
fprintf(fID, 'VISC %1.4f\n', Re);
fprintf(fID, 'MACH %1.6f\n', MACH);

% Add Starting Actions
fprintf(fID, 'ALFA %2.4f\n', AoA(1));
fprintf(fID, 'INIT\n');

% Add Open Polar File
fprintf(fID, 'PACC\n%s\n\n', polar_nr);

% Add Angle of Attack
for ii = 1:length(AoA)
    fprintf(fID, 'ALFA %2.4f\n', AoA(ii));
end

% Add Close Polar File
fprintf(fID, 'PACC\n\n');

% Add Quit Action
fprintf(fID, '\n\n\n\n\n\n\nQUIT\n');

% Close Script File
fclose(fID);

%% Execute XFOIL Application (System-Dependent)
runTimer = tic;
currentDir = pwd;
cd(xfoilDir);
if ispc % Execute DOS-Style Commands and XFOIL.exe
    arg = sprintf('xfoil.exe < %s', infile_nr);
    [exitValue, ~] = dos(arg);
elseif isunix % Execute bash commands and *nix version of XFOIL
    arg = sprintf('./xfoil < %s', infile_nr);
    [exitValue, ~] = unix(arg);
end
    
if exitValue ~= 0
    try
        delete(airfoil_n);
        delete(polar_n);
        delete(infile_n);
    catch deleteException %#ok<NASGU>
    end
    throw(MException('XFOIL:BadExitStatus', 'XFoil returned with exit status %u', exitValue));
end
    
cd(currentDir);
computationInfo.time = toc(runTimer);

%% Output Polars
% Read Polar File
fID    = fopen(polar_n, 'r');
status = fseek(fID, 429, 'bof');
data   = textscan(fID, '%f%f%f%f%f%f%f');
fclose(fID);
if isempty(data) || status ~= 0
    try
        delete(airfoil_n);
        delete(polar_n);
        delete(infile_n);
    catch deleteException %#ok<NASGU>
    end
    str = sprintf('XFOIL returned with status %u, polar cannot be read\n',status);
    throw(MException('XFOIL:BadPolarFile', str));
end

% Convert to Matrix, Remove NaNs
data   = cell2mat(data);
n = any(isnan(data), 2);
data(n,:) = [];

% Organize Polar Data
polar.alpha   = data(:,1);
polar.CL      = data(:,2);
polar.CD      = data(:,3);
polar.CDp     = data(:,4);
polar.CM      = data(:,5);
polar.Top_Xtr = data(:,6);
polar.Bot_Xtr = data(:,7);

%% Delete All Files
recycleState = recycle('off');
try
    if isunix
        %rmdir(xfoilDir, 's');
        arg = sprintf('rm -rf %s', strrep(xfoilDir, ' ', '\ '));
        [status, ~] = unix(arg);
        if status ~= 0
            throw(MException('XFOIL:LockedTempFiles','Failed to Remove Temporary Files'));
        end
        computationInfo.tempFiles = {};
    elseif ispc
        [status1,~] = dos(sprintf('DEL /F "%s"', airfoil_n));
        [status2,~] = dos(sprintf('DEL /F "%s"', polar_n));
        [status3,~] = dos(sprintf('DEL /F "%s"', infile_n));
        if (status1 || status2 || status3)
            throw(MException('XFOIL:LockedTempFiles','Failed to Remove Temporary Files'));
        end
        computationInfo.tempFiles = {XFOIL_Loc; xfoilDir};
    end
catch deleteException %#ok<NASGU>
    temp = ls(xfoilDir);
    toDelete = cell(size(temp, 1) + 1, 1);
    for ii = 1:length(toDelete)-1
        toDelete{ii} = strtrim(temp(ii,:));
    end
    toDelete{end} = {xfoilDir};
    computationInfo.tempFiles = toDelete;
end
recycle(recycleState);