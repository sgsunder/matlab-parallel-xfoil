% Just like the XFOIL command, except that it allows the batch processing
% of multiple airfoils (if airfoils is given as a cell array) or multiple
% operating Reynolds numbers (if Re is given as a vector), or both.
% 
% The polar files are output in a cell array. The first index (row)
% corresponds to the specific airfoil used, and the second index (column)
% corresponds to the operating Reynolds number
% 
% As of now, N_CRIT must be input as a scalar.
% 
% displayProgressBar = if true, shows a progress bar (default: true)
% displayWindow      = if true, shows XFOIL solver window (default: false)
% pool               = define the parallel processing pool (default: gcp)
function polars = XFOILbatch(airfoils, Re, N_CRIT, AoA, displayWindow, pool)
%% Sanitize Input
% Put Airfoils in a one-element cell array, if required
if ~iscell(airfoils)
    airfoils = {airfoils};
end
% Make sure Airfoils is in a column vector
airfoils = reshape(airfoils, numel(airfoils), 1);
% Make sure Re is in a column vector
Re       = Re(:);

%% Begin Timer
batchTimer = tic;

%% Create Design Space of Airfoils, Reynolds Numbers Pairs
designSpace = zeros(length(airfoils)*length(Re), 2);
index = 1;
for ii = 1:length(airfoils)
    for jj = 1:length(Re)
        designSpace(index,:) = [ii,jj];
        index = index + 1;
    end
end

%% Parallel Execute XFOIL to solve for the 
N_RUNS = length(designSpace);
vpolar = cell(N_RUNS, 1);
vstatus = cell(N_RUNS, 1);

% Initialize Waitbar if needed
fprintf('Executing XFOIL in Batch Mode: %u threads on "%s"\n', ...
        pool.NumWorkers, pool.Cluster.Profile);
parfor_progress(N_RUNS);

% Begin Parallel Execution
thisAirfoil = cell(N_RUNS,1);
thisRe      = cell(N_RUNS,1);
for n = 1:N_RUNS
    thisAirfoil{n} = airfoils{designSpace(n,1)};
    thisRe{n}      = Re(designSpace(n,2));
end

parfor n = 1:N_RUNS
    [vpolar{n}, vstatus{n}] = XFOIL(displayWindow, thisAirfoil{n}, thisRe{n}, N_CRIT, AoA);
    parfor_progress;
end

parfor_progress(0);

% Repack Polars
polars = reshape(vpolar, length(airfoils), length(Re));

%% Print Output Results
batchTime  = toc(batchTimer);
threadTime = 0;
for ii = 1:N_RUNS
    threadTime = threadTime + vstatus{ii}.time;
end
fprintf('Completed %u runs in %-6.1f sec (%-7.1f thread-seconds). Speedup is %3.2f\n', ...
    N_RUNS, batchTime, threadTime, threadTime/batchTime);

%% Delete All Files
recycleState = recycle('off');
failedToDelete = '';
for ii = 1:N_RUNS
    filelist = vstatus{ii}.tempFiles;
    if ~isempty(filelist)
        for jj = 1:length(filelist)
            toDelete = filelist{jj};
            try
                if isdir(toDelete)
                    rmdir(toDelete,'s');
                else
                    delete(toDelete);
                end
            catch noException %#ok<NASGU>
                failedToDelete = [failedToDelete toDelete '\n']; %#ok<AGROW>
            end
        end
    end
end
recycle(recycleState);

if ~isempty(failedToDelete)
    warning('Failed to Delete the following files: \n%s', failedToDelete);
end

end