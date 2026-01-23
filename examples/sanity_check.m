%% Run All MFiles
% This script runs all .m files in the current folder and subfolders.
% It will try/catch each file and print errors without stopping execution.
% .m functions will be skipped.

clear; close all; clc

% Get list of all .m files recursively
fileList1 = dir(fullfile(pwd, '**', '*.m'));
fileList2 = dir(fullfile(pwd, '**', '*.mlx'));
fileList = [fileList1; fileList2];

% Store errors in a cell array
errors = {};

for k = 1:numel(fileList)
    filePath = fullfile(fileList(k).folder, fileList(k).name);
    [~, fileName] = fileparts(filePath);

    if strcmpi(fileName,'sanity_check')==1
        continue;
    end

    % if ~isScriptFile(filePath)
    %     fprintf('Skipping non-script file: %s\n', filePath);
    %     continue;
    % end
    
    fprintf('Running: %s ...\n', filePath);
    
    try
        runFileSafely(filePath);
    catch ME
        fprintf(2, 'Error in %s: %s\n', filePath, ME.message);
        errors{end+1} = struct('file', filePath, 'error', ME); %#ok<SAGROW>
    end
end
close all;

% Print summary of errors
if isempty(errors)
    fprintf('\nAll files executed successfully.\n');
else
    fprintf('\nSummary of errors:\n');
    for i = 1:numel(errors)
        fprintf(2, 'File: %s\nError: %s\n\n', errors{i}.file, errors{i}.error.message);
    end
end

% --- Save error report to a text file ---
if ~isempty(errors)
    % Create filename with date/time
    timestamp = char(datetime('now','format','yyyy-MM-dd_HHmmss'));
    logFile = fullfile(pwd, ['ErrorReport_' timestamp '.txt']);
    
    fid = fopen(logFile,'w');
    if fid == -1
        warning('Could not create log file.');
    else
        fprintf(fid,'Error Report (%s)\n\n', datestr(now));
        for i = 1:numel(errors)
            % fprintf(fid,'File: %s\nError: %s\n\n', ...
            %     errors{i}.file, errors{i}.error.message);
            % If you want full stack trace:
            fprintf(fid, '%s\n\n', getReport(errors{i}.error));
        end
        fclose(fid);
        fprintf('\nError report saved to: %s\n', logFile);
    end
end


%% functions

function runFileSafely(filePath)
    evalc(sprintf('run(''%s'')', filePath));
end

function isScript = isScriptFile(filePath)
    % Open file
    fid = fopen(filePath,'r');
    isScript = true; % assume script unless proven otherwise
    
    if fid ~= -1
        while ~feof(fid)
            line = strtrim(fgetl(fid));
            % Skip empty lines and comments
            if isempty(line) || startsWith(line,'%')
                continue;
            end
            % First meaningful line found
            if startsWith(line, 'function') || startsWith(line,'classdef')
                isScript = false;
            end
            break;
        end
        fclose(fid);
    end
end