clc

[filepath, name, ext] = fileparts(mfilename('fullpath'));

addpath(genpath(fullfile(filepath, 'src')));
addpath(genpath(fullfile(filepath, 'external')));
addpath(genpath(fullfile(filepath, 'examples', 'Meshes')));

disp('              _____ _____     ')
disp('  _   _  __ _|  ___| ____|___ ')
disp(' | | | |/ _` | |_  |  _| / __|')
disp(' | |_| | (_| |  _| | |__| (__ ')
disp('  \__, |\__,_|_|   |_____\___|')
disp('  |___/       YetAnotherFEcode')
fprintf('\n\n')