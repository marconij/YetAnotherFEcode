clear; clc; close all;

% Load solution
% The solution is a n x 3 matrix, where n is the number of elements.
% The 1st column corresponds to the element densities.
% The 2nd and 3rd columns are the x and y coordinates, respectively.
sol = readmatrix('sol.txt');

% Extract densities and coordinates
d = sol(:, 1);
coord = sol(:, 2:end);

% Specify the length unit of the coordinates
lengthUnit = [native2unicode(hex2dec({'00' 'b5'}), 'unicode'), 'm'];

% Create model
model = ComsolTO(coord, d);
model.create_model('lengthUnit', lengthUnit, 'meshSize', 3)
model.add_material_properties('2330[kg/m^3]', '148e9[Pa]', '0.23')
model.add_physics('PlaneStress', '24[um]')

% Boundary conditions
model.add_selection_box([0, 0], [1, 1e10], 'left')
model.add_fixed_boundary('box_left')

model.add_selection_box([1000, 0], [1, 1e10], 'right')
model.add_roller('box_right')

model.add_selection_box([500, 500], [0.2*1000, 1], 'load')
model.add_boundary_load('box_load', {'0[uN]', '10[uN]', '0[uN]'})

% Add study
C = model.stationary_study;
f = model.eigenfrequency_study('1[kHz]', 'eigunit', 'kHz');
volFrac = model.compute_volume_fraction;

% Save model
model.model.save(fullfile(pwd, [model.name, '.mph']));
