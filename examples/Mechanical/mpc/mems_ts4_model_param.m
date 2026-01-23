function obj = mems_ts4_model_param(p)
% This function returns the assembly of the TS4 mems structure.
% The function receives a vector p with the parameters of the model.
% The parameters are the following:
%   p(1): dLength drive beams
%   p(2): dLength sense beams
%   p(3): dWidth drive beams
%   p(4): dWidth sense beams
%   p(5): dWidth sense beam connections
%   p(6:10):  amplitudes of the first 5 harmonics applied to the Top Drive 
%             beams to change their shape.
%   p(11:15): amplitudes of the first 5 harmonics applied to the Bottom 
%             Drive beams to change their shape.
% INPUTS:
% p: vector with the parameters of the model (default is zeros(15,1)).
% OUTPUTS:
% obj: the MultiPointConstraint object of the TS4 mems structure, with the
%       parameters set according to the input vector p.

if nargin == 0
    p = zeros(15,1);
end
if length(p)~=15
    error('p must be a vector of dimesion 15x1.')
end

deltaLength(1) = p(1);
deltaLength(2) = p(1);
deltaLength(3) = p(1);
deltaLength(4) = p(1);

deltaLength(5) = p(2);
deltaLength(6) = p(2);
deltaLength(7) = p(2);
deltaLength(8) = p(2);

model = mems_ts4_model( deltaLength , false );

model.nElements = [10*ones(4,1); repmat([10 3 10]',4,1)];
model.width = [(4+p(3))*ones(4,1); 
    repmat([3.4+p(4) 15+p(5) 3.4+p(4)]',4,1)]; % [um] beam widths
model.thickness = 24; % [um] 

model.E   = 148e3;    % [MPa] Young's modulus (polysilicon)
model.rho = 2.33e-15; % [Kg/um³] density (polysilicon)
model.nu  = 0.23;     % [-] Poisson modulus (polysilicon)

model.harmonics = [...
    p(6)   p(7)  p(8)   p(9)  p(10)
    p(6)  -p(7)  p(8)  -p(9)  p(10)
    p(11)  p(12) p(13)  p(14) p(15)
    p(11) -p(12) p(13) -p(14) p(15)
];

model.harmonics = [model.harmonics; zeros(3*4,5)];

obj = mems_assembly(model);
obj.DATA.model = model;