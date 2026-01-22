% damping_rayleigh
%
% Synthax:
% [alpha, beta] = damping_rayleigh(dtype, Q, f0, n)
%
% Description: computes the alpha and beta constants to build the rayleigh
% damping matrix C.
% INPUTS:
%       - dtype: select the damping type between
%           "M"  --> C = alpha * M              (beta is set to 0)
%           "K"  --> C = beta  * K              (alpha is set to 0)
%           "MK" --> C = alpha * M + beta  * K
%       - Q: vector with the quality factors
%       - f: eigenfrequency vector (same order as Q)
%       - plotfigure (optional): draws the graph of xi and Q
%
% OUTPUT: alfa, beta
%
% Jacopo Marconi, PhD at Politecnico di Milano
% created: 13th May 2022

function [alpha, beta] = damping_rayleigh(dtype, Q, f0, plotfigure)

if nargin < 4
    plotfigure = 0;
end

% force inputs to column vectors
Q = reshape(Q, [], 1);
f0 = reshape(f0, [], 1);

csi = 1./(2*Q);

switch upper(dtype)
    case 'M'
        if length(Q)>1
            error(' Only one Q factor is allowed for M-damping type')
        end
        om0 = 2*pi*f0;
        alpha = 2 * om0 * csi;
        beta = 0;
    case 'K'
        if length(Q)>1
            error(' Only one Q factor is allowed for K-damping type')
        end
        om0 = 2*pi*f0;
        alpha = 0;
        beta = 2*csi/om0;
    case {'KM','MK'}
        if isscalar(Q)
            error(' Define at least 2 Q factors and frequencies for MK-damping type')
        end
        om0 = 2*pi*f0;
        AA = [1./(2*om0) om0/2];
        if size(AA,1)==size(AA,2)
            XX  = AA \ csi;
        else
            XX = (AA'*AA) \ (AA'*csi);
        end
        alpha = XX(1);
        beta = XX(2);
end

if plotfigure == 1
    figure
    f = linspace(min(f0)*0.95, max(f0)*1.05, 1000);
    w = 2*pi*f;

    xi = alpha ./ (2*w) + beta * w/2;
    QQ = 1./(2*xi);    

    yyaxis left
    plot(f,xi)
    hold on
    xline(f0,'k:','LineWidth',1)
    xlabel('frequency')
    ylabel('\xi [-]')

    yyaxis right
    plot(f,QQ)
    yline(Q,'r:','LineWidth',1)
    ylabel('Q [-]')

    grid on
    title('Rayleigh Damping')
    set(gca,'FontSize',14)
end