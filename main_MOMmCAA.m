% ----------------------------------------------------------------------- %
% Example of use of the funcion MOMmCAA.m, which performs a Multi-Objective %
% Majority-minority Cellular Automata Algorithm (MOMmCAA).               %
% ----------------------------------------------------------------------- %
%   Author:  Juan Carlos Seck Tuoh Mora                                      %
%   Date:    31/06/2024                                                   %
%   E-mail:  jseck@uaeh.edu.mx                              %
% ----------------------------------------------------------------------- %
%   References:                                                           %
%                    %
% ----------------------------------------------------------------------- %
clear all; clc;

% Multi-objective function
MultiObjFnc = 'CEC2020_MMF4';

switch MultiObjFnc
    case 'CEC2020_MMF1'         % CEC2020_MMF1
        f_{1} = @(x) abs((x(1)-2));
        f_{2} = @(x) 1.0 - sqrt(abs((x(1)-2))) + 2.0*(x(2)-sin(6*pi*abs((x(1)-2))+pi))^2;
        MultiObj.fun = @(x) [f1(x(:,1),x(:,2)); f2(x(:,1),x(:,2))];        
        MultiObj.nVar = 2;
        MultiObj.var_min =  [-1, -1];
        MultiObj.var_max = [3, 3];
    case 'CEC2020_MMF2'         % CEC2020_MMF2
        f_{1} = @(x) x(1);
        f_{2} = @(x) 1.0 - sqrt(x(1)) + 2*((4*(resta_1(x(2))-x(1)^0.5)^2)-2*cos(20*(resta_1(x(2))-x(1)^0.5)*pi/sqrt(2))+2);
        MultiObj.fun = @(x) [f1(x(:,1),x(:,2)); f2(x(:,1),x(:,2))];
        MultiObj.nVar = 2;
        MultiObj.var_min =  [0, 0];
        MultiObj.var_max = [1, 2];
    case 'CEC2020_MMF4'         % CEC2020_MMF4
        f_{1} = @(x) abs(x(1));
        f_{2} = @(x) 1.0 - (x(1))^2 + 2*(resta_1(x(2))-sin(pi*abs(x(1))))^2;
        MultiObj.fun = @(x) [f1(x(:,1),x(:,2)); f2(x(:,1),x(:,2))];
        MultiObj.nVar = 2;
        MultiObj.var_min =  [-1, 0];
        MultiObj.var_max = [1, 2];

end

% Create the problem using the parameters above.
% construct the multiobjective vector function
F = @(x) [];
for i = 1 : MultiObj.nVar
    F = @(x) [F(x);f_{i}(x)];
end
    
% construct the transpose of the function (for MOPSO and NSGA-II)
Ft = @(x) [];
for i = 1 : MultiObj.nVar
    Ft = @(x) [Ft(x),f_{i}(transpose(x))];
end

MultiObj.fun     = Ft;

% Parameters
params.ns                  = 11;    % Number of smart-cells
params.neigh             = 3;    % Number of neighbors 
params.maxgen         = 500; % Maximum number of generations
params.proporcion     = 5;  % Proportion of rule effect
params.n_dec_inf       = 2;     % Lower bound for round rule
params.n_dec_sup      = 5;     % Upper bound for round rule
params.nrep              = 50;   % Number of solutions in the repository
params.num_inter      = 5;   % Number of intervals to classify the position of each dim value in every solution

% MOMmCAA
REP = MOMmCAA(params,MultiObj);


%%%%Function to subtract one if the value is greater than 1, otherwise it returns the same value 
%%%otherwise returns the same value
function c = resta_1(x)
c=x;
if c>1
    c = c-1;
end
end


    