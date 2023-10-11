function [T, P, rho, a, mu, liuID1, liuID2] = ISAfunction(Z)
% ISA function - Calculate atmospheric properties for one or several
% altitudes
%
% Synopsis
%   [T P rho a mu liuid1 liuid2] = ISAfunction(Z)
%
% Inputs
%   (scalar or vector) Z      Geometric altitude in [m]
%
% Outputs
%   (scalar or vector) T      Temperature, K
%   (scalar or vector) P      Pressure, Pa
%   (scalar or vector) rho    Density, kg/m^3
%   (scalar or vector) a      Sound speed, m/s
%   (scalar or vector) mu     Dynamic viscosity, kg/(m*s)
%
% References
%   Lab Reference Document (US Standard Atmosphere 1976)
%
% Authors
    liuID1 = "nikgi434"; % Niklas Gierse
    liuID2 = "leomu719"; % Leonhard Muehlstrasser
%
% License
%   This program is part of an academic exercise for the course TMAL02,
%   Linköping University, year 2023. The program is therefore free for 
%   non-commercial academic use.
%
% Code History
%   https://github.com/ngiersetum/tmal02_lab1

%% The executable code starts here:
% Define constants and reference values
    
    g0 = 9.80665; % [m/s^2]
    r0 = 6356766; % [m]
    R = 8.31432e3; % [(N*m)/(kmol*K)]
    gamma = 1.4; % [-] (heat capacity ratio)
    beta = 1.458e-6; % [kg/(s*m*K^1/2)]
    S = 110.4; % [K] (Sutherland Constant)
    Hb = [0, 11000, 20000, 32000, 47000, 51000, 71000, 84852]; %[m']
    bet = [-0.0065, 0, 0.001, 0.0028, 0, -0.0028, -0.0020]; %[K/m]
    Tb = [288.15, 216.650, 216.650, 228.650, 270.650, 270.560, 214.650, 186.8673]; % [K]
    pb = [101325, 22632, 5474.8, 868.01, 110.9, 66.938, 3.9564, 0.39814]; % [Pa]
    M0 = 28.9644; %[kg/kmol]
    
    
% Determine number of input altitudes
    nInput = length(Z);
    
% Convert input (geometric altitude) to geopotential altitude for
% calculations (eq 18 in document)
    H = (r0.*Z)./(r0 + Z);

% Initialize arrays just so they don't change size every iteration later
    T = zeros(size(Z));
    P = zeros(size(Z));
    rho = zeros(size(Z));
    a = zeros(size(Z));
    mu = zeros(size(Z));
        
% Main loop for all calculations
    for j = 1:nInput
        % Determine which layer we are in. Outside of boundaries all
        % properties will be set to NaN.
        b = 1;
        if H(j) <= Hb(1)        % below sea level
            b = NaN;
        elseif H(j) <= Hb(2)    % b=0
            b = 1;
        elseif H(j) <= Hb(3)    % b=1
            b = 2;
        elseif H(j) <= Hb(4)    % b=2
            b = 3;
        elseif H(j) <= Hb(5)    % b=3
            b = 4;
        elseif H(j) <= Hb(6)    % b=4
            b = 5;
        elseif H(j) <= Hb(7)    % b=5
            b = 6;
        elseif H(j) <= Hb(8)    % b=6
            b = 7;
        else                    % b=7
            b = NaN;
        end

        % Calculate properties. If outside boundaries result is NaN
        if ~isnan(b)
            % Calculate primary properties: temperature and pressure
            T(j) = Tb(b) + bet(b) .* (H(j) - Hb(b));
            if bet(b) ~= 0
                P(j) = pb(b) .* (Tb(b)./(Tb(b) + bet(b) .* (H(j) - Hb(b)))).^((g0 .* M0)/(R .* bet(b)));
            else
                P(j) = pb(b) .* exp((-g0 .* M0 .* (H(j) - Hb(b))) ./ (R .* Tb(b)));
            end

            % Calculate secondary properties: density, speed of sound, and dynamic viscosity
            rho(j) = (P(j).*M0)./(T(j).*R);
            a(j) = sqrt((gamma.*R.*T(j))./(M0));
            mu(j) = (beta.*(T(j).^(1.5)))./(T(j)+S);
        else
            T(j) = NaN;
            P(j) = NaN;
            rho(j) = NaN;
            a(j) = NaN;
            mu(j) = NaN;
        end
    end
end
