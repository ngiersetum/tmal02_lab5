function [parasitic, induced] = dragFunction(altitude, velocities)
% Calculate different drag contributions for one given altitude over one or
% several velocities
% 
% Inputs
%   aircraft (string)               - key in the aircraft database
%   altitude (scalar)               - flight altitude [m]
%   velocities (scalar or vector)   - flight velocities [m/s]
% 
% License
%   This program is part of an academic exercise for the course TMAL02,
%   Linköping University, year 2023. The program is therefore free for 
%   non-commercial academic use.
%
% Code History
%   https://github.com/ngiersetum/tmal02_lab3

    load('data/actable.mat')

    %% Technical data A340-300
    l = table2array(acdata("A340-300", "Fuse_Length")); % length in [m]
    l_cockpit = 6; % [m]
    l_empennage = 12; % [m]
    s_wing = 363.1; % wing reference are in [m^2]
    AR = table2array(acdata("A340-300", "AspectRatio"));
    r_fuse = table2array(acdata("A340-300", "Fuse_Width"))/2; % fuselage outside width in [m]
    slant_height_cockpit = sqrt(r_fuse^2 + l_cockpit^2);
    slant_height_empennage = sqrt(r_fuse^2 + l_empennage^2);

    s_wet_wing = 2.25 * s_wing;
    s_wet_fuse = 2*pi*(l - l_cockpit - l_empennage)*r_fuse + pi * r_fuse * slant_height_cockpit + pi * r_fuse * slant_height_empennage;
    s_wet_nacelle = ((pi*(table2array(acdata("A340-300", "NacelleWidth"))/2)^2) + 2 * (pi*table2array(acdata("A340-300", "NacelleWidth")) * table2array(acdata("A340-300", "NacelleLength"))));
    s_wet_tvert = 2.25 * table2array(acdata("A340-300", "TailVertArea"));
    s_wet_thor = 2.25 * table2array(acdata("A340-300", "TailHorArea"));
    
    %% Parameters for induced drag

    e = 0.776;   % Oswald Efficiency Number

    %% Drag Calculation

    [T, P, rho, a, mu] = ISAfunction(altitude);     % atmospheric conditions
    
    nu = mu/rho;

    for i=1:numel(velocities)

        vel = velocities(i);
        mach = vel/a;
        dynPress = 0.5 * rho * vel*vel;
        Re_pre = vel/nu;
        
        %% PARASITIC DRAG (Component build-up method)

        % characteristic lengths
        lchar_fuse = l;
        lchar_wing = table2array(acdata("A340-300", "MAC"));
        lchar_tvert = 5.8; % [m]
        lchar_thor = 2.9; % [m]
        lchar_nacelle = table2array(acdata("A340-300", "NacelleLength"));

        % chordwise location of maximum thickness point
        xc_wing = 0.4;
        xc_tvert = 0.25;
        xc_thor = 0.25;

        %thickness to chord
        tc_wing = 0.128;
        tc_tvert = 0.11;
        tc_thor = 0.11;

        % sweep angle at maximum thickness line
        phi_wing = deg2rad(30);
        phi_tvert = deg2rad(table2array(acdata("A340-300", "TailVertSweep")));
        phi_thor = deg2rad(table2array(acdata("A340-300", "TailHorSweep")));

        % slenderness ratios
        f_fuse = lchar_fuse / table2array(acdata("A340-300", "Fuse_Width"));
        f_nacelle = lchar_nacelle / table2array(acdata("A340-300", "NacelleWidth"));

        % component Reynolds numbers
        Re_fuse = Re_pre * lchar_fuse;
        Re_wing = Re_pre * lchar_wing;
        Re_tvert = Re_pre * lchar_tvert;
        Re_thor = Re_pre * lchar_thor;
        Re_nacelle = Re_pre * lchar_nacelle;
       
        % Skin Friction Drag
        Csd_fuse = 0.455 / ((log10(Re_fuse))^2.58 * (1 + 0.144 * mach^2)^0.65);
        Csd_wing = 0.455 / ((log10(Re_wing))^2.58 * (1 + 0.144 * mach^2)^0.65);
        Csd_tvert = 0.455 / ((log10(Re_tvert))^2.58 * (1 + 0.144 * mach^2)^0.65);
        Csd_thor = 0.455 / ((log10(Re_thor))^2.58 * (1 + 0.144 * mach^2)^0.65);
        Csd_nacelle = 0.455 / ((log10(Re_nacelle))^2.58 * (1 + 0.144 * mach^2)^0.65);

        % Form Drag
        FF_fuse = (1 + 5/(f_fuse^(1.5)) + f_fuse/400);
        FF_wing = (1 + (0.6/xc_wing) * tc_wing + 100 * tc_wing^4) * (1.34 * mach^(0.18) * (cos(phi_wing))^(0.28));
        FF_tvert = (1 + (0.6/xc_tvert) * tc_tvert + 100* tc_tvert^4) * (1.34 * mach^(0.18) * (cos(phi_tvert))^(0.28));
        FF_thor = (1 + (0.6/xc_thor) * tc_thor + 100* tc_thor^4) * (1.34 * mach^(0.18) * (cos(phi_thor))^(0.28));
        FF_nacelle = 1 + (0.35 / f_nacelle);
        

        % Interference Drag factors
        Q_fuse = 1.0;
        Q_wing = 1.0;
        Q_tvert = 1.06;
        Q_thor = 1.06;
        Q_nacelle = 1.3;


        % Component drag
        Cd_fuse = (Csd_fuse * FF_fuse * Q_fuse * s_wet_fuse)/s_wing;
        Cd_wing = (Csd_wing * FF_wing * Q_wing * s_wet_wing)/s_wing;
        Cd_tvert = (Csd_tvert * FF_tvert * Q_tvert * s_wet_tvert)/s_wing;
        Cd_thor = (Csd_thor * FF_thor * Q_thor * s_wet_thor)/s_wing;
        Cd_nacelle = 4 * (Csd_nacelle * FF_nacelle * Q_nacelle * s_wet_nacelle)/s_wing;
        

        % Overall component drag
        Cd = (Cd_fuse + Cd_wing + Cd_tvert + Cd_thor +  Cd_nacelle) * 1.1; % 10% misc drag

        parasitic(i) = Cd * s_wing * dynPress;
        %% LIFT-INDUCED DRAG

        weight = 0.8 * table2array(acdata("A340-300","MTOW")); % [kg], 80% of MTOW
        lift = weight * 9.80665; % lift is equal to weight in cruise (L = mg)

        cl = lift / (s_wing * dynPress);

        cdi = (cl*cl) / (e*pi*AR);

        induced(i) = dynPress * s_wing * cdi;
    end
end