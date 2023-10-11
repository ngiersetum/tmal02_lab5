function thrust = calculateThrust(static, altitude, mach)
    rho0 = 1.2250;
    [ignore1, ignore2, rho] = ISAfunction(altitude);
    thrust = static * (rho / rho0) * (1 - 0.25*mach);
end