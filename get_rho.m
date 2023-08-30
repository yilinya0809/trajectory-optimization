function rho = get_rho(altitude)
    pressure = 101325 * (1 - 2.25569e-5 * altitude)^5.25616;
    temperature = 288.14 - 0.00649 * altitude;
    rho = pressure / (287 * temperature);
end