% Refractive Index: https://doi.org/10.1038/s41597-023-02898-2

% https://doi.org/10.1364/OE.25.025574
gold.alpha = 72.4; % [1/um]

% https://doi.org/10.1016/j.tsf.2004.02.028
ZrO2.alpha = 1.32e-5; % [1/um]

% https://apps.dtic.mil/sti/citations/ADA158623
graphite.alpha = 15.84; % [1/um]

save("optical_properties.mat", "gold", "ZrO2", "graphite")