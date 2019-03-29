
settingsfile='settingsfile.txt';

% ESTIMATE THE GREEN'S FUNCTION:
egfs = estimate_GF(settingsfile);

% MEASURE THE TIME SHIFT
delays = measure_timeshift(settingsfile);

% PLOT DAILY CROSS CORRELATION AND TIME SHIFTS
h = plot_egf_td(settingsfile);
h = plot_egf_td(settingsfile,'all');


% INVERT TO FIND THE TIME DELY OF EACH STATION
[fdelay fdelayc] = invert_TD(settingsfile,'plot');

% Plot the results with distance, filters and ampliude spectrum
h = plot_distance(settingsfile)
h = apply_filterband(settingsfile,'stack')









