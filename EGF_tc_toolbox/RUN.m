
settingsfile='settingsfile.txt';

% ESTIMATE THE GREEN'S FUNCTION:
%egfs = estimate_GF(settingsfile);

% MEASURE THE TIME SHIFT
delays = measure_timeshift(settingsfile);

% PLOT DAILY CROSS CORRELATION AND TIME SHIFTS
h = plot_egf_td(settingsfile);
h = plot_egf_td(settingsfile,'all');
%h = plot_egf_td(settingsfile,'Daily');
%h = plot_egf_td(settingsfile,'Frequency');

%% INVERT TO FIND THE TIME DELY OF EACH STATION
%[fdelay fdelayc] = invert_TD(settingsfile,'plot');

% Plot the results with distance, filters and ampliude spectrum
h = plot_distance(settingsfile);
h = apply_filterband(settingsfile,'stack');
h = apply_filterband(settingsfile,'daily');
%%
% Give the drift on a specified date:
dato = {'2015-10-01'}
[tdd tddc dayn]= date_select('settingsfile.txt',char(dato));

% Corrcet for measured time shift:
corrected_stat = correct_td(settingsfile);


