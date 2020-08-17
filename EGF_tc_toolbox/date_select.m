 function [tdd tddc dayn stations]= date_select(settingsfile,dato)
% Finds the estimated time delay at specified dates
%
% Input:
%       settingsfile = text file where the input values are defined 
%       dato = the date the specific drift should be given
%
% Output:
%       tdd = the time shift on the given date
%       tddc = the corrected time shift on the given date
%       dayn = the day number of the given date
%
%
%
% Sub-function: read_settings.m
%
% Written by Karina Loviknes 
% 

% Default values from settings file
[network,stations,first_day,last_day,channels,location,num_stat_cc,Fq,filename,fileformat,pz_file,dateformat,deci,missingfiles,bpf,norm,wl,swl,perco] = read_settings(settingsfile,'EGF');
[network,stations,first_day,last_day,channels,location,num_stat_cc,Fq,xaxis,yaxis,titl,bpfp,lag_red,datesm] = read_settings(settingsfile,'PLOT');

validateattributes(stations,{'cell'},{'nonempty'});
nost = length(stations);
nch = length(channels);

% Find the the spesified date:
fd = datetime(first_day);
ld = datetime(last_day);
dates1 = [char(first_day) '-' char(last_day)];
if isempty(datesm)
    dates2 = dates1; 
else
    dates2 = [char(datesm(1)) '-' char(datesm(2))];
    fd2 = datetime(datesm(1));
    ld2 = datetime(datesm(2));
end

for ch=1:nch
    channel = channels(ch);
for jj=1:nost
    % Loop over all the station pairs
    station = char(stations(jj));
    stationN = [char(stations(jj))  '-' channel];
    
    filename1=['FTD_' stationN '_' dates1 '.mat']; 
    filename2=['FTD_' stationN '_' dates2 '.mat']; 
    if exist(filename1,'file') % Check that the file exists
        filename1;
        file2 = load(filename1);
        timedelay = file2.timedelayF.timedelay; % Relative timedelay
        timedelayC = file2.timedelayF.timedelayC; % Correcetd timedelay
        datevector = [fd:ld];

    elseif exist(filename2,'file') % Check that the file exists
        filename2;
        file = load(filename2);
        timedelay = file.timedelayF.timedelay;
        timedelayC = file.timedelayF.timedelayC;
        datevector = [fd2:ld2];
       
    else
        error(['Cannot find a mat.file with an estimated greens function for stationpair ' stationN '. Fileformat must be: FTD_' stationN '_' dates1 '/' dates2 '.mat' ])
    end
    
    num_days = length(datevector); % Number of days
    num_corr = num_days*24/swl;
    dd1 = linspace(1,num_days,num_corr); % The days to be plotted

    dayd = datetime(dato);
    dayn = find(datevector==dayd);
    
    if isempty(dayn)
        error(['No time shift is measured for ' dayn])
    end
    
    tdd(jj) = timedelay(dayn)/Fq;
    tddc(jj) = timedelayC(dayn)/Fq;
    
    disp(['Timedelay for ' stationN ' on ' dato ' is ' num2str(tddc(jj)) ' s']);

end
hold off
end
end
