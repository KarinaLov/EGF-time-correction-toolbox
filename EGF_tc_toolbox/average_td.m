function timedelayF = average_td(settingsfile,pair)         
% Derives the final time delay by averaging the measured relative
% timedelays. This is an alternative to inversion when the strat time is
% unknown. The relative time delays should be measured from cross 
% correlation between assumed errousnes stations and correcetd stations. 
%
% Input:
%       settingsfile = text file where the input values are defined
%       pair = the station pairs that the relative timedelay is measured
%               over
%
% Output:
%       timedelayF = final averaged timedelays
%
%
% Sub-function: read_settings.m, read_daily.m, 
%               rdsac.m (https://se.mathworks.com/matlabcentral/fileexchange/46356-rdsac-and-mksac-read-and-write-sac-seismic-data-file)

%
% Written by Karina Loviknes 
% 

% Default values from settings file:
[network,stations,first_day,last_day,channels,location,num_stat_cc,Fq,...
    fileName,fileformat,pz_file,dateformat,deci,missingfiles,bpf,norm,wl,...
    swl,perco,filenameO,fileformatO,dateformatO,datesm,fit..
    ] = read_settings(settingsfile,'CORRECT');

validateattributes(stations,{'cell'},{'nonempty'});
nost=length(stations);
nch = length(channels);

% Extract the pairs to be used
lp = length(pair);

% Define the dates to run
fd = datetime(first_day);
ld = datetime(last_day);
datevector = fd:ld;
dates1 = [char(first_day) '-' char(last_day)];
if isempty(datesm)
    dates2 = dates1;
else
    dates2 = [char(datesm(1)) '-' char(datesm(2))];
    fd2 = datetime(datesm(1));
    ld2 = datetime(datesm(2));
    datevector = fd2:ld2;
end
num_days = length(datevector);
dd1 = [1:num_days];

for ch = 1:nch
    channel = channels(ch);
for jj = 1:nost
    % Loop over all the station pairs
    station = char(stations(jj)); 
    stationN = [char(stations(jj))  '-' channel]
    
    for ll = 1:lp
        % Extract measured time delays from file
        filename1 = ['TD_' pair{ll} '-Z_' dates2 '.mat']
        if exist(filename1,'file') % Check that the file exists
            file2 = load(filename1);
            timedelay1(ll,:) = -file2.timedelay.timedelay;
            timedelay0(ll,:) = -file2.timedelay.timedelay0;
            linear_td1(ll,:) = -file2.timedelay.linear_td; 
                
            is1(ll) = linear_td1(ll,1)
            fs1(ll) = linear_td1(ll,end)
            dd_fit = polyfit(dd1(1:end),linear_td1(ll,1:end),1);
            slope1(ll)=dd_fit(1)
            
        else
            error(['Cannot find a mat.file with an estimated greens function for stationpair ' pair{ll} '. Fileformat must be: TD_' pair '_' dates2 '.mat' ])
        end
    end
    if strcmp(fit,'linfit')
        timedelay = sum(linear_td1,1)/3;
        stdC = std(linear_td1,1);
    else    
        timedelay = sum(timedelay1,1)/3;
        stdC = std(timedelay1,1);
    end
    
    dd1 = [1:num_days];
    dd_fit = polyfit(dd1(1:end),timedelay(1:end),1)
    dd_fit_eval = polyval(dd_fit,dd1);
    slope=dd_fit(1)/Fq
    rsq = 1 - (sum((timedelay-dd_fit_eval).^2))/((num_days-1)*var(timedelay));
    std1 = std(timedelay)/Fq;
    
    timedelayF = struct('timedelay',timedelay,'timedelayC',timedelay,'pair',stationN,'number_of_days',num_days,'slope',slope,'Rsquared',rsq,'Standard_deviation',std1);
    save(['FTD_' stationN '_' dates2 '.mat'],'timedelayF')

end
end
end
