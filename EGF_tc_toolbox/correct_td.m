%function corrected_stat = correct_td(settingsfile)         
% Correct for the measured time shift and creates new
% corrected files
%
% Input:
%       settingsfile = text file where the input values are defined
%
% Output:
%       corrected_stat = the corrected files
%
%
% Sub-function: read_settings.m, read_daily.m, 
%               rdsac.m (https://se.mathworks.com/matlabcentral/fileexchange/46356-rdsac-and-mksac-read-and-write-sac-seismic-data-file)

%
% Written by Karina Løviknes 
% 

% Default values from settings file:
[network,stations,first_day,last_day,channels,location,num_stat_cc,Fq,...
    fileName,fileformat,pz_file,dateformat,deci,missingfiles,bpf,norm,wl,...
    swl,perco,filenameO,fileformatO,dateformatO,datesm] = read_settings(settingsfile,'CORRECT');

validateattributes(stations,{'cell'},{'nonempty'});
nost=length(stations);
nch = length(channels);
nsp = nch*nost*num_stat_cc/2; % Number of station pair 

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
    datevector2 = fd2:ld2;
end

for ch = 1:nch
    channel = channels(ch);
for jj = 1:nost
    % Loop over all the station pairs
    station = char(stations(jj)); 
    stationN = [char(stations(jj))  '-' channel]
    
    % Extract measured time delays from file
    filename1 = ['FTD_' stationN '_' dates1 '.mat'];
    filename2 = ['FTD_' stationN '_' dates2 '.mat'];
    if exist(filename1,'file') % Check that the file exists
        filename1
        file = load(filename1);
        timedelay = file.timedelayF.timedelay;
        timedelayC = file.timedelayF.timedelayC;
        num_corr = length(timedelay);  
    
    elseif exist(filename2,'file') % Check that the file exists
        filename2
        file = load(filename2);
        timedelay = file.timedelayF.timedelay;
        timedelayC = file.timedelayF.timedelayC;
        num_corr = length(timedelay); 
        datevector = datevector2;
    else
        error(['Cannot find a mat.file with an estimated greens function for stationpair ' stationN '. Fileformat must be: FTD_' stationN '_' dates2 '.mat' ])
    end
    num_days = length(datevector); % Number of days
    
     % Extract the data from all the daily files of station B
    [daily dinfo] = read_daily(network,station,channel,location,datevector,fileName,fileformat,dateformat,Fq,deci,missingfiles);
    
    starttime = dinfo.Starttime;
    starttime = datestr(starttime);
    endtime = dinfo.Endtime;
    count = dinfo.Count;
    Lc = length(daily(1,:));
    
    corrected = zeros(num_days,Lc);
    
        for ii = 1:num_days
            % Make sure the day exist:
            st1 = dinfo(ii).Starttime;
            if st1 ~= 0
                nq = 1;
                tdfs = abs(round(timedelayC(ii)*nq));
                tv = (1:count)';

                 zmsp = zeros(1,tdfs);
                 data11 = daily(ii,:);

                if timedelay(ii) > 0
                    data_2 = [zmsp data11(1,1:end-tdfs)];
                else
                    data_2 = [data11(1,1+tdfs:end) zmsp];
                end
                shift2 = data_2;
                corrected(ii,:) = shift2;

                % Extract the starttime and count of the day to be corrected
                starttime = datestr(dinfo(ii).Starttime);
                count = dinfo(ii).Count;

                if strcmp(fileformatO,'sac')
                    % Create a corrected sac file
                    Header=struct('DELTA',1/Fq,'B',0,'KSTNM',station,'KHOLE',00,'KCMPNM',channel(ch),'KNETWK',network,'NZDTTM',datevec(starttime));
                    sacname = str2filename(filenameO,station,dateformatO,'network',network,'channels',channels,'location',location,'datevector',datetime(starttime))
                    mksac(sacname,shift2,datenum(starttime),Header);

                elseif strcmp(fileformatO,'miniseed')
                    % Create a corrected miniseed file

                    mseedname = [network '.' station '.00.HH' channel(ch) '.D'];
                    mkmseed(mseedname,shift2,datenum(starttime),Fq);

                    % Change the filename of the outputname spesified in the settingsfile: 
                    mseedname1 = [network '.' station '.00.HH' channel(ch) '.' char(datetime(datestr(starttime),'Format','yyyy')) '.' char(datetime(starttime,'Format','DDD'))];
                    mseedname2 = str2filename(filenameO,station,dateformatO,'network',network,'channels',channels,'location',location,'datevector',datetime(starttime))
                    movefile(mseedname1,mseedname2);
                end
            end  
        end
    corrected_stat(jj) = struct('Corrected',corrected,'Name',stationN,'Dates',dates2);   
end
end
%end

