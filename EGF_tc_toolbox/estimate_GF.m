function egfs = estimate_GF(settingsfile)
% Estimates the green's function from noise recorded on stationA and
% stationB between first_day and last_day
%
% Input:
%       settingsfile = text file where the input values are defined
%
% Output:
%       egfs = struct containing the estimdecimaated Green's function, lag time, 
%               the number of the days, and name of station pair
%
%
% Sub-function: read_settings.m, str2filename.m, gen_response.m, 
% read_daily.m, prepros.m and cross_conv.m
%
% Written by Karina L??viknes 
% 

% Default values from settings file:
[network, stations, first_day, last_day, channels, location,...
    num_stat_cc, Fq, filename, fileformat, pz_filename, dateformat,...
    deci, missingfiles, bpf, norm, wl, swl, perco] = read_settings(...
    settingsfile,'EGF');

validateattributes(stations,{'cell'},{'nonempty'});
nost = length(stations);
nch = length(channels);

fd = datetime(first_day);
ld = datetime(last_day);
datevector = [fd:ld];
num_days = length(datevector); % Number of days
num_corr = num_days*24/swl; % Number of correlations

tspd = Fq*60*60*24; % Total samples per day
lcc = 2*tspd/(24/wl)-1; % Length of cross correlation function

% Design the filter using the given cutoff frquencies and designfilt
df1 = designfilt('bandpassiir', 'FilterOrder', 4, ...
    'HalfPowerFrequency1', bpf(1), 'HalfPowerFrequency2', bpf(2), ...
    'SampleRate', Fq, 'DesignMethod', 'butter');

% Loop over all the channels:
sp = 0; % Count the station pairs
for ch = 1:nch
    channel = channels(ch);
    ii = 0;   

    % EXTRACT THE DAILY FILES FOR EACH STATION:
    % Preallocate for speed:
    resp = zeros(nost,tspd);
    Sdata = zeros(num_corr,tspd,nost);
    for jj = 1:nost
        % Loop over all the stations
        station = char(stations(jj))

        % Pole zero file for the stations: 
        pz_file = str2filename(pz_filename, station, dateformat,...
            'channels', channel, 'network', network)
        resp1 = gen_response(tspd,Fq,pz_file);
        resp(jj,:) = resp1;

        % Extract the data from the daily sac file (the sac-file needs to 
        % havethe format 'stationname-ch.yyyy-mm-dd.sac'
        Sdata1 = read_daily(network, station, channel, location,...
            datevector, filename, fileformat, dateformat, Fq, deci,...
            missingfiles); 
        Sdata(:,:,jj) = Sdata1; 
    end

    % LOOP OVER ALL THE STATIONPAIRS:
    for jj = 1:nost-1
        stationA = char(stations(jj))  
        respA = resp(jj,:);
        SAdata = Sdata(:,:,jj);

        for kk = 1:num_stat_cc-ii
            % check that we're not running out of stations on the list
            if jj+kk > nost
                continue
            end
            
            sp = sp+1;

            stationB = char(stations(jj+kk))
            respB = resp(jj+kk,:);
            SBdata = Sdata(:,:,jj+kk);

            pair = [stationA '-' stationB '-' channel]
            dates = [char(first_day) '-' char(last_day)];

            % Preallocate for speed:
            EGF = zeros(num_corr,lcc);

            % ESTIMATE THE GREEN'S FUNCTION:
            nk=24/swl;
            for d = 1:num_days

                % Preprocess the data for each station:
                SAprosd = prepros(SAdata(d,:),Fq,df1,respA,channel,norm);
                SBprosd = prepros(SBdata(d,:),Fq,df1,respB,channel,norm);

                % Cross correlations:
                [EGF1 lag] = cross_conv(SAprosd,SBprosd,Fq,wl,swl,perco);
                k=d*nk;

                EGF(k-(nk-1):k,:)=EGF1;

                % Save the daily cross correlations as SAC-files:
                % Header=struct('DELTA',1/Fq,'B',lag(1)/Fq,'E',...
                % lag(end)/Fq,'KSTNM', pair,'KHOLE',00,'KCMPNM',...
                % channels,'KNETWK',network,'NZDTTM',...
                % datevec(datevector(1)));
                % mksac(['Egf_' pair '_' datevector(d) '.SAC'],stack,...
                % datenum(first_day),Header)

            end
            stack = sum(EGF);

            estimatedGF = struct('EGF', EGF, 'lag', lag,...
                'number_of_days', num_days, 'pair', pair);
            egfs(sp) = estimatedGF;
            save(['Egf_' pair '_' dates '.mat'],'estimatedGF','-v7.3')

            Header=struct('DELTA', 1/Fq, 'B', lag(1)/Fq, 'E',...
                lag(end)/Fq, 'KSTNM', pair, 'KHOLE', 00,...
                'KCMPNM', channel, 'KNETWK', network,...
                'NZDTTM',datevec(datevector(1)));
            mksac(['Egf_' pair '_' dates '.SAC'], stack,...
                datenum(first_day),Header)

        end
        % Make sure the stations are cross correlted with the rigth number
        % of stations: 
        if ii >= num_stat_cc
            ii = 0; 
        else
            ii = ii + 1;
        end
    end
end
end




