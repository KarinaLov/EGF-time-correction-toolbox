function h = apply_filterband(settingsfile,spesification)
% Applies different bandpass filters to the noise correlation and plots the
% result 
%
% Input: 
%       settingfile = textfile where the input values are defined 
%       spesification:
%               stack = the bandpass filters will be is applied to the 
%                       stacked noise correlation 
%               daily = the bandpass filters will be is applied to the 
%                       daily noise correlation 
%
% Output:
%       h = figure
%
%
% Sub-function: read_settings.m
%
% Written by Karina Løviknes 
%

ord = 4; % Filter order, default

% Default values from settings file
[network,stations,first_day,last_day,channels,location,num_stat_cc,Fq,xaxis,lag_red,fc1,fc2] = read_settings(settingsfile,'FILTER');

validateattributes(stations,{'cell'},{'nonempty'});
nost = length(stations);

nfc = floor(length(fc1)); % Minimum number of filterbands

% Loop over all the station pairs:
sp = 0; % Count the station pairs
for jj=1:nost   
    stationA = char(stations(jj));    

for kk = 1:num_stat_cc-(jj-1)
    sp = sp+1; % Count the station pairs
    
    stationB = char(stations(jj+kk));

    pair = [stationA '-' stationB];
    dates = [char(first_day) '-' char(last_day)];
    
    filename1 = ['Egf_' pair '_' dates '.mat'];
    if exist(filename1,'file') % Check that the file exists
        file1 = load(filename1);
        EGF = file1.estimatedGF.EGF;
        lag = file1.estimatedGF.lag;
        num_days = file1.estimatedGF.number_of_days;
        num_corr = length(EGF(:,1));
    else
        error(['Cannot find a mat.file with an estimated greens function for stationpair ' pair '. Fileformat must be: egf_' pair '_' dates '.mat' ])
    end

    % Reduce computional effort by only using +-lag_red time lag
    zerolag = find(lag==0);
    egf = EGF(:,zerolag-lag_red:zerolag+lag_red);
    cl = lag(zerolag-lag_red:zerolag+lag_red);
    lcc = length(egf(1,:));

    dd1 = [1 num_corr]; % The days to be plotted
    cld = [cl(1)/Fq:1/Fq:cl(end)/Fq]; % I divide lag with Fq to plot the ccs in s
    
    if strcmp(spesification,'stack')
        stack = sum(egf);   
        filt = zeros(nfc,length(stack));    
        h=figure;
        for ll=1:nfc

            % Design the filter defined by the input cutoff frquencies:
            df = designfilt('bandpassiir','FilterOrder',ord, ...
            'HalfPowerFrequency1',fc1(ll),'HalfPowerFrequency2',fc2(ll), ...
            'SampleRate',Fq,'DesignMethod','butter');

            % filtfilt filters the signal in both directions
            filt(ll,:) = filtfilt(df,stack);

            subplot(nfc,1,ll)
            plot(cld,filt(ll,:))
            xlim(xaxis)
            title([pair ': ' num2str(fc1(ll)) '-' num2str(fc2(ll))])
            hold on 
        end
        hold off
    
    elseif strcmp(spesification,'daily')        
        % Preallocate for speed:
        egff = zeros(num_corr,lcc);
        egfn = zeros(num_corr,lcc);

        h = figure;
        for ll=1:nfc
        % Design the filter using the given cutoff frquencies and designfilt
        dfpll = designfilt('bandpassiir','FilterOrder',4, ...
            'HalfPowerFrequency1',fc1(ll),'HalfPowerFrequency2',fc2(ll), ...
            'SampleRate',Fq,'DesignMethod','butter');

            for d=1:num_corr
                % Filter and normalize the daily cross correlations:
                egff(d,:) = filtfilt(dfpll,egf(d,:));
                egfn(d,:) = egff(d,:)/max(abs(egff(d,:)));
            end
                subplot(2,nfc/2,ll)
                imagesc(cld,dd1,egfn)
                title([pair ' - filtered ' num2str(fc1(ll)) '-' num2str(fc2(ll))  ' Hz'],'FontSize', 12)
                axis([xaxis 1 num_corr])
                colorbar
                xlabel('Time (s)','FontSize', 14), ylabel('Days','FontSize', 14)
                hold on
        end
        hold off
    end
end
end
end
