function h = apply_filterband(settingsfile,specification)
% Applies different bandpass filters to the noise correlation and plots the
% result 
%
% Input: 
%       settingfile = textfile where the input values are defined 
%       specification:
%               ref = the bandpass filters will be is applied to the 
%                       reference used to measure the time shift 
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
[network,stations,first_day,last_day,channels,location,num_stat_cc,Fq,xaxis,lag_red,fc1,fc2,datesm,titl] = read_settings(settingsfile,'FILTER');

validateattributes(stations,{'cell'},{'nonempty'});
nost = length(stations);
nch = length(channels);

fd = datetime(first_day);
ld = datetime(last_day);
datevector = [fd:ld];
dates1 = [char(first_day) '-' char(last_day)];

% Define the dates to run the td
if isempty(datesm)
    dates2 = dates1;
    fdi = 1;
    ldi = length(datevector);
else
    dates2 = [char(datesm(1)) '-' char(datesm(2))];
    fd2 = datetime(datesm(1));
    ld2 = datetime(datesm(2));
    fdi = find(datevector==fd2);
    ldi = find(datevector==ld2);
end

nfc = floor(length(fc1)); % Minimum number of filterbands
alph = 'abcdefghijklmnopqrstuvwxyz';
% Loop over all the station pairs:
sp = 0; % Count the station pairs
ii = 0;
for ch = 1:nch
    channel = channels(ch);
for jj=1:nost   
    stationA = char(stations(jj));    

for kk = 1:num_stat_cc-(jj-1)
    sp = sp+1; % Count the station pairs
    
    stationB = char(stations(jj+kk));

    pair = [stationA '-' stationB '-' channel];
    
    filename1 = ['Egf_' pair '_' dates1 '.mat'];
    if exist(filename1,'file') % Check that the file exists
        file1 = load(filename1);
        EGF = file1.estimatedGF.EGF(fdi:ldi,:);
        lag = file1.estimatedGF.lag;
        num_corr = length(EGF(:,1));
    else
        error(['Cannot find a mat.file with an estimated greens function for stationpair ' pair '. Fileformat must be: egf_' pair '_' datesm '.mat' ])
    end

    % Reduce computional effort by only using +-lag_red time lag
    zerolag = find(lag==0);
    egf = EGF(:,zerolag-lag_red:zerolag+lag_red);
    cl = lag(zerolag-lag_red:zerolag+lag_red);
    lcc = length(egf(1,:));

    dd1 = [1 num_corr]; % The days to be plotted
    cld = [cl(1)/Fq:1/Fq:cl(end)/Fq]; % I divide lag with Fq to plot the ccs in s
    
     if strcmp(specification,'ref')
        filename2=['TD_' pair '_' dates2 '.mat']; 
        if exist(filename2,'file') % Check that the file exists
            file2 = load(filename2);
            ref = file2.timedelay.reference;
        else
            warning(['Cannot find a mat.file with an estimated greens function for stationpair ' pair '. Fileformat must be: TD_' pair '_' datesm '.mat' ])
        end
        
        tit = [pair ' - unfiltered'];
        if strcmp(titl,'alphabet')
            tit = alph(1);
        end

        filt = zeros(nfc,length(ref));    
        h=figure;
        subplot(nfc,1,1)
        plot(cld,ref)
        xlim(xaxis)
        title(tit,'FontSize', 12)
        hold on 
        for ll=2:nfc

            % Design the filter defined by the input cutoff frquencies:
            df = designfilt('bandpassiir','FilterOrder',ord, ...
            'HalfPowerFrequency1',fc1(ll),'HalfPowerFrequency2',fc2(ll), ...
            'SampleRate',Fq,'DesignMethod','butter');

            % filtfilt filters the signal in both directions
            filt(ll,:) = filtfilt(df,ref);
                   
            tit = [pair ' - filter: ' num2str(fc1(ll)) '-' num2str(fc2(ll))  ' Hz'];
            if strcmp(titl,'alphabet')
                tit = alph(ll);
            end

            subplot(nfc,1,ll)
            plot(cld,filt(ll,:))
            xlim(xaxis)
            title(tit,'FontSize', 12)
            hold on 
        end
        hold off
        
     elseif strcmp(specification,'stack')

        ref=sum(egf);
        filt = zeros(nfc,length(ref));    
        h=figure;
        for ll=1:nfc

            % Design the filter defined by the input cutoff frquencies:
            df = designfilt('bandpassiir','FilterOrder',ord, ...
            'HalfPowerFrequency1',fc1(ll),'HalfPowerFrequency2',fc2(ll), ...
            'SampleRate',Fq,'DesignMethod','butter');

            % filtfilt filters the signal in both directions
            filt(ll,:) = filtfilt(df,ref);
            
            tit = [pair ' - filter: ' num2str(fc1(ll)) '-' num2str(fc2(ll))  ' Hz'];
            if strcmp(titl,'alphabet')
                tit = alph(ll);
            end

            subplot(nfc,1,ll)
            plot(cld,filt(ll,:))
            xlim(xaxis)
            title(tit,'FontSize', 12)
            hold on 
        end
        hold off
    
    elseif strcmp(specification,'daily')        
        % Preallocate for speed:
        egff = zeros(num_corr,lcc);
        egfn = zeros(num_corr,lcc);

        for d=1:num_corr
        % Filter and normalize the daily cross correlations:
            egfn(d,:) = egf(d,:)/max(abs(egf(d,:)));
        end
        
        tit = [pair ' - unfiltered'];
        if strcmp(titl,'alphabet')
            tit = alph(1);
        end
            
        h = figure;
        subplot(1,nfc,1)
        imagesc(cld,dd1,egfn)
        title(tit,'FontSize', 12)
        axis([xaxis 1 num_corr])
        xlabel('Time (s)','FontSize', 14), ylabel('Days','FontSize', 14)         
        level = 50; 
        n = ceil(level/2);
        cmap1 = [linspace(1,1,n); linspace(0,1,n); linspace(0,1,n)]';
        cmap2 = [linspace(1,0,n); linspace(1,0,n); linspace(1,1,n)]';
        cmap = [cmap1; cmap2(2:end,:)];
        colormap(cmap);
        hold on
        for ll=2:nfc
        % Design the filter using the given cutoff frquencies and designfilt
        dfpll = designfilt('bandpassiir','FilterOrder',4, ...
            'HalfPowerFrequency1',fc1(ll),'HalfPowerFrequency2',fc2(ll), ...
            'SampleRate',Fq,'DesignMethod','butter');
        
            tit = [pair ' - filter: ' num2str(fc1(ll)) '-' num2str(fc2(ll))  ' Hz'];
            if strcmp(titl,'alphabet')
                tit = alph(ll);
            end

            for d=1:num_corr
                % Filter and normalize the daily cross correlations:
                egff(d,:) = filtfilt(dfpll,egf(d,:));
                egfn(d,:) = egff(d,:)/max(abs(egff(d,:)));
            end
            
            subplot(1,nfc,ll)
            imagesc(cld,dd1,egfn)
            title(tit,'FontSize', 12)
            axis([xaxis 1 num_corr])
            xlabel('Time (s)','FontSize', 14), ylabel('Days','FontSize', 14)         
            level = 50; 
            n = ceil(level/2);
            cmap1 = [linspace(1,1,n); linspace(0,1,n); linspace(0,1,n)]';
            cmap2 = [linspace(1,0,n); linspace(1,0,n); linspace(1,1,n)]';
            cmap = [cmap1; cmap2(2:end,:)];
            colormap(cmap);
            hold on
        end
        hold off
    end
end
% Make sure the stations are cross correlted with the rigth number of stations: 
if ii >= num_stat_cc
    ii = 0; 
else
    ii = ii + 1;
end
end
end
end
