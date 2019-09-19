function delay = measure_timeshift(settingsfile)         
% Calculate timing difference for a stationpair by comparing daily cross
% correlations to a reference trace stacked over a specified time period
%
% Input:
%       settingsfile = text file where the input values are defined
%
% Output:
%       delay = the measured time delays
%
%
% Sub-function: read_settings.m, make_reference.m and cross_conv.m
%
% Written by Karina Løviknes 
% 

% Default values from settings file
[network, stations, first_day, last_day, channels, location,...
    num_stat_cc, Fq, datesm, bpf, iterations, lag_red, stackperiod,...
    signal_part, thr] = read_settings(settingsfile, 'TD');

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
    datevector2 = datevector;
else
    dates2 = [char(datesm(1)) '-' char(datesm(2))];
    fd2 = datetime(datesm(1));
    ld2 = datetime(datesm(2));
    fdi = find(datevector==fd2);
    ldi = find(datevector==ld2);
    datevector2 = fd2:ld2;
end
num_days=length(datevector2);

% If specified, filter the daily cross correlations:
if ~isempty(bpf)
    % Design the filter using the given cutoff frquencies and designfilt
    ord = 4; % Filter order, default
    df = designfilt('bandpassiir','FilterOrder',ord, ...
        'HalfPowerFrequency1',bpf(1),'HalfPowerFrequency2',bpf(2), ...
        'SampleRate',Fq,'DesignMethod','butter');
end

sp = 0; % Count the station pairs
for ch = 1:nch
    ii = 0;
    channel = channels(ch);
    for jj = 1:nost
        % Loop over all the station pairs
        stationA = char(stations(jj));    

        for kk=1:num_stat_cc-ii
             sp=sp+1;

            stationB = char(stations(jj+kk));

            pair = [stationA '-' stationB '-' channel]

            % Extract information about the reference
            stackp1 = stackperiod{1};
            for k = 2:length(stackperiod)
                stackdays{k-1} = stackperiod{k};
            end   
            if strcmp(stackp1,'increasing')
                stackp = 'oneday';
            else
                stackp = stackp1;
            end

            filename = ['Egf_' pair '_' dates1 '.mat'];
            if exist(filename,'file') % Check that the file exists
                file = load(filename);
                EGF = file.estimatedGF.EGF(fdi:ldi,:);
                lag = file.estimatedGF.lag;
                num_corr = length(EGF(:,1));
            else
                error(['Cannot find a mat.file with an estimated ',...
                    'Greens function for stationpair ' pair,...
                    '. Fileformat must be: Egf_' pair '_' dates1 '.mat' ])
            end

            % Reduce computional effort by only using +-lag_red time lag
            zerolag = find(lag==0);
            egf = EGF(:,zerolag-lag_red:zerolag+lag_red);
            cl = lag(zerolag-lag_red:zerolag+lag_red);
            Lc = length(cl);  
            zero_lag = find(cl==0);

            dd1 = linspace(1,num_days,num_corr); % The days to be plotted

            % Determine signal and noise area for calculating SNR:
            narrp = zero_lag+400:zero_lag+500;
            narrn = zero_lag-400:zero_lag-500;
            sarr = zero_lag-350:zero_lag+350;

            dayshift = egf;

            % Preallocate for speed:   
            delay_dyn = zeros(1,num_corr);
            delay_dyn0 = zeros(1,num_corr);
            delay_pos = zeros(1,num_corr);
            delay_neg = zeros(1,num_corr);

            for it=1:iterations
                % Calculate refernce traces based on specification:
                [ref,numref]=make_reference(dayshift,stackp,stackdays);
                % Empty vector for specifying how the time delay is
                % measured
                type = []; 

                k = 0; % Count the days
                for j=1:numref

                    % If specified, filter the reference
                    if ~isempty(bpf)
                        % Filter the refrence
                        reff = filtfilt(df,ref(j,:));
                    else
                        reff = ref(j,:);
                    end

                    % Separate the time lags:
                    zero_lag = find(cl==0);
                    r_neg1 = flip(reff(1,1:zero_lag));
                    r_pos1 = reff(1,zero_lag:end);

                    % Normalize
                    r_neg = r_neg1/max(abs(r_neg1));
                    r_pos = r_pos1/max(abs(r_pos1));

                    % Find the static timing error (the unsymmetry of the
                    % reference):
                    delay_stat(it) = 0; %cl(find(abs(reff)==max(abs(reff))))
                    % Measure over the entire waveform
                    refwhn = reff/max(abs(reff)); 

                    % Loop over the daily cross correlations
                    for d=1:num_corr
                        k = k+1;

                        % Calculate the SNR:
                        snrm(d) = max(abs(egf(d,sarr)))/std(egf(d,narrp)); 

                        % If specified, filter the daily cross correlations:
                        if ~isempty(bpf)
                            egff(k,:) = filtfilt(df,egf(k,:));
                        else
                            egff(k,:) = egf(k,:);
                        end          

                        % Daily trace
                        s_neg1 = flip(egff(k,1:zero_lag));
                        s_pos1 = egff(k,zero_lag:end);

                        % Normalize
                        s_neg = s_neg1/max(abs(s_neg1));
                        s_pos = s_pos1/max(abs(s_pos1));
                        % Measure over the entire waveform
                        crconvfn = egff(k,:)/max(abs(egff(k,:)));  

                        % Calculate timing difference:
                        [td_neg1,lag_neg] = cross_conv(r_neg,s_neg,Fq);            
                        [td_pos1,lag_pos] = cross_conv(r_pos,s_pos,Fq);
                        [td_wh,lag_wh] = cross_conv(refwhn,crconvfn,Fq);

                        % Normalize the correlations:
                        au_nr = cross_conv(r_neg,r_neg,Fq);
                        au_pr = cross_conv(r_pos,r_pos,Fq);
                        au_ns = cross_conv(s_neg,s_neg,Fq);
                        au_ps = cross_conv(s_pos,s_pos,Fq);

                        au_crc = cross_conv(crconvfn,crconvfn,Fq);
                        au_ref = cross_conv(refwhn,refwhn,Fq);

                        td_neg = td_neg1/max(abs(sqrt(au_nr.*au_ns)));
                        td_pos = td_pos1/max(abs(sqrt(au_pr.*au_ps)));

                        td_whn = td_wh/max(abs(sqrt(au_ref.*au_crc)));

                        % Calculate the maximum amplitude of the normalized
                        % signal
                        maxtdn = max(td_neg);
                        maxtdp = max(td_pos);

                        maxtwh = max(td_whn);

                        % Only use signals with correlation coefficient 
                        % above the specified threrhold
                        if maxtdp>thr && maxtdn>thr && (strcmp(...
                                signal_part, 'separated') ||...
                                strcmp(signal_part, 'all'))
                            % Find the time shift of each side of the
                            % signal:
                            delay_neg(k) = lag_neg(find(td_neg==max(...
                                td_neg)));
                            delay_pos(k) = lag_pos(find(td_pos==max(...
                                td_pos)));

                            % Check that the delay is of the same size at
                            % both sides of the signal, but with opposite 
                            % direction:
                            if (-(delay_neg(k)+2))<=delay_pos(k) &&...
                                    delay_pos(k)<=(-(delay_neg(k)-2))
                                % Calculate the average delay:
                                delay_dyn0(k) = (delay_pos(k)-...
                                    delay_neg(k))/2;
                                delay_dyn(k) = 0;

                                type=[type 's'];
                            else
                                % If not the delay is set to zero or the 
                                %delay of the previous day:
                                delay_dyn0(k) = NaN;                   
                                if k>1
                                    % The delay can not measured, and is 
                                    % therefore set as the same as the 
                                    % previous day
                                    delay_dyn(k) = delay_dyn(k-1); 
                                else
                                    % If it is the first day the delay is 
                                    % set to zero
                                    delay_dyn(k) = NaN;
                                end                   
                                type = [type '0'];

                            end
                        elseif maxtwh>thr &&...
                                (strcmp(signal_part, 'whole') ||...
                                strcmp(signal_part, 'all'))
                            % Measure the delay over the entire waveform       
                            delay_dyn0(k) = lag_wh(find(td_whn==...
                                max(td_whn)));
                            delay_dyn(k) = delay_dyn0(k);

                            type = [type 'w'];

                        elseif maxtdp>0.4 && maxtdn<0.4 &&...
                                (strcmp(signal_part, 'positive') ||...
                                strcmp(signal_part, 'all'))
                            % The positive amplitudes are higher than the 
                            % negative, only use the positive side of the 
                            % signal:
                            delay_pos(k) = lag_pos(find(...
                                td_pos==max(td_pos)));
                            delay_dyn0(k) = delay_pos(k);    
                            delay_dyn(k) = delay_dyn0(k);

                            type = [type 'p'];

                        elseif maxtdn>0.4 && maxtdp<0.4 &&...
                                (strcmp(signal_part, 'negative') ||...
                                strcmp(signal_part, 'all'))
                            % The negative amplitudes are higher than the 
                            % positive, only use the negative side of the 
                            % signal:
                            delay_neg(k) = lag_neg(find(...
                                td_neg==max(td_neg)));
                            delay_dyn0(k) = -delay_neg(k);
                            delay_dyn(k) = delay_dyn0(k);

                            type = [type 'n'];

                        else
                            % No part of the signal have high enough 
                            % quality (correlation coefficient) to be used,
                            % the delay is set to zero
                            delay_dyn0(k) = NaN;

                            if k>1
                                % The delay can not measured, and is 
                                % therefore set as the same as the previous
                                % day
                                delay_dyn(k) = delay_dyn(k-1);  
                            else
                                % If it is the first day the delay is set
                                % to zero
                                delay_dyn(k) = 0;
                            end

                            type = [type '0'];
                        end
                        delay_dyn0(k) = delay_dyn0(k) + delay_stat(it);
                        delay_dyn(k) = delay_dyn(k) + delay_stat(it);
                    end
                end

                % Fit the line:
                dd_fit = polyfit(dd1,delay_dyn,1);
                dd_fit_eval = polyval(dd_fit,dd1);

                for kk = 1:k
                    % Correct for found timing errors:
                    t01 = -dd_fit_eval(kk)/Fq;
                    omgc = exp([0:Lc-1]*Fq/Lc*-1i*2*pi*t01);
                    shift1 = fft(egf(kk,:)).*omgc;
                    shift2 = ifft(shift1,'symmetric');
                    dayshift(kk,:) = shift2;
                end 

                % Change the refernce trace for next iteration:
                if strcmp(stackp1,'increasing') 
                    stackp = 'firstdays';
                    stackdaysi = str2num(stackdays{2})*it;
                    if num_corr>stackdaysi
                        stackdays{2} = num2str(stackdaysi);
                    end
                end

            timedelay = struct('timedelay', delay_dyn,...
                'timedelay0', delay_dyn0, 'linear_td', dd_fit_eval,...
                'reference', reff, 'SNR', snrm, 'pair', pair,...
                'number_of_days', num_days, 'type', type);
            delay(sp) = timedelay;
            save(['TD_' pair '_' dates2 '.mat'],'timedelay')
            end
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
