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
[network,stations,first_day,last_day,channels,location,num_stat_cc,Fq,bpf,iterations,lag_red,stackperiod] = read_settings(settingsfile,'TD');

validateattributes(stations,{'cell'},{'nonempty'});
nost = length(stations);

% Design the filter using the given cutoff frquencies and designfilt
ord = 4; % Filter order, default
df = designfilt('bandpassiir','FilterOrder',ord, ...
    'HalfPowerFrequency1',bpf(1),'HalfPowerFrequency2',bpf(2), ...
    'SampleRate',Fq,'DesignMethod','butter');

sp = 0; % Count the station pairs
for jj = 1:nost
    % Loop over all the station pairs
    stationA = char(stations(jj));    

for kk=1:num_stat_cc-(jj-1)
     sp=sp+1;
    
    stationB = char(stations(jj+kk));

    pair = [stationA '-' stationB]
    dates = [char(first_day) '-' char(last_day)];

    filename = ['Egf_' pair '_' dates '.mat'];
    if exist(filename,'file') % Check that the file exists
        file = load(filename);
        EGF = file.estimatedGF.EGF;
        lag = file.estimatedGF.lag;
        num_days = file.estimatedGF.number_of_days;
        num_corr = length(EGF(:,1));
    else
        error(['Cannot find a mat.file with an estimated Greens function for stationpair ' pair '. Fileformat must be: Egf_' pair '_' dates '.mat' ])
    end

    % Reduce computional effort by only using +-lag_red time lag
    zerolag = find(lag==0);
    egf = EGF(:,zerolag-lag_red:zerolag+lag_red);
    cl = lag(zerolag-lag_red:zerolag+lag_red);
    Lc = length(cl);
    
    dayshift = egf;

    % Preallocate for speed:   
    delay_dyn = zeros(1,num_corr);
    delay_dyn0 = zeros(1,num_corr);
    delay_pos = zeros(1,num_corr);
    delay_neg = zeros(1,num_corr);
    
    for it=1:iterations
        % Calculate refernce traces based on spesification:
        [ref,numref]=make_reference(dayshift,stackperiod);

        type = []; % Empty vector for spesifying how the time delay is measured

        k = 0; % Count the days
    for j=1:numref

        % Filter the refrence
        reff = filtfilt(df,ref(j,:));

        % Separate the time lags:
        zero_lag = find(cl==0);
        r_neg1 = flip(reff(1,1:zero_lag));
        r_pos1 = reff(1,zero_lag:end);

        % Normalize
        r_neg = r_neg1/max(abs(r_neg1));
        r_pos = r_pos1/max(abs(r_pos1));
        
        % Find the static timing error (the unsymmetry of the negative and positive side of the reference):
        [td_ref,lag_ref] = cross_conv(r_neg,r_pos,Fq);
        delay_stat = lag_ref(find(td_ref==max(td_ref)));
        
        refwhn = reff/max(abs(reff)); % Measure over the entire waveform
        
        % Loop over the daily cross correlations
        for d=1:num_corr
            k = k+1;

            % Filter the daily cross correlations
            egff(k,:) = filtfilt(df,egf(k,:));

            % Daily trace
            s_neg1 = flip(egff(k,1:zero_lag));
            s_pos1 = egff(k,zero_lag:end);

            % Normalize
            s_neg = s_neg1/max(abs(s_neg1));
            s_pos = s_pos1/max(abs(s_pos1));

            crconvfn = egff(k,:)/max(abs(egff(k,:)));  % Measure over the entire waveform

            % Calculate timing difference:
            [td_neg1,lag_neg] = cross_conv(r_neg,s_neg,Fq);            
            [td_pos1,lag_pos] = cross_conv(r_pos,s_pos,Fq);
            [td_wh,lag_wh] = cross_conv(crconvfn,refwhn,Fq);

            % Normalize the correlations:
            au_nr = cross_conv(r_neg,r_neg,Fq);
            au_pr = cross_conv(r_pos,r_pos,Fq);
            au_ns = cross_conv(s_neg,s_neg,Fq);
            au_ps = cross_conv(s_pos,s_pos,Fq);

            au_crc = cross_conv(crconvfn,crconvfn,Fq);
            au_ref = cross_conv(refwhn,refwhn,Fq);

            td_neg = td_neg1/max(abs(sqrt(au_nr.*au_ns)));
            td_pos = td_pos1/max(abs(sqrt(au_pr.*au_ps)));

            td_whn = td_wh/max(abs(sqrt(au_crc.*au_ref)));

            % Calculate the maximum amplitude of the normalized signal
            maxtdn = max(td_neg);
            maxtdp = max(td_pos);

            maxtwh = max(td_whn);

            % Only use signals with correlation coeffisient above the spesified thershold
            if maxtdp>0.4 && maxtdn>0.4
                % Find the time shift of each side of the signal:
                delay_neg(k) = lag_neg(find(td_neg==max(td_neg)));
                delay_pos(k) = lag_pos(find(td_pos==max(td_pos)));

                % Check that the delay is of the same size at both sides of the
                % signal, but with opposite direction:
                if (-(delay_neg(k)+2))<=delay_pos(k) && delay_pos(k)<=(-(delay_neg(k)-2))
                    % Calculate the average delay:
                    delay_dyn0(k) = (delay_pos(k)-delay_neg(k))/2;
                    delay_dyn(k) = 0;

                    type=[type 's'];
                else
                    % If not the delay is set to zero or the delay of the
                    % previous day:
                    delay_dyn0(k) = NaN;
                    
                    if k>1
                        % The delay can not measured, and is therefore set 
                        % as the same as the previous day
                        delay_dyn(k) = delay_dyn(k-1); 
                    else
                        % If it is the first day the delay is set to zero
                        delay_dyn(k) = NaN;
                    end
                    
                    type = [type '0'];
                end

            elseif maxtdp>0.4 && maxtdn<0.4
                % The positive amplitudes are higher than the negative, only 
                % use the positive side of the signal:
                delay_pos(k) = lag_pos(find(td_pos==max(td_pos)));
                delay_dyn0(k) = delay_pos(k);    
                delay_dyn(k) = delay_dyn0(k);

                type = [type 'p'];
                
            elseif maxtdn>0.4 && maxtdp<0.4
                % The negative amplitudes are higher than the positive, only 
                % use the negative side of the signal:
                delay_neg(k) = lag_neg(find(td_neg==max(td_neg)));
                delay_dyn0(k) = -delay_neg(k);
                delay_dyn(k) = delay_dyn0(k);

                type = [type 'n'];

            elseif maxtwh>0.4
                % Measure the delay over the entire waveform       
                delay_dyn0(k) = lag_wh(find(td_whn==max(td_whn)));
                delay_dyn(k) = delay_dyn0(k);
            
                type = [type 'w'];

            else
                % No part of the signal have high enough quality (correlation
                % coeffisient) to be used, the delay is set to zero
                delay_dyn0(k) = NaN;
                
                if k>1
                    % The delay can not measured, and is therefore set
                    % as the same as the previous day
                    delay_dyn(k) = delay_dyn(k-1);  
                else
                    % If it is the first day the delay is set to zero
                    delay_dyn(k) = 0;
                end
                
                type = [type '0'];
            end

            delay_dyn0(k) = delay_dyn0(k) + delay_stat;
            delay_dyn(k) = delay_dyn(k) + delay_stat;
            
            % Correct for found timing errors:
            if isnan(delay_dyn0(k))
                dayshift(k,:) = egff(k,:);
            else
                t01 = -(delay_dyn0(k)/Fq);
                omgc = exp([0:Lc-1]*Fq/Lc*-1i*2*pi*t01);
                shift1 = fft(egff(k,:)).*omgc;
                shift2 = ifft(shift1,'symmetric');
                dayshift(k,:) = shift2;
            end          
        end 
        
    end
    end
    
timedelay = struct('timedelay',delay_dyn,'timedelay0',delay_dyn0,'reference',reff,'pair',pair,'number_of_days',num_days,'type',type);
delay(sp) = timedelay;
save(['TD_' pair '_' dates '.mat'],'timedelay')

%Header=struct('DELTA',1/Fq,'B',lag(1)/Fq,'E',lag(end)/Fq,'KSTNM',pair,'KHOLE',00,'KCMPNM',channels,'KNETWK',network,'NZDTTM',datevec(datevector(1)));
%mksac(['Egf_' pair '_' dates '.SAC'],stack,datenum(first_day),Header)

end
end
end

