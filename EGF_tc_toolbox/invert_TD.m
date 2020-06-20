function [fdelay fdelayc] = invert_TD(settingsfile,varargin)
% The delays found for each station pair is inverted to find the final
% time delay of each station
%
% Input:
%       settingsfile = text file where the input values are defined
%       varargin = define wheter to plot the inverted shift
%
% Output:
%     fdelay = the final delay after inversion
%     fdelayc = the final delay after inversion and correction
%
%
% Sub-function: read_settings.m
%
% Written by Karina Løviknes
%

% Default values from settings file
[network,stations,first_day,last_day,channels,location,num_stat_cc,Fq,...
    fit,fitperiod,RCS,yaxis,datesm] = read_settings(settingsfile,'INVERT');

validateattributes(stations,{'cell'},{'nonempty'});
nost=length(stations);
nch = length(channels)

% Define the dates to run
dates1 = [char(first_day) '-' char(last_day)];
if isempty(datesm)
    dates2 = dates1;
else
    dates2 = [char(datesm(1)) '-' char(datesm(2))];
end

nsp = nost*num_stat_cc/2; % Number of station pair 

sp = 0; % Count the stationpairs
G = zeros(nsp,nost); % Matrix for inversion:

for ch=1:nch
    channel = channels(ch);
for jj = 1:nost
    % Loop over all the station pairs
    stationA = char(stations(jj));    

for kk = 1:num_stat_cc-(jj-1)
    sp = sp+1;
    
    % Fill in the inversion matrix with the station pair combinations:
    G(sp,jj) = 1; 
    G(sp,kk+jj) = -1;
    
    stationB = char(stations(jj+kk));

    pair = [stationA '-' stationB '-' channel]
    
    % Extract measured time delays from file
    filename = ['TD_' pair '_' dates2 '.mat'] 
    if exist(filename,'file') % Check that the file exists
        file = load(filename);
        timedelay(sp,:) = file.timedelay.timedelay;
        timedelay0(sp,:) = file.timedelay.timedelay0;
        linear_td(sp,:) = file.timedelay.linear_td;
        num_days = file.timedelay.number_of_days;
        num_corr = length(timedelay(sp,:));       
    else
        error(['Cannot find a mat.file with an estimated greens function for stationpair ' pair '. Fileformat must be: TD_' pair '_' dates2 '.mat' ])
    end
    dd1 = linspace(1,num_days,num_corr); % The days to be plotted
    
    if strcmp(fit,'linfit')
    % Define the linear fit:
        if isempty(fitperiod)
            dd_fit = polyfit(dd1,timedelay,1);
            dd_fit_eval = polyval(dd_fit,dd1);

        elseif length(fitperiod)==3 && trcmp(statf,'all')
            dd_fit = polyfit(dd1(dfp1:dfp2),timedelay(dfp1:dfp2),1);
            dd_fit_eval = polyval(dd_fit,dd1);

        elseif length(fitperiod) >= 3        
            llp = 0;
            for ll = 1:lfp/3
                statf = fitperiod{1+llp};
                dfp1 = str2num(fitperiod{2+llp});
                dfp2 = str2num(fitperiod{3+llp});

                % Fit the line:
                if strcmp(stationA,statf) || strcmp(stationB,statf)
                    dd_fit = polyfit(dd1(dfp1:dfp2),timedelay(dfp1:dfp2),1);
                    dd_fit_eval = polyval(dd_fit,dd1);
                elseif trcmp(statf,'all')
                    dd_fit = polyfit(dd1(dfp1:dfp2),timedelay(dfp1:dfp2),1);
                    dd_fit_eval = polyval(dd_fit,dd1);
                else
                    dd_fit = polyfit(dd1,timedelay,1);
                    dd_fit_eval = polyval(dd_fit,dd1);
                end
                llp =llp+3;
            end
        end
        timeshift = dd_fit_eval;
    else
        timeshift = timedelay;
    end
end
end

% Inversion matrix:
G

% Calculate and plot the final time delay using the pseudo inverse
fdelay = pinv(G)*timeshift;
    
% Correct to absolute time delay by removing the shift of the given reliable station:
if isempty(RCS)
    warning('A reliable station is not given, the relative time is not corrected for: fdelayc = fdelay')
    fdelayc = fdelay;
else
    realm = find(strcmp(stations,cellstr(RCS)));
    if isempty(realm)
        warning('A reliable station is not found, the relative time is not corrected for: fdelayc = fdelay')
        fdelayc = fdelay;
    else  
        for i = 1:nost
            fdelayc(i,:) = fdelay(i,:) - fdelay(realm,:);
        end
    end
end

if strcmp(varargin,'plot')
    % Plot the corrected inverted timedelays 
    figure
    for ff = 1:nost    
        char(stations(ff)) 
        dd_fit = polyfit(dd1,fdelayc(ff,:),1);
        dd_fit_eval = polyval(dd_fit,dd1);
        
        slope=dd_fit(1)/Fq;

        rsq = 1 - (sum((fdelayc(ff,:)-dd_fit_eval).^2))/((num_days-1)*var(fdelayc(ff,:)));
        std1 = std(fdelayc(ff,:))/Fq;
        
        subplot(nost,1,ff)
        plot(dd1,fdelayc/Fq)
        axis([0 num_corr yaxis])
        title(['Corrected timedelay for ' char(stations(ff)) '-' channel])
        xlabel('days')
        hold on 
        
        timedelayF = struct('timedelay',fdelay(ff,:),'timedelayC',...
            fdelayc(ff,:),'pair',[char(stations(ff)) '-' channel],...
            'number_of_days',num_days,'slope',slope,'Rsquared',rsq,...
            'Standard_deviation',std1);
        save(['FTD_' char(stations(ff)) '-' channel '_' dates2 '.mat'],'timedelayF')
    end
    hold off 
end
end
end
