function h = plot_spectrum(settingsfile)
% Plot the estimated Green's function and amplitude spectrum
%
% Input:
%       settingsfile = textfile containin where the input values are defined
%
% Output:
%      h=the figure with the plots
%

% Default values from settings file
[network,stations,first_day,last_day,channels,location,filename,fileformat,num_stat_cc,Fq,xaxis,yaxis,bpfp,lag_red] = read_settings(settingsfile,'PLOT');

validateattributes(stations,{'cell'},{'nonempty'});
nost = length(stations);

filtt = '';
sp = 0; % Count the station pairs
for jj = 1:nost
    % Loop over all the station pairs
    stationA = char(stations(jj));    

for kk = 1:num_stat_cc-(jj-1)
     sp = sp+1;
    
    stationB = char(stations(jj+kk));

    pair = [stationA '-' stationB];
    dates = [char(first_day) '-' char(last_day)];
    
    filename1 = ['Egf_' pair '_' dates '.mat'];
    if exist(filename1,'file') % Check that the file exists
        file1 = load(filename1);
        EGF = file1.estimatedGF.EGF;
        lag = file1.estimatedGF.lag;
        num_days = file1.estimatedGF.number_of_days;
    else
        error(['Cannot find a mat.file with an estimated greens function for stationpair ' pair '. Fileformat must be: egf_' pair '_' dates '.mat' ])
    end

    % Reduce computional effort by only using +-lag_red time lag
    zerolag = find(lag==0);
    egf = EGF(:,zerolag-lag_red:zerolag+lag_red);
    cl = lag(zerolag-lag_red:zerolag+lag_red);
    lcc = length(cl);
        
    egfn = zeros(num_corr,lcc); % Preallocate for speed:
    for d = 1:num_days
        % Normalize the daily cross correlations:
        egfn(d,:) = egf(d,:)/max(abs(egf(d,:)));
    end
    
    stack = sum(EGF);
    
    % Amplitude spectrum:
    freq = linspace(0,Fq,length(stack));
    aft = abs(fft(stack));
     
    dd1 = [1 num_days]; % The days to be plotted
    cld = [cl(1)/Fq:1/Fq:cl(end)/Fq]; % I divide lag with 10 to plot the ccs in s

    h = figure;
    subplot(1,2,1)
    imagesc(cld,dd1,egfn)
    title(['Daily cross correlation of ' pair ' between ' dates filtt],'FontSize', 15)
    axis([xaxis 1 num_days])
    xlabel('Time (s)','FontSize', 16), ylabel('Days','FontSize', 16)
    % Set the colorbar to red-white-blue:
    colorbar
    level = 50; 
    n = ceil(level/2);
    cmap1 = [linspace(1,1,n); linspace(0,1,n); linspace(0,1,n)]';
    cmap2 = [linspace(1,0,n); linspace(1,0,n); linspace(1,1,n)]';
    cmap = [cmap1; cmap2(2:end,:)];
    colormap(cmap);
        
    xlabel('Time (s)','FontSize', 16), ylabel('Days','FontSize', 16)
    subplot(2,2,2)
    plot(cld,stack(zerolag-lag_red:zerolag+lag_red))
    title(['Corrected reference trace used to measure timedelays'],'FontSize', 15)
    xlim([xaxis])
    subplot(2,2,4)
    plot(freq,aft)
    xlim([0 Fq/2])
    title(['Amplitude spectrum for ' pair ' between ' dates],'FontSize', 15)
    xlabel('Days','FontSize', 16), ylabel('Time delay (s)','FontSize', 16)
end
end
end
