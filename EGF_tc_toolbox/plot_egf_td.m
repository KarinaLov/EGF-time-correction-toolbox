 function h = plot_egf_td(settingsfile,varargin)
% Plot the estimated Green's function and/or the measured time delays
%
% Input:
%       settingsfile = text file where the input values are defined 
%       varargin: specification about what to plot
%
% Output:
%      h = the figure with the plots
%
%
% Sub-function: read_settings.m
%
% Written by Karina Løviknes 
% 

% Default values from settings file
[network,stations,first_day,last_day,channels,location,num_stat_cc,Fq,xaxis,yaxis,bpfp,lag_red] = read_settings(settingsfile,'PLOT');

validateattributes(stations,{'cell'},{'nonempty'});
nost = length(stations);

nsp = nost*num_stat_cc/2; % Number of station pair 

if strcmp(varargin,'all')
    h=figure;
end

sp=0; % Count the station pairs
for jj=1:nost
    % Loop over all the station pairs
    stationA=char(stations(jj));    

for kk=1:num_stat_cc-(jj-1)
     sp=sp+1;
    
    stationB=char(stations(jj+kk));

    pair=[stationA '-' stationB];
    dates=[char(first_day) '-' char(last_day)];
    
    filename1=['Egf_' pair '_' dates '.mat'];
    if exist(filename1,'file') % Check that the file exists
        file1=load(filename1);
        EGF=file1.estimatedGF.EGF;
        lag=file1.estimatedGF.lag;
        num_days = file1.estimatedGF.number_of_days;
        num_corr = length(EGF(:,1));
    else
        error(['Cannot find a mat.file with an estimated greens function for stationpair ' pair '. Fileformat must be: Egf_' pair '_' dates '.mat' ])
    end

    filename2=['TD_' pair '_' dates '.mat']; 
    if exist(filename2,'file') % Check that the file exists
        file2 = load(filename2);
        timedelay = file2.timedelay.timedelay;
        timedelay0 = file2.timedelay.timedelay0;
        ref = file2.timedelay.reference;
    else
        error(['Cannot find a mat.file with an estimated greens function for stationpair ' pair '. Fileformat must be: TD_' pair '_' dates '.mat' ])
    end
    
    % Reduce computional effort by only using +-lag_red time lag
    zerolag = find(lag==0);
    egf = EGF(:,zerolag-lag_red:zerolag+lag_red);
    cl = lag(zerolag-lag_red:zerolag+lag_red);
    lcc = length(cl);
    
    % Normalize and filter the daily cross correlations:
    egfn = zeros(num_corr,lcc); % Preallocate for speed:
    for d = 1:num_corr
        % Normalize the daily cross correlations:
        egfn(d,:) = egf(d,:)/max(abs(egf(d,:)));
        
        % If spesified, filter the daily cross correlations:
        if ~isempty(bpfp)
            % Design the filter using the given cutoff frquencies and designfilt
            dfp = designfilt('bandpassiir','FilterOrder',4, ...
                'HalfPowerFrequency1',bpfp(1),'HalfPowerFrequency2',bpfp(2), ...
                'SampleRate',Fq,'DesignMethod','butter');
            egfn(d,: ) = filtfilt(dfp,egf(d,:));
        end
    end
    
    dd1 = linspace(1,num_days,num_corr); % The days to be plotted
    cld = [cl(1)/Fq:1/Fq:cl(end)/Fq]; % I divide lag with Fq to plot the ccs in s
        
    if isempty(varargin)
        % Default: Plot the daily Green's function and the timedelays  
        h=figure;
        subplot(1,2,1)
        imagesc(cld,dd1,egfn)
        title(['Cross correlations of ' pair ' between ' dates ],'FontSize', 15)
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
        
        subplot(2,2,2)
        plot(cld,ref)
        title(['Corrected reference trace used to measure timedelays'],'FontSize', 15)
        xlim([xaxis])
        subplot(2,2,4)
        plot(dd1,timedelay/Fq,dd1,timedelay0/Fq)
        legend('Continous timedelay','Measured timedelay')
        axis([1 num_days yaxis])
        title(['Relative timedelay for ' pair ' between ' dates],'FontSize', 15)
        xlabel('Days','FontSize', 16), ylabel('Time delay (s)','FontSize', 16)

     elseif strcmp(varargin,'EGF')
        % Plot ONLY the daily Green's function 
        h=figure;
        imagesc(cld,dd1,egfn)
        title(['Daily EGFs for ' pair ' between ' dates ],'FontSize', 14)
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

     elseif strcmp(varargin,'TD')
        % Plot ONLY the measured timedelays 
        h=figure;
        subplot(2,1,1)
        plot(cld,ref)
        title(['Corrected reference trace used to measure timedelays'],'FontSize', 15)
        xlim([-100 100])
        subplot(2,1,2)
        plot(dd1,timedelay/Fq,dd1,timedelay0/Fq)
        legend('Continous timedelay','Measured timedelay')
        axis([1 num_days yaxis])
        title(['Relative timedelay for ' pair ' between ' dates],'FontSize', 15)
        xlabel('Days','FontSize', 16), ylabel('Time delay (s)','FontSize', 16)

    elseif strcmp(varargin,'Daily')
        % Plot the daily Green's function as signals not amplitude 
        pp=1;
        h=figure;
        for j=1:num_corr
            egfp = egfn(j,:) + pp;
            
            plot(cld,egfp,'k')
            set(gca,'Ydir','reverse')
            title(['Daily cross correlation of ' pair ' between ' dates],'FontSize', 15)
            axis([xaxis -0.1 num_corr])
            xlabel('Time (s)','FontSize', 16), ylabel('Days','FontSize', 16)
            hold on
            
            pp = pp + 1.1;
        end
        hold off            
    elseif strcmp(varargin,'all')
        % Plot all the station pairs
        subplot(2,nsp/2,sp)
        imagesc(cld,dd1,egfn)
        title(['Daily EGFs for ' pair ],'FontSize', 14)
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
        hold on
    end
end
end
hold off
end
