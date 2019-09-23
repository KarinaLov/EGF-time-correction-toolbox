 function h = plot_egf_td(settingsfile,varargin)
% Plot the estimated Green's function and/or the measured time delays
%
% Input:
%       settingsfile = text file where the input values are defined 
%       varargin: 
%               specification about what to plot
%               specification about wheter to save the figure
%
% Output:
%      h = the figure with the plots
%
%
% Sub-function: read_settings.m
%
% Written by Karina L??viknes 
% 

% Default values from settings file
[network, stations, first_day, last_day, channels, location,...
    num_stat_cc, Fq, xaxis, yaxis, titl, bpf, lag_red, datesm] =...
    read_settings(settingsfile, 'PLOT');

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
    datevector2=datevector;
else
    dates2 = [char(datesm(1)) '-' char(datesm(2))];
    fd2 = datetime(datesm(1));
    ld2 = datetime(datesm(2));
    fdi = find(datevector==fd2);
    ldi = find(datevector==ld2);
    datevector2 = [fd2:ld2];
end
num_days=length(datevector2); % Number of days
nsp = nch*nost*num_stat_cc/2; % Number of station pair 

h1=1;
if ~isempty(varargin) && strncmp(varargin{1},'all',3) 
    h=figure;
    h1=0;
end

alph = 'abcdefghijklmnopqrstuvwxyz';

% If spesified, filter the daily cross correlations:
filttxt = '';
if ~isempty(bpf)
    % Design the filter using the given cutoff frquencies and designfilt
    dfp = designfilt('bandpassiir', 'FilterOrder', 4,...
        'HalfPowerFrequency1', bpf(1), 'HalfPowerFrequency2', bpf(2),...
        'SampleRate', Fq, 'DesignMethod', 'butter');
    filttxt = [' - filt: ' num2str(bpf(1)) '-' num2str(bpf(2))];
end

sp = 0; % Count the station pairs
ii = 0;
for jj=1:nost-1
    % Loop over all the station pairs
    stationA=char(stations(jj));    

    for kk=1:num_stat_cc-ii
        % check that we're not running out of stations on the list
        if jj+kk > nost
            continue
        end
            
        sp=sp+1;
        spp=0;
        %h=figure
        for ch = 1:nch
            channel = channels(ch);
            spp=spp+1;

            stationB = char(stations(jj+kk));

            pair = [stationA '-' stationB '-' channel];

            filename1=['Egf_' pair '_' dates1 '.mat'];
            if java.io.File(filename).exists  % Check that the file exists                
                file1=load(filename1);
                EGF=file1.estimatedGF.EGF(fdi:ldi,:);
                lag=file1.estimatedGF.lag;
                %num_days = file1.estimatedGF.number_of_days;
                num_corr = length(EGF(:,1));
            else
                error(['Cannot find a mat.file with an estimated ',...
                    'greens function for stationpair ' pair,...
                    '. Fileformat must be: Egf_' pair '_' dates1 '.mat' ])
            end

            filename2=['TD_' pair '_' dates2 '.mat'];
            if java.io.File(filename2).exists  % Check that the file exists                
                file2 = load(filename2);
                timedelay = file2.timedelay.timedelay;
                timedelay0 = file2.timedelay.timedelay0;
                linear_td = file2.timedelay.linear_td;
                ref = file2.timedelay.reference;
            else
                warning(['Cannot find a mat.file with a measured time ',...
                    'delay for stationpair ' pair '. Fileformat must ',...
                    'be: TD_' pair '_' dates2 '.mat' ])
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
                if ~isempty(bpf)
                    egfn(d,: ) = filtfilt(dfp,egf(d,:));
                end
            end
            
            % The days to be plotted
            dd1 = linspace(1,num_days,num_corr); 
            % I divide lag with Fq to plot the ccs in s
            cld = [cl(1)/Fq:1/Fq:cl(end)/Fq]; 
            % I divide lag with Fq to plot the ccs in s
            lagd = [lag(1)/Fq:1/Fq:lag(end)/Fq]; 

            if isempty(varargin)
                % Default: Plot the daily Green's function and the
                % timedelays  
                h=figure;
                subplot(1,2,1)
                imagesc(cld,dd1,egfn)
                axis([xaxis 1 num_days])
                xlabel('Time (s)','FontSize', 16), ylabel('Days',...
                    'FontSize', 16)
                % Set the colorbar to red-white-blue:
                colorbar
                level = 50; 
                n = ceil(level/2);
                cmap1 = [linspace(1,1,n); linspace(0,1,n);...
                    linspace(0,1,n)]';
                cmap2 = [linspace(1,0,n); linspace(1,0,n);...
                    linspace(1,1,n)]';
                cmap = [cmap1; cmap2(2:end,:)];
                colormap(cmap);
                if strcmp(titl,'alphabet')
                    title(alph(1),'FontSize', 15)
                else
                    title(['Cross correlations of ' pair ' between ',...
                        dates2 ],'FontSize', 15)
                end

                subplot(2,2,2)
                plot(cld,ref)
                xlim([xaxis])
                if strcmp(titl,'alphabet')
                    title(alph(2),'FontSize', 15)
                else
                    title(['Corrected reference trace used to measure ',...
                        'timedelays'],'FontSize', 15)
                end

                subplot(2,2,4)
                plot(dd1,linear_td/Fq,dd1,timedelay0/Fq,'.')
                legend('Linearly fitted timedelay','Measured timedelay')
                axis([1 num_days yaxis])
                xlabel('Days','FontSize', 16), ylabel('Time delay (s)',...
                    'FontSize', 16)
                if strcmp(titl,'alphabet')
                    title(alph(3),'FontSize', 15)
                else
                    title(['Relative timedelay for ' pair ' between ',...
                        dates2],'FontSize', 15)
                end

            elseif strcmp(varargin{1},'save')
                % Default: Plot the daily Green's function and the
                % timedelays  
                h=figure;
                subplot(1,2,1)
                imagesc(cld,dd1,egfn)
                axis([xaxis 1 num_days])
                xlabel('Time (s)','FontSize', 16), ylabel('Days',...
                    'FontSize', 16)
                % Set the colorbar to red-white-blue:
                colorbar
                level = 50; 
                n = ceil(level/2);
                cmap1 = [linspace(1,1,n); linspace(0,1,n);...
                    linspace(0,1,n)]';
                cmap2 = [linspace(1,0,n); linspace(1,0,n);...
                    linspace(1,1,n)]';
                cmap = [cmap1; cmap2(2:end,:)];
                colormap(cmap);
                if strcmp(titl,'alphabet')
                    title(alph(1),'FontSize', 15)
                else
                    title(['Cross correlations of ' pair ' between ',...
                        dates2 ],'FontSize', 15)
                end

                subplot(2,2,2)
                plot(cld,ref)
                xlim([xaxis])
                if strcmp(titl,'alphabet')
                    title(alph(2),'FontSize', 15)
                else
                    title(['Corrected reference trace used to measure ',...
                        'timedelays'],'FontSize', 15)
                end

                subplot(2,2,4)
                plot(dd1,linear_td/Fq,dd1,timedelay0/Fq,'.')
                legend('Linearly fitted timedelay','Measured timedelay')
                axis([1 num_days yaxis])
                xlabel('Days','FontSize', 16), ylabel('Time delay (s)',...
                    'FontSize', 16)
                if strcmp(titl,'alphabet')
                    title(alph(3),'FontSize', 15)
                else
                    title(['Relative timedelay for ' pair ' between ',
                        dates2],'FontSize', 15)
                end

                tit= ['EGFs, reference and TD of ' pair ' between ',...
                    dates1 ];       
                set(h, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
                print(['figures/' tit ],'-dpng','-r0')
                close h

             elseif strcmp(varargin{1},'EGF')
                % Plot ONLY the daily Green's function 
                h=figure;
                imagesc(cld,dd1,egfn)
                tit = ['Daily EGFs for ' pair ' between ' dates2 ];
                title(tit,'FontSize', 14)
                axis([xaxis 1 num_days])
                xlabel('Time (s)','FontSize', 16), ylabel('Days',...
                    'FontSize', 16)
                % Set the colorbar to red-white-blue:
                colorbar
                level = 50; 
                n = ceil(level/2);
                cmap1 = [linspace(1,1,n); linspace(0,1,n);...
                    linspace(0,1,n)]';
                cmap2 = [linspace(1,0,n); linspace(1,0,n);...
                    linspace(1,1,n)]';
                cmap = [cmap1; cmap2(2:end,:)];
                colormap(cmap);

             elseif strcmp(varargin{1},'TD')
                % Plot ONLY the measured timedelays 
                h=figure;
                subplot(2,1,1)
                plot(cld,ref)
                xlim(xaxis)       
                if strcmp(titl,'alphabet')
                    title(alph(1),'FontSize', 15)
                else
                   title(['Corrected reference trace used to measure ',...
                       'timedelays'],'FontSize', 15) 
                end

                subplot(2,1,2)
                plot(dd1,linear_td/Fq,dd1,timedelay0/Fq,'.')
                legend('Linearly fitted timedelay','Measured timedelay')
                axis([1 num_days yaxis])
                xlabel('Days','FontSize', 16), ylabel('Time delay (s)',...
                    'FontSize', 16)        
                if strcmp(titl,'alphabet')
                    title(alph(2),'FontSize', 15)
                else
                   title(['Relative timedelay for ' pair ' between ',...
                       dates2],'FontSize', 15)
                end

                tit = ['Reference and TD for ' pair ' between ' dates2];

            elseif strcmp(varargin{1},'Frequency')
                % Plot ONLY the daily Green's function 
                freq = linspace(0,Fq,lcc);
                h=figure;
                imagesc(freq,dd1,abs(fft(egfn,lcc,2)))
                tit = ['Daily EGFs for ' pair ' between ' dates2 ];
                title(tit,'FontSize', 14)
                axis([0 Fq/2 1 num_days])
                % Set the colorbar to red-white-blue:
                colorbar

            elseif strcmp(varargin{1},'Daily')
                % Plot the daily Green's function as signals not amplitude 
                pp = 1;
                h = figure;
                for j=1:num_corr
                    egfp = egfn(j,:) + pp;

                    plot(cld,egfp)
                    set(gca,'Ydir','reverse')
                    tit = ['Daily cross correlation of ' pair,...
                        ' between ' dates2];
                    title(tit,'FontSize', 15)
                    axis([xaxis -0.1 num_corr+2])
                    xlabel('Time (s)','FontSize', 16), ylabel('Days',...
                        'FontSize', 16)
                    hold on

                    pp = pp + 1.1;
                end
                hold off      

            elseif strcmp(varargin{1},'Stack')
                % Plot the daily Green's function as signals not amplitude 
                stack = sum(EGF);

                h=figure;      
                plot(lagd,stack)
                tit=['Stacked cross correlation of ' pair ' between ',...
                    dates2];
                title(tit,'FontSize', 15)
                xlim(xaxis)
                xlabel('Time (s)','FontSize', 16), ylabel('Days',...
                    'FontSize', 16)

            elseif strcmp(varargin{1},'all')  
                    % Plot all the station pairs with TD
                    subplot(2,nsp,sp)
                    imagesc(cld,dd1,egfn)
                    axis([xaxis 1 num_days])
                    xlabel('Time (s)','FontSize', 10), ylabel('Days',...
                        'FontSize', 10)
                    % Set the colorbar to red-white-blue:
                    level = 50; 
                    n = ceil(level/2);
                    cmap1 = [linspace(1,1,n); linspace(0,1,n);...
                        linspace(0,1,n)]';
                    cmap2 = [linspace(1,0,n); linspace(1,0,n);...
                        linspace(1,1,n)]';
                    cmap = [cmap1; cmap2(2:end,:)];
                    colormap(cmap);
                    colorbar

                    if strcmp(titl,'alphabet')
                        title(alph(sp),'FontSize', 15)
                    else
                       title(pair,'FontSize', 15)
                    end

                    subplot(2,nsp,sp+nsp)
                    plot(dd1,linear_td/Fq,'b-',dd1,timedelay0/Fq,'r.')
                    legend('Linearly fitted timedelay',...
                        'Measured timedelay')
                    axis([1 num_days yaxis])
                    xlabel('Days','FontSize', 10), ylabel(...
                        'Time delay (s)','FontSize', 10)
                    hold on      

            elseif strcmp(varargin{1},'allC')  
                    % Plot all the components for each station pairs
                    subplot(2,nch,spp)
                    imagesc(cld,dd1,egfn)
                    axis([xaxis 1 num_days])
                    xlabel('Time (s)','FontSize', 10), ylabel('Days',...
                        'FontSize', 10)
                    % Set the colorbar to red-white-blue:
                    level = 50; 
                    n = ceil(level/2);
                    cmap1 = [linspace(1,1,n); linspace(0,1,n);...
                        linspace(0,1,n)]';
                    cmap2 = [linspace(1,0,n); linspace(1,0,n);...
                        linspace(1,1,n)]';
                    cmap = [cmap1; cmap2(2:end,:)];
                    colormap(cmap);
                    colorbar

                    if strcmp(titl,'alphabet')
                        title(alph(sp),'FontSize', 15)
                    else
                       title(pair,'FontSize', 15)
                    end

                    subplot(2,nch,spp+nch)
                    plot(dd1,linear_td/Fq,'b-',dd1,timedelay0/Fq,'r.')
                    legend('Linearly fitted timedelay',...
                        'Measured timedelay')
                    axis([1 num_days yaxis])
                    xlabel('Days','FontSize', 10), ylabel(...
                        'Time delay (s)','FontSize', 10)
                    hold on           

            elseif strcmp(varargin{1},'allEGF')
                    % Plot all the station pairs
                    subplot(1,round(nsp),sp)
                    imagesc(cld,dd1,egfn)
                    axis([xaxis 1 num_days])
                    xlabel('Time (s)','FontSize', 15), ylabel('Days',...
                        'FontSize', 15)
                    % Set the colorbar to red-white-blue:
                    colorbar
                    level = 50; 
                    n = ceil(level/2);
                    cmap1 = [linspace(1,1,n); linspace(0,1,n);...
                        linspace(0,1,n)]';
                    cmap2 = [linspace(1,0,n); linspace(1,0,n);...
                        linspace(1,1,n)]';
                    cmap = [cmap1; cmap2(2:end,:)];
                    colormap(cmap);  

                    if strcmp(titl,'alphabet')
                       title(alph(sp),'FontSize', 15)
                    else
                       title(pair,'FontSize', 15)
                    end

                    hold on
            end
            if length(varargin)>1 && h1~=0 && strcmp(varargin{2},'save')
                set(h, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
                print(['figures/' tit ],'-dpng','-r0')
                close h
            end
            end
        hold off
        if length(varargin)>1 && h1==0 && strcmp(varargin{2},'save')
            tit = [pair 'all_components'];
            set(h, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
            print(['figures/' tit ],'-dpng','-r0')
            close h
        end
        % Make sure the stations are cross correlted with the rigth number 
        % of stations:
        end
        if ii >= num_stat_cc
            ii = 0; 
        else
            ii = ii + 1;
        end
    end
end
