function varargout = read_daily(network, stationname, channels, location,...
    datevector, fileName, fileformat, dateformat, Fq, deci, missingfiles)
% Reads the daily sac- or mseed -files and checks that the start- and endtime 
% is correct. Decimates the data, if requested
%
% Input:
%       network = network of the station to be read
%       stationname = name of the station to be read
%       channels = channel
%       location = location
%       datevector = vector containing the dates to be read 
%       fileName = filename format for the files to be read
%       fileformat = fileformat, either sac or miniseed
%       dateformat = the date format used in the filename, default is 'yyyy-mm-dd'
%       Fq = sampling frequency
%       deci = specifies if the data should be downsampled and by how much
%       missingfiles = maximum number of allowed missingfiles on a row
%
% Output:
%       dataout = the raw data
%       varargout: daily = a struct giving the raw data and information about
%                       the day and station
%       The daily data-information is also written into a textfile with nameformat: 
%       [stationname '_' startdate '-' enddate '.txt']
%
%
% Sub-function: str2filename.m, 
%   rdsac.m 
%   (https://se.mathworks.com/matlabcentral/fileexchange/46356-rdsac-and-mksac-read-and-write-sac-seismic-data-file)
%   ReadMSEEDfast.m
%   (https://se.mathworks.com/matlabcentral/fileexchange/46532-readmseedfast-filename)
%
% Written by Karina L??viknes 
% 

Name = {[network '-' stationname '-' channels]};
num_days = length(datevector);
tspd = Fq*60*60*24; % Total samples per day

fe = 0; % Count the days when the file does not exist
for d = 1:num_days    
    % Extract the filename
    filename = str2filename(fileName, stationname, dateformat,...
        'network', network, 'channels', channels, 'location', location,...
        'datevector',datevector(d));
    disp(['Requesting file ', filename])

    Warning_msg = {[]};

    if java.io.File(filename).exists  % Check that the file exists        
        fe = 0;
        if strcmp(fileformat,'sac')
            % Retrieve a sac file
            filename
            file = rdsac(filename);

            data = file.d; % The data vector
            time_vector = file.t; % Time vector
            Count = file.HEADER.NPTS; % Number of points  
            delta = file.HEADER.DELTA; % Sampling time intervall
            
            if ~isempty(deci)
                % Downsample the input file:
                data = downsample(data,deci);
                time_vector = decimate(time_vector,deci);
                Count = length(data);
                delta = round(1/(delta/deci));
            end
       
            % Check that the given samplingrate matches the samplingrate of 
            % the file:    
            if round(1/delta)~=Fq
                error('The given samplingrate is not correct') 
            end

        elseif strcmp(fileformat,'miniseed')
            % Retrive miniseed files
            file = ReadMSEEDFast(filename);
            
            data = double(file(1).data); % The data vector
            Count = file(1).sampleCount; % Number of points
            sps = file(1).sampleRate; % The sampling frequency

            % Extrct the time vector if exists
            time_vector0 = file(1).matlabTimeVector; % Time vector 
            if isempty(time_vector0)
                time_vector_str = datetime(file(1).dateTimeString,...
                    'InputFormat','yyyy/MM/dd HH:mm:ss.SSS') + seconds(...
                    (0:Count)/Fq);
                time_vector = datenum(time_vector_str');
                warning('The time vector is empty');
                Warning_msg = ['The time vector is empty'];

            else
                time_vector = time_vector0(:,1); % Time vector
            end

            % Extract the starting time:
            ts0(1) = time_vector(1); % Start of recording
            te0(1) = time_vector(end); % Start of recording
            
            % loop over parts of a daily file
            for j = 2: length(file)
                
                data1 = double(file(j).data); % The data vector
                Count1 = file(j).sampleCount; % Number of points

                % Extrct the time vector if exists
                time_vector01 = file(j).matlabTimeVector; % Time vector 
                if isempty(time_vector0)
                    time_vector_str1 = datetime(file(j).dateTimeString,...
                        'InputFormat','yyyy/MM/dd HH:mm:ss.SSS') +...
                        seconds((0:Count1)/Fq);
                    time_vector1 = datenum(time_vector_str1');
                else
                    time_vector1 = time_vector01(:,1); % Time vector
                end

                % Extract the starting time:
                ts0(j) = time_vector1(1); % Start of recording
                te0(j) = time_vector1(end); % Start of recording

                % Find the time difference between the end of the last file and
                % the beginning of the current file:
                [Ys0,Ms0,Ds0,Hs0,MNs0,Ss0] = datevec(ts0(j));
                [Ye0,Me0,De0,He0,MNe0,Se0] = datevec(te0(j-1));

                tsdiff = ((Hs0-He0)*60*60)+((MNs0-MNe0)*60)+(Ss0-Se0);
                fsdiff = round(tsdiff * Fq);
                if fsdiff == 1
                    data = [data; data1];
                    time_vector = [time_vector; time_vector1];

                elseif fsdiff > 1
                    zsp = zeros(round(fsdiff),1);
                    data = [data; zsp; data1];
                    time_vector = [time_vector; zsp; time_vector1];

                elseif (fsdiff-tsdiff*Fq) ~= 0 && Fq<1000
                    % Interpolate:
                    msdiff = round(tsdiff *1000);
                    disp(['Time difference between end of last and start ',...
                        'of current file: ', num2str(msdiff)])
                    nq = 1000/Fq;
                    tv = (1:Count1)';
                    tvq = (1:Count1*nq)';

                    zmsp = zeros(msdiff,1);
                    data_intp11 = interp1(tv,data1,tvq);
                    data_intp2 = [zmsp; data_intp11];
                    data = [data; decimate(data_intp2,nq)];

                    newtime1 = [datenum(datetime(datestr(time_vector1(1),...
                        'yyyy-mm-dd HH:MM:SS.FFF'),'InputFormat',...
                        'yyyy-MM-dd HH:mm:ss.SSS') - flip(milliseconds(1:msdiff)))';...
                        datenum(datetime(datestr(time_vector1(1),...
                        'yyyy-mm-dd HH:MM:SS.FFF'),'InputFormat',...
                        'yyyy-MM-dd HH:mm:ss.SSS') + milliseconds(0:Count1*nq))'];
                    time_vector1 = downsample(newtime1,nq);
                end
                Count = length(data);        
            end
            
            % Check if downsampling
            if ~isempty(deci)
                % Downsample the input file:
                data1 = decimate(data,deci);
                time_vector = decimate(time_vector,deci);
                Count = length(data);
                sps = sps/deci;
            end    

            % Check that the given samplingrate mach the samplingrate of the file:    
            if sps ~= Fq
                error('The given samplingrate is not correct') 
            end
            delta = 1/Fq; % Sampling time interval  
                
        end
        
        % Extract the starting time:
        ts1 = time_vector(1,1); % Start of recording
        te1 = time_vector(end); % Start of recording
        Starttime = datestr(ts1,'mmmm dd, yyyy HH:MM:SS.FFF');
        Endtime = datestr(te1,'mmmm dd, yyyy HH:MM:SS.FFF');
        disp(['Starttime: ', Starttime, ' --> Endtime: ', Endtime])

        hour1 = str2num(datestr(ts1,'HH'));
        min1 = str2num(datestr(ts1,'MM'));
        sec1 = str2num(datestr(ts1,'ss'));
        ms1 = str2num(datestr(ts1,'FFF'));

        if hour1==0 && min1==0 && sec1==0 && ms1==000
            % The day segment starts at correct time (00.00.00)
           timestart = ts1;
        else
            % Time gap, samples are missing
             missing_p = floor(Fq*((hour1*60*60)+(min1*60)+sec1+(ms1/1000)));

             % Fill the gaps with zeros and correct the count
             % number and time vector
             mzp = zeros(missing_p,1);
             data = [mzp; data];
             Count = length(data);
             time_vector = [datenum(datetime(datestr(time_vector(1),...
                 'yyyy-mm-dd HH:MM:SS.FFF'), 'InputFormat',...
                 'yyyy-MM-dd HH:mm:ss.SSS') - flip(seconds(...
                 (1:missing_p)/Fq)))'; time_vector];

             timestart = time_vector(1);

             Starttime = datestr(ts1,'mmmm dd, yyyy HH:MM:SS.FFF');
             warning([num2str(missing_p) ' zeros added added to beginning of trace'])
             Warning_msg = [Warning_msg num2str(missing_p) ' zeros added '];
             
             milisec = (ms1/1000)*Fq;
             if  ~isa(milisec,'integer') & Fq<1000
                % Interpolate to get the correct startingtime:
                nq = 1000/Fq;
                tv = (1:Count)';
                tvq = (1:1/nq:Count)';

                zms = zeros(ms1,1);
                data_intp1 = interp1(tv,data,tvq,'spline');
                data_newstarttime = [zms; data_intp1];
                newtime = [datenum(datetime(datestr(time_vector(1),...
                    'yyyy-mm-dd HH:MM:SS.FFF'), 'InputFormat',...
                    'yyyy-MM-dd HH:mm:ss.SSS') - flip(milliseconds(1:ms1)))';...
                    datenum(datetime(datestr(time_vector(1),...
                    'yyyy-mm-dd HH:MM:SS.FFF'), 'InputFormat',...
                    'yyyy-MM-dd HH:mm:ss.SSS') + milliseconds(0:Count*nq))'];
                data = downsample(data_newstarttime,nq);
                time_vector = downsample(newtime,nq);
                Count = length(data);
                
                timestart = time_vector(1);
                
                new_startime = datestr(time_vector(1),'mmmm dd, yyyy HH:MM:SS.FFF')
                Starttime = new_startime;
                warning('Changed starttime')
                Warning_msg = [Warning_msg 'changed starttime'];
             end

        end

         tspd=tspd;
         if Count>=tspd
             % The numer of points is more or equal the total number of
             % samples

             num_po = Count-tspd;
             te1 = time_vector(tspd,1); % The end of recording

             hour2 = str2num(datestr(te1,'HH'));  
             min2 = str2num(datestr(te1,'MM'));
             sec2 = str2num(datestr(te1,'ss'));
             ms2 = str2num(datestr(te1,'FFF'));

             if hour2==23 && min2==59 && sec2==59 && ms2>=1000*(1-delta)...
                     && ms2<=999
                % Correct endtime
               data = data(1:tspd);
               time_vector = time_vector(1:tspd);
               Count = tspd;
               timend = te1;

             else
                 Endtime = datestr(te1,'mmmm dd, yyyy HH:MM:SS.FFF');
                 warning('Doublecheck endtime');
                 Warning_msg = [Warning_msg 'Wrong endtime'];

                 data = data(1:tspd);
                 time_vector = time_vector(1:tspd);
                 Count = tspd;
                 timend = time_vector(end);
             end
         else
             % Count is less than total samples pr day
             te1 = time_vector(Count,1);

             hour2 = str2num(datestr(te1,'HH'));  
             min2 = str2num(datestr(te1,'MM'));
             sec2 = str2num(datestr(te1,'ss'));
             ms2 = str2num(datestr(te1,'FFF'));

             if hour2==23 && min2==59 && sec2==59 && ms2>=1000*(1-delta)...
                     && ms2<=999
                 % The endtime is correct, count number must be wrong
                 timend=te1;

                 Endtime = datestr(te1,'mmmm dd, yyyy HH:MM:SS.FFF')
                 warning('Dobbelcheck count')
                 Warning_msg = [Warning_msg 'Wrong count'];
                                  
                 % Time gap at the end of the day segment
                 num_po = tspd-Count;

                 mzp = zeros(num_po,1);
                 data = [data; mzp];
             else
                 % Time gap at the end of the day segment
                 num_po = tspd-Count;

                 mzp = zeros(num_po,1);
                 data = [data; mzp];
                 Count = length(data);
                 newtime = [time_vector; datenum(datetime(datestr(...
                     time_vector(1), 'yyyy-mm-dd HH:MM:SS.FFF'),...
                     'InputFormat','yyyy-MM-dd HH:mm:ss.SSS') + seconds(...
                     (0:num_po)/Fq))'];

                 timend = time_vector(end);

                 Endtime = datestr(te1,'mmmm dd, yyyy HH:MM:SS.FFF')
                 warning([num2str(num_po) 'zeros added to end of trace'])
                 Warning_msg = [Warning_msg num2str(num_po) 'zeros added e'];
             end
         end 
         Warning_message = {Warning_msg};
         {Warning_message{1,1:end}};
         daily(d,:) = table(Name,Count,Fq,{Starttime},{Endtime},...
             {Warning_message{1,1:end}});  
         dinfo(d) = struct('Name',Name,'Count',Count,'Frequency',Fq,...
             'Starttime',Starttime,'Endtime',Endtime);
    else
        fe = fe+1;
        disp(['Missing files in a row: ', num2str(fe)])
        if fe > missingfiles
            error('Too many missing files on row')
        end
        data = zeros(tspd,1);
        warning('The file does not exist')
        date = datevector(d);
        Count = tspd;
        Starttime = 0;
        Endtime = 0;
        Warning_msg = [Warning_msg 'Missing file'];

        Warning_message = {Warning_msg};
        daily(d,:) = table(Name,Count,Fq,{Starttime},{Endtime},...
            {Warning_message{1,1:end}}); 
        dinfo(d) = struct('Name', Name, 'Count', Count, 'Frequency', Fq,...
            'Starttime', Starttime, 'Endtime', Endtime);
    end
    lddd=size(data);
    dataout(d,:) = data(1:tspd)';
end
if nargout==1
    varargout{1} = dataout;
elseif nargout==2
    varargout{1} = dataout;
    varargout{2} = dinfo;
end
% Write the daily data-information into a textfile:
writetable(daily,[stationname '_' datestr(datevector(1)) '-' datestr(...
    datevector(end)) '.txt']);
end

