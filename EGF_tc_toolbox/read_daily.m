function dataout = read_daily(network,stationname,channels,location,datevector,fileName,fileformat,dateformat,Fq,deci,missingfiles)
% Reads the daily sac- or mseed -files and checks that the start- and endtime 
% is correct. Decimates the data, if requested
%
% Input:
%       network = network of the station to be read
%       stationname = name of the station to be read
%       channels = channel
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
%       varargout: daily=a struct giving the raw data and information about
%                       the day and station
%
%
% Sub-function: str2filename.m, 
%   rdsac.m 
%   (https://se.mathworks.com/matlabcentral/fileexchange/46356-rdsac-and-mksac-read-and-write-sac-seismic-data-file)
%   ReadMSEEDfast.m
%   (https://se.mathworks.com/matlabcentral/fileexchange/46532-readmseedfast-filename)
%
% Written by Karina Løviknes 
% 

Name = {[network '-' stationname '-' channels]};
num_days = length(datevector);
tspd = Fq*60*60*24; % Total samples per day

fe = 0; % Count the days when the file does not exist
for d = 1:num_days    
    % Extract the filename
    filename = str2filename(fileName,stationname,dateformat,'network',network,...
        'channels',channels,'location',location,'datevector',datevector(d));

    Warning_msg = {[]};

    if exist(filename,'file') % Check that the file exists
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
       
            %Check that the given samplingrate mach the samplingrate of the file:    
            if round(1/delta)~=Fq
                error('The given samplingrate is not correct') 
            end

        elseif strcmp(fileformat,'miniseed')
            % Retrive miniseed files
            file = ReadMSEEDFast(filename);

            data = double(file.data); % The data vector
            Count = file.sampleCount; % Number of points
            sps = file.sampleRate; % The sampling frequency
            
            % EXtrct the time vector if exists
            time_vector1 = file.matlabTimeVector; % Time vector 
            if isempty(time_vector1)
                time_vector_str = datetime(file.dateTimeString,'InputFormat','yyyy/MM/dd HH:mm:ss.SSS') + seconds((0:Count)/Fq);
                time_vector = datenum(time_vector_str');
                warning('The time vector is empty');
                Warning_msg = ['The time vector is empty'];
                
            else
                time_vector = time_vector1(:,1); % Time vector
            end
           
            % Check if downsampling
            if ~isempty(deci)
                % Downsample the input file:
                data=decimate(data,deci);
                time_vector=decimate(time_vector,deci);
                Count=length(data);
                sps=sps/deci;
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
             time_vector = [datenum(datetime(datestr(time_vector(1),'yyyy-mm-dd HH:MM:SS.FFF'),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS') - flip(seconds((1:missing_p)/Fq)))'; time_vector];

             timestart = time_vector(1);

             Starttime = datestr(ts1,'mmmm dd, yyyy HH:MM:SS.FFF')
             warning([num2str(missing_p) ' zeros added added to beginning of trace'])
             Warning_msg = [Warning_msg num2str(missing_p) ' zeros added '];
             
             milisec = (ms1/1000)*Fq;
             if  ~isa(milisec,'integer')
                % Interpolate to get the correct startingtime:
                nq = 1000/Fq;
                tv = (1:Count)';
                tvq = (1:Count*nq)';

                zms = zeros(ms1,1);
                data_intp1 = interp1(tv,data,tvq,'spline');
                data_newstarttime = [zms; data_intp1];
                newtime = [datenum(datetime(datestr(time_vector(1),'yyyy-mm-dd HH:MM:SS.FFF'),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS') - flip(milliseconds(1:ms1)))'; datenum(datetime(datestr(time_vector(1),'yyyy-mm-dd HH:MM:SS.FFF'),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS') + milliseconds(0:Count*nq))'];
                data = decimate(data_newstarttime,nq);
                time_vector = downsample(newtime,nq);
                Count = length(data);
                
                timestart = time_vector(1);
                
                new_startime = datestr(time_vector(1),'mmmm dd, yyyy HH:MM:SS.FFF')
                Starttime = new_startime;
                warning('Changed starttime')
                Warning_msg = [Warning_msg 'changed starttime'];
             end

         end

         if Count>=tspd
             % The numer of points is more or equal the total number of
             % samples

             num_po = Count-tspd;
             te1 = time_vector(tspd,1); % The end of recording

             hour2 = str2num(datestr(te1,'HH'));  
             min2 = str2num(datestr(te1,'MM'));
             sec2 = str2num(datestr(te1,'ss'));
             ms2 = str2num(datestr(te1,'FFF'));

             if hour2==23 && min2==59 && sec2==59 && ms2>=1000*(1-delta) && ms2<=999
                % Correct endtime
               data = data(1:tspd);
               time_vector = time_vector(1:tspd);
               Count = tspd;
               timend = te1;

             else
                 Endtime = datestr(te1,'mmmm dd, yyyy HH:MM:SS.FFF')
                 warning('Dobbelcheck endtime');
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

             if hour2==23 && min2==59 && sec2==59 && ms2>=1000*(1-delta) && ms2<=999
                 % The endtime is correct, count number must be wrong
                 timend=te1;

                 Endtime = datestr(te1,'mmmm dd, yyyy HH:MM:SS.FFF')
                 warning('Dobbelcheck count')
                 Warning_msg = [Warning_msg 'Wrong count'];
             else
                 % Time gap at the end of the day segment
                 num_po = tspd-Count;

                 mzp = zeros(num_po,1);
                 data = [data; mzp];
                 Count = length(data);
                 newtime = [time_vector; datenum(datetime(datestr(time_vector(1),'yyyy-mm-dd HH:MM:SS.FFF'),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS') + seconds((0:num_po)/Fq))'];

                 timend = time_vector(end);

                 Endtime = datestr(te1,'mmmm dd, yyyy HH:MM:SS.FFF')
                 warning([num2str(num_po) 'zeros added to end of trace'])
                 Warning_msg = [Warning_msg num2str(num_po) 'zeros added e'];
             end
         end 
         Warning_message = {Warning_msg};
        daily(d,:) = table(Name,Count,Fq,{Starttime},{Endtime},Warning_message{1});  
    else
        fe = fe+1
        if fe > missingfiles
            error('Too many missing files on row')
        end
        data = zeros(tspd,1);
        warning('The file does not exist')
        date = datevector(d)
        timestart = 0;
        timend = 0;
        Warning_msg = [Warning_msg 'Missing file']

        Warning_message = {Warning_msg};
        daily(d,:) = table(Name,Count,Fq,{Starttime},{Endtime},Warning_message{1}); 
    end
    dataout(d,:) = data';
end
% Write the daily data-information into a textfile:
writetable(daily,[stationname '_' datestr(datevector(1)) '-' datestr(datevector(end)) '.txt']);
end

