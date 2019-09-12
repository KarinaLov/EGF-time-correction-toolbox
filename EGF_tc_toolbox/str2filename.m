function filename =str2filename(filename,stationname,dateformat,varargin)
% Converts filename format from string to actual filename by replacing
% stationname, date and other variables by actual given values
%
% Input: 
%       filename = filename format
%       stationname = name of the station
%       dateformat = the format of the date given in the filename
%       varargin: Other variables included in the filename
%            
% Written by Karina Løviknes 
% 

% Stationname must be included!
snm = strfind(filename,'stationname');
if isempty(snm)
    error('Stationname must be included in the filename (or foldername)')
end

filename = strrep(filename,'stationname',stationname);


for j=1:2:length(varargin)
    var1 = char(varargin{j});
    var2 = char(varargin{j+1});
    
    % Fix the dates:
    if strcmp(var1,'datevector') 
        [yr,mm,dy,HH,MN,SS] = datevec(var2);
        % Replace the word DATE with the actual date in the right format
        filename = strrep(filename,'(DATE)',datestr(var2,dateformat));
        
        % Replace 'yyyy' with the actual year
        filename = strrep(filename,'yyyy',datestr(var2,'yyyy'));
        
         % Replace 'ddd' with the actual day
         doy = day(datetime(var2),'dayofyear');
         filename = strrep(filename,'ddd',num2str(doy,'%03d'));
        
        % If the next day is included in the name:
        dp = datetime(var2)+days(1);
        filename = strrep(filename,'(DATE+1)',datestr(dp,dateformat));
    else
        filename = strrep(filename,var1,var2);
    end
end
end