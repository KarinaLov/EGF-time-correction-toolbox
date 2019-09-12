function [cross,varargout] = cross_conv(seq1,seq2,Fq,varargin)
  % Cross correlation function using flip and conv (convolution)
  %
  % Input: 
  %     seq1,se2 = the input segments to be cross correlated
  %     Fq = sampling frequency of the input signals
  %     varargin:
  %          wl = window length to cross correlate over in hours 
  %          swl = number of hour to stack the cross correlations over 
  %          perco = percent overlap
  %
  % Output:
  %     cross = the cross correlation of seg1 and seg 2
  %     varargout:
  %         lag = time lapse of the cross correlation
  %         faulty = vector or array giving the days were the mean are
  %                  NaN (not a number)
  %
  % Written by Karina Løviknes
  %
  
  m = length(seq1);
  n = length(seq2);
  faulty = 0;
  
  if isempty(varargin)
    % Default: the input segments are cross correlated over the entire day 
    % using flip and conv  
    s1 = flip(seq1);
    s2 = seq2;
        
    if (sum(s1)==0 & mean(s1)==0) | (sum(s2)==0 & mean(s2)==0)
        % One of the segments are completly zero and the cross
        % correlation is not necessary  
        warning('One of the segments are completly zero, the CC is set to zero')
    else
        cross = conv(s1,s2);
    end
    
    % Make the lag vector start from negative and have length (m+n)-1:
    k = max(m,n);
    lag = zeros(1,(m+n)-1);
    for i = 1:(m+n)-1
        lag(i) = -k+i;
    end
    
  else
    wl = varargin{1};
    swl = varargin{2};
    perco = varargin{3};
    % The input segments are cut to shorter time windows before cross 
    % correlation to save time
        
    po = perco/100;
    nwl = wl*60*60*Fq; % Window length in samples
    nuwl = (swl/wl)*round(1/po); % Number of windows per day 
    if po<1
        nuwl = nuwl-1;
    end
    
    % Preallocate for speed:
    cross_seg = zeros(nuwl,2*nwl-1); 
    cross = zeros(24/swl,2*nwl-1);
    
    i=0;
    ii=1;
    for l = 1:24/swl
    for j = 1:nuwl
        s1 = seq1(1,((nwl*i)+1):(nwl*ii));
        s2 = seq2(1,((nwl*i)+1):(nwl*ii));

        if (sum(s1)==0 & mean(s1)==0) || (sum(s2)==0 & mean(s2)==0)
            % One of the segments are completly zero and the cross
            % correlation is not necessary
        else
            sf1 = flip(s1);
            cross_seg(j,:) = conv(sf1,s2);     
        end
        
        i = i+po;
        ii = ii+po;
    end
        cross_seg_m = cross_seg(2:end-1,:);
        cross(l,:) = sum(cross_seg_m);
          % If the mean of the cross correlation is NaN the array is set to zero to avoid afffecting the final stack:
          if isnan(mean(cross(l,:)))
              cross(l,:) = zeros(1,length(cross)); 
              faulty = 1;
          end
    end
    
    % Make the lag vector start from negative and have length (m+n)-1:
    mm = length(s1);
    nn = length(s2);
    k = max(mm,nn);
    lag = zeros(1,(mm+nn)-1);
    for i = 1:(mm+nn)-1
        lag(i) = -k+i;
    end
      
  end 
  
  % Outputs:
  if nargout == 2
      varargout{1} = lag;
  elseif nargout == 3
      varargout{1} = lag;
      varargout{2} = faulty;
  end  
end