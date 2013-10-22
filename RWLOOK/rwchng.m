function y=rwchng(x,k)
% rwchng: scaled values or scaled first-difference of a time series
% rwchng(x,k);
% Last revised 5-28-03
% 
% Scaled values or scaled first-difference of a time series, as needed by skelcrn.m (see notes)
% The scaled first-difference ring-width or index time series is the
% change in the series from the previous year expressed as a decimal fraction
% of the local level of the series.  The local level is the values of the series
% smoothed by a 9-weight Hamming filter.  The steps in the transformation:<P>
%
%  1. 9-wt hamming low-pass filter the ring widths or indices<BR>
%  2. compute departures from low-pass series curve <BR>
%  3. compute first difference of departures to accentuate high-frequency change<BR>
%  4. take ratio of first diff of departures to value of the low-pass series in previous year<BR>
%
%
%*** INPUT ARGS
%
% x (mx x 1) input time series, a ring-width series or index series
% k (1 x 1)i  option for first differencing of time series (see notes)
%   ==1 no first differencing ("local" mode in skelcrn)
%   ==2 first differencing 
%
%*** OUTPUT ARGS
%
% y (my x 1)  scaled first difference of departures
%     my=mx-1
%
%*** REFERENCES -- none
%
%*** UW FUNCTIONS CALLED -- none
%
%*** TOOLBOXES NEEDED
% signal processing
%
%*** NOTES
%
%The series is expressed as the decimal fraction of the local level as defined
% by filtering.  The local level is the value of the series
% smoothed by a 9-weight Hamming filter. Depending on input argument k, 
% the series may be first differenced at the start.  The scaled first-difference ring-width or index time 
% series is the change in the series from the previous year expressed as a decimal fraction
% of the local level of the series.  
%
% filtfilt.m is used to Hamming-filter the series.  The series is filtered in both the
% forward and backwards direction to avoid phase shift.
%
% The transformation in rwchng.m is useful to build a "crossdating" series representing
% high-frequency variations adjusted to remove effects of differing level of the series.
% With k==2, extremely high frequency variation is emphasized.
%
% k = option for first differencing.  Optional.  If not included, k assumed equal to 2, which calls for first differencing.

mx=length(x);
% Design low-pass Hamming Filter
b=fir1(8,.1);  % returns a 9-weigth (n+1) filter with cutoff frequency .05 (20 yrs) 

% Filter the ring width.  Using filtfilt this way is the same as 
% convoluting b with itself and then filtering with filter and 
% having to adjust for phase shift.  filtfilt also automatically 
% handles end effects (using reflection across ends).
g=filtfilt(b,1,x);  

% Departures from smooth curves
d=x-g;


if nargin==1;
    k=2;
else
    if ~any(k==[1 2]);
        error('k not equal to 1 or 2');
    end;
end;


if k==2; %  first-differencing
    % First differences of departures
    f=diff(d);
    
    % Standardized to value of smooth curve in previous year
    gs=g(1:mx-1);
else; % no first-differencing
    f=d;
    gs=g;
    
end;


y=f ./ gs;
t=(1:length(y))';
%plot(t,x(2:length(x)),t,y,t,gs)';
