function [p,per] = monotspl(yrv,yrvn,xvn,N,kopt)
%
% Fit non-increasing spline with smallest possible p-parameter.  This
% is the most flexible possible spline that is still non-increasing.
%
% D. Meko 6-20-95
%
%**********  IN ARGS 
%
% yrv -- years at which smoothed values are needed
% yrvn -- years at which original time series has data;  some might
%	be NaN
% xvn -- original time series;  some might be NaN
%
% cvx -- values of smooth curve at years yrv
%
% N -- length of original full-length time series
%
% kopt -- 1== display iteration data
%         2==suppress
%
%***************  OUT ARGS
%
% p -- spline parameter
% per -- period at which the returned spline has ampl freq resp of 0.5 
%
%
%**********  METHOD  
%
% An initial fit is carried out with the 0.5 Amplitude freq response (AFR)
% at 2N, where N is the length of the core time series.  The smooth curve
% from this fit will be monotonic, but likely too stiff.  The factor 
% multiplying N is then changed by an amount specified in one of the
% the elements in step to give a more flexible spline (e.g., 1.5N). 
% The smooth curve is tested to determine whether it is non-increasing over
% all parts.  The curve is successively made more flexible until
% the test for non-increase with time if failed.  The parameters are
% then re-set to the previous values (giving the last non-increasing fit),
% the step size is decreased, and the fitting loop is covered again.
% The final fit is accepted when increasing the spline parameter p
% by decreasing the multiplier of N using the smallest step size would
% result in a trend line that would fail the test for non-increasingness.
%
% If 0.5 AFR spline with period 2N is monotonic increasing, no
% further iteration is done.



yrvn(isnan(xvn))=[];  % delete rows with xvn as NaN
xvn(isnan(xvn))=[];

nobs=length(yrvn);
ihalf=round(nobs/2);


% Hard-code the steps for changing the fraction of sample length.
step = [0.5 0.2 0.1 .05 .02];
ns = length(step);

k1=1; %while control for still fitting
a = 2;  % begin with spline with ampl freq respons of 0.5 at 
	% a=2 times the sample length
i=1;
s = step(i);  % First will change the parameter in large steps

ifirst=0;  % flag that means have not yet got an ititial fit satisfying conditions

while k1;  % while still trying to find suitable fit
  p = splinep(a*N,0.5); % set spline parameter
  c = csaps(yrvn,xvn,p,yrv); % fit spline
  d =  diff(c);  %first difference of smoothed curve 
  d2=diff(d); % second difference of smoothed curve
  L1 = (d <=0.0); % will be 1 if curve never  increasing with time
  L2=d2(ihalf:(nobs-2))>-1e-5; % 1 if change in slope during second half of 
    % positive or imperceptibly negative
  if all(L1) & all(L2);  % curve non-increasing and slope becoming less steep with
     % time over second half of record
     
    a = a - s;  %decrease the spline parameter (make more flexible)
    per = a* N;  % period of 0.5 amp.resp
    ifirst=1;
    txt=[num2str(per),'  ',num2str(a)];
    if kopt==1;
       disp(txt);
    end;
    
  else;  % curve increases somewhere (too flexible), or possible
    % everywhere (monotonic increasing growth), or slope steepens in second half of record
    if ~any(L1);  % curve is monotonic increasing
       per = a*N;
       
       if kopt==1;
          disp('MONOTONIC INCREASING GROWTH TREND');
          disp(['Spline has 0.5 AFR at ',int2str(a*N),' years']);
          disp('(which is twice the series length)');
       end;
       return
    end
    if s > step(ns);  % Have not yet iterated to smallest step size in p
       a = a+ s;  % go back to previous (non-increasing) spline p
       if ifirst==1; % if have already had a fit satisfying conditions
          i = i +1;  % decrease step size
       else;
          % No change to step size until need to back up to less flexible curve
          
       end;
       s = step(i);
    else;  % bingo! Using smallest step size, have gone to too-flexible spl
	 	a = a + s; % go back a step
	 	k1 = 0;  % exit the while loop
    end
  end
  % Handle case where non-increasing conditions have been met but
  % have iterated to less than 0 times sample length
  if a<= 0;
	  a = a + s;
	  i = i+1;
     s = step(i);
  else
  end
end % of while k1

per = a*N;
p = splinep(per,0.5);
if kopt==1;
    disp(['p=',num2str(p)]);
    txt=[num2str(per),'  ',num2str(a)];
    disp(txt)
else;
end;
