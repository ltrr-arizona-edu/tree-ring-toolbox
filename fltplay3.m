function [b,p50,y,yry]=fltplay3(x,yr,prevwind,kopt)
% fltplay3:  trial and error fits of low-pass binomial filter
% CALL: [b,p50,y,yry]=fltplay3(x,yr,prevwind,kopt);
%
% D. Meko 10-14-00
%
%*******  INPUT ARGS
%
% x time series, col vect
% yr corresp years
% prevwind previous window to preserve in calling fltplay2
% kopt(1 x 1)i   options
%   kopt(1)  type of filter
%     ==1 Gaussian
%     ==2 Binomial
%
%**********  OUTPUT ARGS
%
% b   (mb x 1)  filter weights of mb-length filter
% p50  (1x1) period (yr) of 50% freq response of filter
% y  (my x 1) filtered time series
% yry (my x 1)  years for filtered time series
%
%************  COMMENTS **********************
%
% Usually run before filter1.m or firdmm2.m to arrive at settings for 
% period of 50% amplitude freq response and number of weights in filter
%
% Running instructions:
%  -Give the function call
%	-If first filter tried, choose "Design or Revise filter"
%	-Give target period (e.g., 10 (years)) that you want the filter to have
%		ampl of frequency response equal to 0.5;  and give desired number of
%		weights (e.g., 15) in filter
%	-Function computes filter weights, frequency response, etc.
%	-Choose button " Review... ", then click on figure windows
%		to see (1) time series plot with filtered version overplotted, 
%		(2) plot of magnitude of frequency response, and (3) a plot of the
%		filter weights
%	-If satisfied, choose "Accept", and function will have returned the
%		output arguments for that filter (see above)
%	-Or choose "Design and revise filter" and go round again.  

switch kopt(1);
case 1; % Gaussian
   ftype='Gaussian';
case 2; % Binomial
   ftype='Binomial';
otherwise;
   error('kopt(1) must be 1 or 2');
end;

pdfirst=10; % initial Hard code setting for period of 50% response
kfirst=1;

k1=1;
while k1~=3;
	k1=menu('Choose One',['Design or Revise ' ftype ' Filter'],...
		'Review Filter',...
		'Accept Filter');

   if k1==1; % Try a new filter
      
      
      prompt={'Period (yr) at which Amp of Frequency Response is 0.5 :'};
      if kfirst==1;
         deflt={num2str(pdfirst)};
      else;
         deflt={num2str(pdprev)};
      end;
      titdlg1='Enter Desired Filter Properties';
      lineNo=1;
      answer=inputdlg(prompt,titdlg1,lineNo,deflt);
      p50 = str2num(answer{1});
      kfirst=0;
      pdprev=p50;
      %nl = str2num(answer{2});
      if p50<3;
         error('Period of 50% response must be at least 3 yr');
      end;
      
      % Compute filter weights
      switch ftype;
      case 'Gaussian';
         b=wtsgaus(p50); % compute filter weights
      case 'Binomial';
         b=wtsbinom(p50,1);
      end;
      
      nl=length(b);  % length of filter
          
      Wn=2.0 / p50;  % 50% response frequency on scale such that 0 is 
      % 0 per year, and Wn=1 is 0.5 per year (Nyquist).  At
      %  frequency Wn, response of filter b is 0.5
      
      koptf=1;  % filter1 option saying not data extension 
      [y,yry]=filter1(x,yr,b,koptf); % apply filter to series
      
      %y1=filter(b,1,x);  % 
      
      % Get length of filtered series
      ny=length(y);
      
      % Cull period of overlap of original and filtered series
      % for variance comparison
      L3=yr>=yry(1) & yr<=yry(length(yry));
      z=x(L3);
      
      % Compute variance ration
      vary=std(y)*std(y); % variance of smoothed series
      varz=std(z)*std(z); % variance of original series
      
      %pause
      ratvar=vary/varz; % variance ratio
      
   elseif k1==2; % review plots and filter weights
      figure(1+prevwind); % reserve for time series plots
      h4p=plot(yr,x,yry,y);
      %title(['Specified ' int2str(nl) ' weights, 50% response at ',...
      %     num2str(p50) ' years']);
      title('Original and Filtered Time Series');
      pltext(.1,.95,12,'\itx\rm = original,  \ity\rm = filtered');
      pltext(.1,.9,12,['var(\itx\rm) = ',num2str(varz)]);
      txtrat = sprintf('%5.3f',ratvar);
      pltext(.1,.8,12,['var(\ity\rm)/var(\itx\rm)  = ',txtrat]);
      pltext(.1,.85,12,['var(\ity\rm = ',num2str(vary)]);
      set(h4p(2),'Linewidth',1.5);
      
      % Plot freq response
      figure(2+prevwind);
      [w1,A,pv]=freqres1(b,100);
      text(.22,.75,['Final ' int2str(nl),'-Weight Filter']);
      text(.20,.95,'Specifications');
      %text(.25,.90,['Number of weights = ' int2str(nl)]);
      text(.25,.90,['50% Response at ' num2str(p50) ' years']);
      title(['Frequency Response of ' ftype  ' Filter'])
      
      % Plot weights
      figure(3+prevwind);
      plot(b);
      title(['Filter Weights; ',int2str(nl),' Weights']);
      %clc
      %disp('Filter Weights');
      %disp(' ')
      %disp(b);
      %disp('Switch to Figure Windows for Graphics')
   elseif k1==3;
      break;
   end ;  % of if
   
end;  % of while
   