function [b,p50,y,yry]=fltplay2(x,yr,prevwind)
% fltplay2:  trial and error fits of FIR filters to be used with firdmm2.m or filter1.m 
% CALL: [b,p50,y,yry]=fltplay2(x,yr,prevwind);
%
% D. Meko 2-26-93
%
%*******  INPUT ARGS
%
% x time series, col vect
% yr corresp years
% prevwind previous window to preserve in calling fltplay2
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

kfirst=1;
pdfirst=10; % initial period of 50% response
nfirst=15; % initial length of filter


k1=1;
while k1~=3;
	k1=menu('Choose One','Design or Revise a Filter',...
		'Review Filter',...
		'Accept Filter');

   if k1==1; % Try a new filter
      
      
      prompt={'Period (yr) at which Amp of Frequency Response is 0.5 :',...
            'Number of Weights in Filter (odd):'};
      if kfirst==1;
         deflt={num2str(pdfirst),num2str(nfirst)};
      else;
         deflt={num2str(pdprev),num2str(nprev)};
      end;
      
      titdlg1='Enter Desired Filter Properties';
      lineNo=1;
      answer=inputdlg(prompt,titdlg1,lineNo,deflt);
      p50 = str2num(answer{1});
      nl = str2num(answer{2});
      kfirst=0;
      pdprev=p50;
      nprev=nl;
      
      if p50<4 & nl<9;
         error('For 50% response <4yr, Choose number of weights at least 9');
      end;
      if p50<=5 & nl<7;
         error('For 50% response of near 5 yr, choose number of weights at least 7');
      end;
      %if nl<=p50;
       %  error('Choose number of weights greater than desired 50% response period');
     % end;
      if rem((nl+1),2)~=0; 
         error('Number of weights should be odd');
      end;
      
      
         
      n=nl-1;  
      Wn=2.0 / p50;  % Cutoff frequency on scale such that 0 is 
      % 0 per year, and Wn=1 is 0.5 per year (Nyquist).  At
      % this frequency, response of filter b is 0.5
      b=fir1(n,Wn); % Compute the filter weights
      y1=filter(b,1,x);  % filtered series,
      %before shifting and truncating

			% Compute length 
			ny=length(yr)-length(b)+1;
			yrinc=(0:ny-1)';
			yr1=yr(1)+(length(b)-1)/2;
			yry=yr1+yrinc;

			y=y1;
			y(1:length(b)-1)=[];


			L3=yr>=yry(1) & yr<=yry(length(yry));
			z=x(L3);
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
         text(.25,.90,['Number of weights = ' int2str(nl)]);
         text(.25,.85,['50% Response at ' num2str(p50) ' years']);
         title(['Frequency Response of Filter (Windowing-Method; Hamming)'])
		
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