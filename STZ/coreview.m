function coreview
% coreview: graphical summary of core indices
% CALL: coreview;
%
% Meko 4-12-99
%
%******** IN
%
% No args
% User prompted for name of .mat file with core2.m output
%
%******** OUT 
% 
% No args, just graphics in Figure Windows:
%  1) Ring-width with overplotted fitted age curve 
%  2) Index by ratio method and by difference method
%
%******** NOTES
%
% Ring-width with overplotted age trend curve.  Curve may have been fit to
% subset of available years of ring-width. Plotted smooth curve covers
% only the period fit.
%
%********** UW files needed
%

close all;

% Get name of .mat file with corei.m output
flmat=uigetfile('*.mat','Input .mat file with corei.m output');
flold=flmat; % Will need this file name later
eval(['load ',flmat]);


if ~exist('S') | ~exist('nms') | ~exist('yrs') | ~exist('X'),
	error('S, nms, yrs, or X not in .mat file')
end

kwh1 = 1; % while control

while kwh1==1;
   kmen1 = menu('Choose one','Time series plots for single core','Quit');
   switch kmen1;
   case 1; % time series plots for single core
      % Choose core
      kcore = menu('Choose core',cellstr(nms));
          
      % Pull ringwidth series
      yrgo1 = yrs(kcore,1);  yrsp1 = yrs(kcore,2);  % start, end year for ring width
      nyr1 = yrsp1-yrgo1+1; % number of years of ringwidt
      nmcore =nms(kcore,:); % core name
      igo1= yrs(kcore,3); % start row of ring-width in X
      isp1= igo1 + nyr1 -1;
      x1 = X(igo1:igo1+nyr1-1); % ring width
      yr1 = (yrs(kcore,1):yrs(kcore,2))'; % year vector for ring-width
      
      % Plot ring-width
      figure(1);
      hp1 = plot(yr1,x1/100);
      title(nmcore);
      xlabel('Year');
      ylabel('Ring Width (mm)');
      grid;
      
      % Overlay fitted growth curve
      kcurve=questdlg('Overlay growth curve?');
      switch kcurve;
      case 'Yes';
         nfit1=S(kcore,2); % curve-type for detrending
         
         % Compute year vector for growth curve
         igo2 = igo1 + S(kcore,10)-yrgo1; % index to X of first year used in fitting
         isp2 = isp1 - (yrsp1 - S(kcore,11)); % and to last year
         yrgo2 = yrgo1 + (igo2-igo1);
         yrsp2 = yrsp1 - (isp1-isp2);
         yr2 = (yrgo2:yrsp2)';
         nyr2 = isp2-igo2+1;
         x2 = X(igo2:isp2);
         
         if nfit1==9; % Spline -- in this case, wil re-fit 
            g = (cfspl(S(kcore,3),yr2,yr2,x2))';
         elseif nfit1==1; % negative exponential -- generate using stored parameters
            kne=S(kcore,3);
            ane=S(kcore,4);
            bne=S(kcore,5);
            t = yr2-yr2(1) + 1;
            g = (kne + ane * exp(-bne*t));
            % Fits data  to the equation
            % g(t) = k + a * exp(-b*t),
            %   where t is the shifted time variable t = yrvn-yrvn(1)+1
            %   In other words, t is same length as yrvn after 
            %   dropping NaNs
            
         elseif nfit1==4; % Horizontal line through mean
            g= repmat(S(kcore,3),nyr2,1); 
         elseif nfit1==2; % straight line, any slope
            t=yr2-yr2(1)+1;
            Xtemp = [ones(length(t),1) t];
            aconst = S(kcore,3);
            aregr = S(kcore,4);
            g  = Xtemp * [aconst aregr]';
         else;
            error('Unnacceptable value for nfit1');
            
         end
         
         hold on;
         hp2=plot(yr2,g/100);
         set(hp2,'LineWidth',2,'Color',[.5 .5 .5]);
         hold off
         
         kindex = questdlg('Index Plot in Window 2?');
         switch kindex;
         case 'Yes';
            wrat=x2 ./g;
            
            % Compute 'difference' index, forcing same stand dev as ratio index and
            % same mean;
            wsd = std(wrat);
            wmean =nanmean(wrat);
            wtemp = x2-g;
            wtemp = zscore(wtemp);
            wdif = wtemp * wsd + wmean;
            figure(2);
            plot(yr2,wrat,yr2,wdif);
            grid;
            title(nmcore);
            legend('Ratio','Difference');
            xlabel('Year');
            ylabel('Tree-Ring Index');
         case 'No';
         case 'Cancel';
         end; % switch kindex
      case 'No';
      case 'Cancel';
      end; % switch
   case 2; % Quit
      kwh1=0;
   end; % switch kmen1;
end; % while kwh1==1


