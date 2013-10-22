function slidets
% slidets: slide a floating time series along a reference series to find best match
% slidets;
% Last revised 7-24-99
%
% A selected segment of a "floating" time series is slid year-by-year along a reference 
% time series and the correlation coefficient is computed at each iteration. The best 
% match is defined as that with the highest correlation. An optional bootstrap analysis 
% can be used to check how likely this maximum correlation is by chance alone.<P>
%
% Slidets.m was written to help in crossdating remnant-wood ring-width series against 
% very long tree-ring chronologies from a nearby site.<P>
%
% You get to select the segment of the floating series to slide interactively, with 
% guidance from time series plots of the floating series and variance of the floating 
% series.<P>
%
%*** INPUT
%
% Prompts for names of files containing the reference series. 
% Prompt for width of window for sliding segment of floating series. 
% Prompt for ending year of sliding segment.
% Prompt for bootstrap or not
% Prompt for number of bootstrap samles
%
%*** OUTPUT
%
% Screen plots summarizing results.  Key windows are graph of correlation coefficient 
% against year in reference series, and boxplots showing distribution of correlations.  
% Also a boxplot showing distribution of maximum correlations from bootstrapping.
%
%*** REFERENCES -- none
%*** UW FUNCTIONS CALLED 
% corrone
% crn2vec2
% rwchng
% rwread3
%
%*** TOOLBOXES NEEDED
% Signal processing
%
%*** NOTES
%
% To emphasize similarity in year-to-year changes as opposed to trends, both series
% are filtered at the start into scaled first differences.  The series are filtered 
% backwards and forwards with a 9-weight hamming filter to derive smoothed time 
% series.  The first difference of the original time series is divided by the 
% corresponding value of the smoothed series to get the scaled first difference. 
%
% You select the segment of the floating series to be slid along the reference 
% series.  This series is slid one year at a time along the reference series and
% the correlation is computed at each step.  The correlation is plotted as a function 
% of ending year (in the reference series) of the alignment.  A boxplot of the 
% correlations gives some idea of variability of the correlation coefficient for 
% various alignments.  
%
% A plot of variance of the floating series in a sliding window is given to help in 
% choosing the segment of the floating series to test.  The notion is that a period of 
% high variability is a good candidate for getting a definitive match.
%
% Bootstrapped significance level of maximum correlation is available on option.  
% This option takes a while to run.  The bootstrapping can be done once you 
% have found the maximum sliding correlation of the floating segment with the 
% reference series.  A user-specified number of bootstrap samples are drawn from the 
% floating segment (e.g., 1000 samples). Each bootstrap sample is slid along the 
% reference series and the maximum correlation is recorded.  The distribution of 
% maximum bootstrap correlations can then be compared to the maximum correlation 
% found for the original (non-bootstrapped) data. 



% Notation
% variable 1: floating series
% variable 2: reference series
% u full length original series
% v windowed original series
% w windowed scaled first-diff series
% s standard-deviation of rw change, floating series

close all;

% Hard Code
defwidth=50; % default window width
dwin=defwidth; % initialize window width

%--- INPUT THE TIME SERIES

%-- Floating Series -- specify type of file
ksrc=menu('Choose Type of Source File for Floating Series',...
   'Index in ITRDB-formatted .crn file',...
   'Ring width in single .rw file','Two-column time series mtx, year as col 1');
if ksrc==1;
   src = 'crn';
elseif ksrc==2;
   src = 'rw';
elseif ksrc==3; 
   src = 'dat';
else;
   error('ksrc must be 1, 2 or 3');
end
clear ksrc;

%---Floating Series, get the file name
if strcmp(src,'crn');
   [file1,path1]=uigetfile('*.crn','Infile with chronology index');
elseif strcmp(src,'rw');
   [file1,path1]=uigetfile('*.rw','Infile with measured ring width series');
elseif strcmp(src,'dat');
   [file1,path1]=uigetfile('*.dat','Infile of 2-col time series matrix');
end;
pf1=[path1 file1]; % merge path and filename

%--- Floating Series, read the input time series;  
if strcmp(src,'crn');
   [x,s,yr]=crn2vec2(pf1); %x is index, s is sample size, yr is year
elseif strcmp(src,'rw');
   [X,guy,day,pf1]=rwread3(path1,file1); % X has year in col 1, data in col 2
   yr=X(:,1); x=X(:,2);
   clear guy day X;
elseif strcmp(src,'dat');
   eval(['load ' pf1  ' -ascii;']);
   filepref = strtok(file1,'.');
   eval(['x = ' filepref '(:,2);']);
   eval(['yr = ' filepref '(:,1);']);
   clear filepref;
end;

%--- Floating series: store data, name  & year for later use
u1=x; yru1=yr; name1=file1;


%-- Reference Series -- specify type of file
ksrc=menu('Choose Type of Source File for Reference Series',...
   'Index in ITRDB-formatted .crn file',...
   'Ring width in single .rw file','Two-column time series mtx, year as col 1');
if ksrc==1;
   src = 'crn';
elseif ksrc==2;
   src = 'rw';
elseif ksrc==3; 
   src = 'dat';
else;
   error('ksrc must be 1, 2 or 3');
end


%--- Reference Series, get the file name
if strcmp(src,'crn');
   [file2,path2]=uigetfile('*.crn','Infile with chronology index');
elseif strcmp(src,'rw');
   [file2,path2]=uigetfile('*.rw','Infile with measured ring width series');
elseif strcmp(src,'dat');
   [file2,path2]=uigetfile('*.dat','Infile of 2-col time series matrix');
end;
pf2=[path2 file2]; % merge path and filename

%--- Reference Series, read the input time series;  
if strcmp(src,'crn');
   [x,s,yr]=crn2vec2(pf2); %x is index, s is sample size, yr is year
elseif strcmp(src,'rw');
   [X,guy,day,pf2]=rwread3(path2,file2); % X has year in col 1, data in col 2
   yr=X(:,1); x=X(:,2);
   clear guy day X;
elseif strcmp(src,'dat');
   eval(['load ' pf2  ' -ascii;']);
   filepref = strtok(file2,'.');
   eval(['x = ' filepref '(:,2);']);
   eval(['yr = ' filepref '(:,1);']);
   clear filepref;
end;

%--- Reference series: store data, name  & year for later use
u2=x; yru2=yr; name2=file2;

%--- Cleanup
clear x yr src ksrc

%--- COMPUTE SCALED FIRST DIFF SERIES
w1=rwchng(u1);  
w2=rwchng(u2);
yrw1=yru1;  yrw1(1)=[];
yrw2=yru2;  yrw2(1)=[];
nw1=length(w1);
nw2=length(w2);

%--- MAKE WINDOWS WITH TIME SERIES PLOTS OF FULL ORIGINAL SERIES
figure(1); % floating series
plot(yru1,u1); 
title(['Floating Time Series, ' name1]);
xlabel('Year');
grid;
zoom xon;

figure(2); % reference series
plot(yru2,u2); 
title(['Reference Time Series, ' name2]);
xlabel('Year');
grid;
zoom xon;

%--- MAKE WINDOWS WITH FULL LENGTH SCALED FIRST DIFFERENCES
figure(3); % floating series
plot(yrw1,w1); 
title([name1 ': Scaled First-Difference']);
xlabel('Year');
grid;
zoom xon;

figure(4); % reference series
plot(yrw2,w2); 
title([name2 ': Scaled First-Difference']);
xlabel('Year');
grid;
zoom xon;



%--- COMPUTE AND PLOT SLIDING STANDARD DEVIATION OF FLOATING SERIES
% use defwidth-year window

% Compute ending year or each period
iend=nw1:-5:dwin;
iend=(fliplr(iend))'; % to col vector, increasing
nper = length(iend); % number of standard devs
istart = iend-defwidth+1; % start year of period
Wtemp = repmat(NaN,defwidth,nper); % to hold sub-period data
for n = 1:nper;
   i1 = (istart(n):iend(n))';
   Wtemp(:,n)=w1(i1);
end;
s1 = (std(Wtemp))'; % standard deviations for sub-periods
yrs1 = (yrw1(iend));

figure(5);
plot(yrs1,s1);
grid;
title(['Moving Window of Standard Deviation, Scaled First-Difference of ' name1]);
xlabel(['Ending Year of ' int2str(dwin) '-Year Period']);
ylabel('Standard Deviation');
clear Wtemp i1 iend istart nper ;


%--- ANALYSIS LOOP

kwh1 = 1; % while control for level-1 menu
while kwh1==1; 
   
   %--- Set window width
   kmen1 = menu('Choose','Set window width', 'Quit');
   if kmen1 == 1;
      kmen4=menu('Choose',['Use window width of ' int2str(dwin)],'Input new width');
      if kmen4==1; % use current value of window width
         windwid=dwin;
      else;
         prompt={'Enter width'};
         def={int2str(dwin)};
         dlgTitle='Input for width of sliding window (yr)';
         lineNo=1;
         answer=inputdlg(prompt,dlgTitle,lineNo,def);
         windwid=str2num(answer{1});
         dwin=windwid; % re-set default window
      end;
      
      %--- Main analysis 
      
      % Pick test segment
      kmen1a=menu('Choose','Pick end year of test segment','Quit');
      if kmen1a==1; % pick end year of test segment
         kmen1a1=menu('Choose','Graphically--click on end year','In text window');
         if kmen1a1==1; 
            % bring plot of scaled first diffs, floating series, forward
            figure(3);
            xlimits=get(gca,'XLim');
            [t1,w1t1]=ginput(1);
            t1 = round(t1);
            if t1<xlimits(1);
               t1=xlimits(1);
            else;
            end;
            if t1>xlimits(2);
               t1=xlimits(2);
            else;
            end;
            yrend=t1;
            clear t1 xlimits;
         else;
            prompt={'Enter end "year"'};
            def={int2str(max(yru1))}; % default is last year of floating series
            dlgTitle='Input last year of test segment (as numbered in floating series)';
            lineNo=1;
            answer=inputdlg(prompt,dlgTitle,lineNo,def);
            yrend=str2num(answer{1});
         end; % if kmen1a1==1
         
         %--- Check that selected segment valid
         if yrend<yrw1+windwid-1;
            error([int2str(windwid) '-yr period ending in ' int2str(yrend) ' impossible']);
         end;
         yrstart = yrend-windwid+1;
         str1 = sprintf('%5.0f-%5.0f',yrstart,yrend);
         
                 
         kwh2=1; % while control for different test segment
         while kwh2==1;
            % Pull segment of floating series
            L1=yrw1>=yrstart & yrw1<=yrend;
            v1=w1(L1); yrv1=yrw1(L1);
            
            % Make matrix of reference series
            nH=nw2-windwid+1; % number of segments
            Z= repmat(NaN,windwid,nH);
            j1=(0:(windwid-1))';  % column vector
            J1 = repmat(j1,1,nH); 
            j2 = 1:nH; % row vector
            J2 = repmat(j2,windwid,1);
            J3 = J1 + J2; % each col has row indices into Z
            for n = 1:nH;
               j3=J3(:,n);
               Z(:,n)=w2(j3);
            end;
            yrZ = ((yrw2(1)+windwid-1):yrw2(nw2))'; % col vect of ending years of segments
            
            % Compute correlations
            r = corrone(v1,Z);
            
            % Correlations plot
            figure(6);
            plot(yrZ,r);
            grid;
            zoom xon; 
            title(['Correlation of segment ' str1 ' of ' name1]);
            xlabel(['End Year of Segment in ' name2]);
            ylabel('Correlation coefficient');
            
            % Box plot
            figure(7);
            boxplot(r,1,'+',1,1)
            title(['Boxplot of Correlation Coefficients']);
            
            % Top 10 alignments
            figure(8);
            [rsort,isort]=sort(r);
            n10=min([length(isort) 10]);
            rsort=flipud(rsort');
            isort=flipud(isort');
            S1 = ['Correlations with ' name2 ' segments ending in year:'];
            for n = 1:n10;
               str2=sprintf('%5.0f   %6.2f',yrZ(isort(n)), rsort(n));
               S1=strvcat(S1,str2);
            end;
            rwow=rsort(1);
            yrwow=yrZ(isort(1));
            yrbingo=yrZ(isort(1));
            S1=cellstr(S1);
            axes;
            xlims = get(gca,'XLim');
            ylims=get(gca,'YLim');
            xrange=abs(diff(xlims));
            yrange=abs(diff(ylims));
            text(xlims(1)+xrange/10,ylims(2)-yrange/10,S1,...
               'VerticalAlignment','top');
            set(gca,'XTick',[]);
            title([name1 ' Segment ' str1 ' Correlation Summary']);
            
            % Time series plots for period of best match
            figure(9);
            subplot(2,1,1);
            % recall that yrend is the last year in floating segment a
            % recall that yrbingo is the best match year in reference series
            Lu1 = yru1>=yrstart & yru1<=yrend;
            Lw1 = yrw1>=yrstart & yrw1<=yrend;
            yrborn = yrbingo-windwid+1;
            Lu2 = yru2>=yrborn & yru2<=yrbingo;
            Lw2 = yrw2>=yrborn & yrw2<=yrbingo;
            
            % Raw-data plots
            rr=corrcoef([u1(Lu1) u2(Lu2)]);
            rr=rr(1,2);
            strcor1 = sprintf('r = %6.2f',rr);
            plot(yru2(Lu2),zscore(u1(Lu1)),yru2(Lu2),zscore(u2(Lu2)));
            grid;
            title(['Raw Data For Best Match Period: ' strcor1]);
            xlabel('Year');
            ylabel('Z-scores');
            legend(name1,name2);
            zoom xon;
            
            % Scaled First-diff plots
            rr=corrcoef([w1(Lw1) w2(Lw2)]);
            rr=rr(1,2);
            strcor2 = sprintf('r = %6.2f',rr);
            subplot(2,1,2);
            plot(yrw2(Lw2),zscore(w1(Lw1)),yrw2(Lw2),zscore(w2(Lw2)));
            grid;
            title(['Scaled First Differences for Best Match Period: ' strcor2]);
            xlabel('Year');
            ylabel('Z-scores');
            legend(name1,name2);
            zoom xon;
            
            % Bring correlation ts to front
            figure(6);
            
            kmen1a2=menu('Choose','Satisfied','Different test segment?','Quit');
            if kmen1a2==1;
               kmen7=menu('Choose','Bootstrap Analysis','Skip Bootstrapping');
               if kmen7==1; % do bootstrap
                   prompt={'Enter number of samples:'};
                   def={'1000'};
                   dlgTitle='Input desired number of bootstrap samples';
                   lineNo=1;
                   answer=inputdlg(prompt,dlgTitle,lineNo,def);
                   nboot=str2double(answer{1});
                   rboot = repmat(NaN,nboot,1); % to hold bootstrapped r
                   % Generate nboot bootstrap samples of the floating segment
                   % Recall that v1 is the floating segment, with years yrv1
                   [bstat,bsam]=bootstrp(nboot,'mean',v1); % each col of bsam has indices to v1
                   for n = 1:nboot;
                      iboot = bsam(:,n);
                      xboot=v1(iboot);
                      rtemp =max(corrone(xboot,Z)); % highest correl with any segment of ref series
                      rboot(n)=rtemp;
                   end;
                   figure(10);
                   boxplot(rboot);
                   ytemp=get(gca,'YLim');
                   ytemp(2)=max([max(rboot)  min([rwow+.1 1.0])]);
                   set(gca,'YLim',ytemp);
                   text(1,rwow,'o');
                   text(1.1,rwow,'<--- Observed');
                   title(['Maximum r for ' int2str(nboot) ' bootstrapped samples']);
                else; % skip bootstrapping
                end;
             elseif kmen1a2==2; % pick diffent segment
                close(6); close(7); close(8); close(9);
                kwh2=0;
             elseif kmen1a2==3; % Quit 
                kwh2=0;
             end; % kmen1a2==1
          end; % kwh2==1
          
       else;
       end; % if kmen1a==1
       
       
    else;
       kmen1b=menu('Choose','Close windows before quitting','Keep windows open');
       if kmen1b ==1
          close all;
       else;
          % no action- keep windows open
       end;
       kwh1=0;
    end; % if kmen1==1
 end; % while kwh1==1
 