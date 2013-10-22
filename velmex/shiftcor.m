function shiftcor(jpick)
%
% YET:  1) use ratio for high-pass if ring-width, 2) code for other types of file input besides choice 1
% shiftcor: check crossdating of undated time series against a master by high-pass filtering and correlation
% shiftcor(jpick);
% Last revised 2006-9-27
%
% Check crossdating of undated time series against a master by high-pass filtering and correlation. 
% The user selects a part of the undated time series to be slid year-by-year along the master time series.
% A correlation coefficient is computed at each iteration. The best  match is defined as that with the 
% highest correlation. An optional bootstrap analysis can be used to check how likely this maximum correlation
% is by chance alone.<P>
%

%
%*** INPUT
%
% Prompts for names of files containing the master series. 
% Prompt for width of window for sliding segment of undated series. 
% Prompt for ending year of sliding segment.
% Prompt for bootstrap or not
% Prompt for number of bootstrap samles
%
%*** OUTPUT
%
% Screen plots summarizing results.  Key windows are graph of correlation coefficient 
% against year in master series, and boxplots showing distribution of correlations.  
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
% shiftcor.m was written to help in crossdating remnant-wood ring-width series against 
% very long tree-ring chronologies from a nearby site. t 
%
% You get to select the segment of the undated series to slide interactively, with 
% guidance from time series plots of the undated series and variance of the undated 
% series.
%
% To emphasize similarity in year-to-year changes as opposed to trends, series
% are high-pass filtered at the start with a 20-year splines.   
%
% You select the segment of the undated series to be slid along the master 
% series.  This series is slid several years along the master series and
% the correlation is computed at each step.  The correlation is plotted as a function 
% of ending year (in the master series) of the alignment.  A boxplot of the 
% correlations gives some idea of variability of the correlation coefficient for 
% various alignments.  
%
% A plot of variance of the undated series in a sliding window is given to help in 
% choosing the segment of the undated series to test.  The notion is that a period of 
% high variability is a good candidate for getting a definitive match.
%
% Bootstrapped significance level of maximum correlation is available on option.  
% This option takes a while to run.  The bootstrapping can be done once you 
% have found the maximum sliding correlation of the undated segment with the 
% master series.  A user-specified number of bootstrap samples are drawn from the 
% undated segment (e.g., 1000 samples). Each bootstrap sample is slid along the 
% master series and the maximum correlation is recorded.  The distribution of 
% maximum bootstrap correlations can then be compared to the maximum correlation 
% found for the original (non-bootstrapped) data. 
%
% Notation
% variable 1: undated series
% variable 2: master series
% u full length original series
% v windowed original series
% w windowed scaled high-pass series
% s standard-deviation of rw change, undated series

close all;

% Hard Code
wind_width=50; % default window width
rthresh = 2/sqrt(wind_width); % critical threshold
wind_shift=5; % initialize window width
spline_period=20; % wavelength spline has frequency response 0.
p_spline = splinep(spline_period,0.50); % spline parameter
str_spline1 = [num2str(spline_period) ' spline'];
str_spline2 = ['High Pass (' num2str(spline_period) ' spline)'];

if nargin==0;
    jpick=[];
end

%--- PROMPT FOR TYPES OF FILES WITH INPUT

filemodes={'Undated series from a .rwl file, master from a crn',...
    'Undated series from a .rwl file, master from specified cores from same .rwl file',...
    'Undated series from a .crn file, master from a different .crn file',...
    'Undated series from a .crn file, master from a two-column ascii file'};
filesuffs1 ={'.rwl','.rwl','.crn','.crn'};
filesuffs2={'.crn','.rwl','.crn','.txt'};
kinput = menu('Choose form of input file(s)',filemodes);

if kinput~=1;
    error('Not yet coded for other choices of filemodes')
end

%--- GET FILENAMES OF FILES WITH UNDATED AND MASTER SERIES

if kinput==1; % undated from rwl, master from crn
    [file1,path1]=uigetfile('*.rwl','rwl infile with undated series');
    pf1=[path1 file1];
    [file2,path2]=uigetfile('*.crn','crn Infile with master series');
    pf2=[path2 file2];
    
elseif kinput==2; % undated from rwl, master from specified subset of cores in same rwl
    [file1,path1]=uigetfile('*.rwl','rwl infile with undated and master');
    pf1=[path1 file1];
    pf2=pf1;
    file2=file1;
elseif kinput==3; % undated from a crn, master from a different crn
    [file1,path1]=uigetfile('*.crn','crn Infile with undated series');
    [file2,path2]=uigetfile('*.crn','crn Infile with master series');
    pf1=[path1 file1];
     pf2=[path2 file2];
elseif kinput==4; % undated from a crn, master from a two-column ascii
    [file1,path1]=uigetfile('*.crn','crn Infile with undated series');
    pf1=[path1 file1];
    [file2,path2]=uigetfile('*.txt','txt infile with master series');
    pf2=[path2 file2];
end;


%--- INPUT DIALOG TO SPECIFY SPINE WAVELENGTH, WINDOW WIDTH

    prompt={'Enter the spline wavelength:','Enter the sliding window width:','Enter the offset:'};
    name='Program settings';
    numlines=1;
    defaultanswer={num2str(spline_period),num2str(wind_width),num2str(wind_shift)};
    answer=inputdlg(prompt,name,numlines,defaultanswer);
    spline_period=str2num(answer{1});
    wind_width=str2num(answer{2});
    wind_shift = str2num(answer{3});
    p_spline = splinep(spline_period,0.50); % spline parameter



%--- INPUT THE  TIME SERIES

if kinput==1; % undated from rwl, master from crn
    pf1a=rwlinp(pf1);
    eval(['load ' pf1a ';']); % X, nms, yrs are key
    if size(jpick,1)==1;
        jpick=jpick';
    end;
    if isempty(jpick);
        jpick=(1:size(nms,1))';
    end
    Xall=X;
    nmsall=nms;
    yrsall=yrs;
    YRS = yrsall(jpick,1:2); % start and end year of undated series
    
    % Put undate series in time series matrix
    [X,yrX,nms]=sov2tsm3(X,yrs,nms,jpick,[]); % x is tsm of rw series, yrx is year vector; nms is col-cell of strings
    [mX,nX]=size(X);
    ntest = nX; % this many ring width series to be tested
    
    % Build selection menu and select undated seires
    pL = repmat('(',nX,1);
    pR = repmat(')',nX,1);
    S=[num2str(jpick) repmat(' ',nX,1) char(nms)  repmat(' ',nX,1) pL num2str(YRS(:,1)) repmat('-',nX,1) num2str(YRS(:,1))   pR];
    kpick = menu('Choose',cellstr(S)); % choose undated series
    x =X(:,kpick);
    [x,yr]=trimnan(x,yrX);
    file1=[file1 '--' nms{kpick}];
    u1=x; yru1=yr; name1=file1; % store undated series 
    
    % Get master
    [x,s,yr]=crn2vec2(pf2); %x is index, s is sample size, yr is year
    u2=x; yru2=yr; name2=file2;
    
elseif kinput==2; % undated from rwl, master from specified subset of cores in same rwl
    
elseif kinput==3; % undated from a crn, master from a different crn
    
elseif kinput==4; % undated from a crn, master from a two-column ascii
    
end;

% 
% 
% 
% %--- Undated Series, read the input time series;
% if strcmp(src,'crn');
%    [x,s,yr]=crn2vec2(pf1); %x is index, s is sample size, yr is year
% elseif strcmp(src,'rwl');
%   
% elseif strcmp(src,'txt');
%    eval(['load ' pf1  ' -ascii;']);
%    filepref = strtok(file1,'.');
%    eval(['x = ' filepref '(:,2);']);
%    eval(['yr = ' filepref '(:,1);']);
%    clear filepref;
% end;
% 
% %--- Undated series: store data, name  & year for later use
% % Note that u1 and x can be a vector or matrix, depending on whether src is crn, txt or rwl
% u1=x; yru1=yr; name1=file1;
% 
% 
% 
% %--- Reference Series, read the input time series;  
% if strcmp(src,'crn');
%    [x,s,yr]=crn2vec2(pf2); %x is index, s is sample size, yr is year
% elseif strcmp(src,'rwl');
%     pf2a=rwlinp(pf2);
%     eval(['load ' pf2a ';']); % X, nms, yrs are key
%     jpick=[];
%     tends=[];
%     [x,yr,nms]=sov2tsm3(X,yrs,nms,jpick,tends); % x is tsm of rw series, yrx is year vector; nms is col-cell of strings
%     kpick = menu('Choose ref series',nms);
%     x =x(:,kpick);
%     [x,yr]=trimnan(x,yr);
%     file2=[file2 '--' nms{kpick}];
% elseif strcmp(src,'txt');
%    eval(['load ' pf2  ' -ascii;']);
%    filepref = strtok(file2,'.');
%    eval(['x = ' filepref '(:,2);']);
%    eval(['yr = ' filepref '(:,1);']);
%    clear filepref;
% end;
% 
% %--- Reference series: store data, name  & year for later use
% u2=x; yru2=yr; name2=file2;

%--- Cleanup
clear x yr src ksrc



%--- Check whether "test" series starts with year "1".  If so, consider this series to be test anywhere
% relative to the master series.  If test series starts with a real year (other than 1), consider the
% test series tentatively dated, and DO NOT allow picking of a test segment that would not be jointly fully
% covered by both the test and master series.
if yru1==1;
    true_float=1;
    startOk=-inf;
    endoOK=inf;
else;
    true_float=0;
    startOK = max([yru1(1) yru2(1)]); % the slab cannot start earlier than this
    endOK = min([yru1(end) yru2(end)]); % the slab cannot end later than this
end



%--- COMPUTE  spline, and ratio or departures in std units

p_spline = splinep(spline_period,0.50);
[s1,p] = csaps(yru1,u1,p_spline,yru1);
[s2,p] = csaps(yru2,u2,p_spline,yru2);

% Ratio or difference depending on data source
if kinput==1; % Undated series from a .rwl file, master from a crn
    if any(s1<=0);
        error('Spline for undated goes zero or negative, ratio index explodes');
    end
    w1 = u1 ./ s1;
    w2 = u2 - s2;
elseif kinput==2; % Undated series from a .rwl file, master from specified cores from same .rwl file
    if any(s1<=0) | any(s2<=0);
        error('Spline for undated or master goes zero or negative, ratio index explodes');
    end
    w1 = u1 ./ s1;
    w2 = u2 ./ s2;
elseif kinput==3; % Undated series from a .crn file, master from a different .crn file
     w1 = u1 - s1;
     w2 = u2 - s2;
else; Undated series from a .crn file, master from a two-column ascii file
     w1 = u1 - s1;
     w2 = u2 - s2;
end

w1mn =nanmean(w1);
w2mn = nanmean(w2);
w1std = nanstd(w1);
w2std = nanstd(w2);
w1 = (w1-w1mn)/w1std;
w1 = (w1-w1mn)/w1std;

yrw1=yru1;
yrw2=yru2;
nw1=length(w1);
nw2=length(w2);

%--- MAKE WINDOWS WITH TIME SERIES PLOTS OF FULL ORIGINAL SERIES

figure(1); % test series
hp1=plot(yru1,u1,yru1,s1); 
set(hp1(2),'Color',[1 0 0],'LineWidth',2);
legend('Original',str_spline1);
title(['Undated Time Series, ' name1]);
xlabel('Year');
grid;
zoom xon;

figure(2); % master series
hp1=plot(yru2,u2,yru2,s2); 
set(hp1(2),'Color',[1 0 0],'LineWidth',2);
legend('Original',str_spline1);
title(['Reference Time Series, ' name2]);
xlabel('Year');
grid;
zoom xon;

%--- MAKE WINDOWS WITH FULL LENGTH HIGH PASS
figure(3); % test series
plot(yrw1,w1); 
title([name1 ': HIGH PASS']);
xlabel('Year');
grid;
zoom xon;

figure(4); % master series
plot(yrw2,w2); 
title([name2 ': HIGH PASS']);
xlabel('Year');
grid;
zoom xon;


%--- COMPUTE AND PLOT SLIDING STANDARD DEVIATION OF FLOATING SERIES
% use defwidth-year window

% Compute ending year or each period
iend=nw1:-5:wind_width;
iend=(fliplr(iend))'; % to col vector, increasing
nper = length(iend); % number of standard devs
istart = iend-wind_width+1; % start year of period
Wtemp = repmat(NaN,wind_width,nper); % to hold sub-period data
for n = 1:nper;
   i1 = (istart(n):iend(n))';
   Wtemp(:,n)=w1(i1);
end;
s1 = (nanstd(Wtemp))'; % standard deviations for sub-periods
yrs1 = (yrw1(iend));

figure(5);
plot(yrs1,s1);
grid;
title(['Moving Window of Standard Deviation, high-pass of ' name1]);
xlabel(['Ending Year of ' int2str(wind_width) '-Year Period']);
ylabel('Standard Deviation');
clear Wtemp i1 iend istart nper ;


%--- ANALYSIS LOOP

kwh1 = 1; % while control for level-1 menu
while kwh1==1; 
   
   %--- Set window width
   kmen1 = menu('Choose','Continue', 'Quit');
   if kmen1 == 1;
%       kmen4=menu('Choose',['Use window width of ' int2str(wind_width)],'Input new width');
%       if kmen4==1; % use current value of window width
%          windwid=wind_width;
%       else;
%          prompt={'Enter width'};
%          def={int2str(wind_width)};
%          dlgTitle='Input for width of sliding window (yr)';
%          lineNo=1;
%          answer=inputdlg(prompt,dlgTitle,lineNo,def);
%          windwid=str2num(answer{1});
%          wind_width=windwid; % re-set default window
%       end;
      endOK1=startOK+wind_width-1; % earliest slab cannot end before this
      endOK2 = endOK; % latest slab cannot end aftr this
      str_endOK = ['End year of slab must be in range ' num2str(endOK1) '-' num2str(endOK2)];
      
      %--- Main analysis 
      
      % Pick test segment
      kmen1a=menu(['Choose slab: ' str_endOK],'Pick end year of test slab of ''test'' series ','Quit');
      if kmen1a==1; % pick end year of test segment
         kmen1a1=menu('Choose','Graphically--click on end year','In text window');
         if kmen1a1==1; 
            % bring plot of high-pass filtered series, test series, forward
            figure(3);
            xlimits=get(gca,'XLim');
            kwh3=1;
            while kwh3==1
                [t1,w1t1]=ginput(1);
                t1 = round(t1);
                if t1<endOK1 | t1>endOK2
                    uiwait(msgbox(str_endOK,'Pick again','modal'));
                else
                    kwh3=0;
                end
            end
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
             kwh3=1;
             while kwh3==1;
                 prompt={['Enter end "year" (' str_endOK ')']};
                 def={int2str(max(yru1))}; % default is last year of test series
                 dlgTitle='Input last year of test segment (as numbered in test series)';
                 lineNo=1;
                 answer=inputdlg(prompt,dlgTitle,lineNo,def);
                 yrend=str2num(answer{1});
                 if yrend<endOK1 | yrend>endOK2
                    uiwait(msgbox(str_endOK,'Pick again','modal'));
                else
                    kwh3=0;
                end
             end;


         end; % if kmen1a1==1
         
         %--- Check that selected segment valid
         if yrend<yrw1+wind_width-1;
            error([int2str(wind_width) '-yr period ending in ' int2str(yrend) ' impossible']);
         end;
         yrstart = yrend-wind_width+1;
         str1 = sprintf('%5.0f-%5.0f',yrstart,yrend);
         
                 
         kwh2=1; % while control for different test segment
         while kwh2==1;
            % Pull segment of test series
            L1=yrw1>=yrstart & yrw1<=yrend;
            v1=w1(L1); yrv1=yrw1(L1);
            
            % Make matrix of master series
            nH=nw2-wind_width+1; % number of segments
            Z= repmat(NaN,wind_width,nH);
            j1=(0:(wind_width-1))';  % column vector
            J1 = repmat(j1,1,nH); 
            j2 = 1:nH; % row vector
            J2 = repmat(j2,wind_width,1);
            J3 = J1 + J2; % each col has row indices into Z
            for n = 1:nH;
               j3=J3(:,n);
               Z(:,n)=w2(j3);
            end;
            yrZ = ((yrw2(1)+wind_width-1):yrw2(nw2))'; % col vect of ending years of segments
            
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
            % recall that yrend is the last year in test segment a
            % recall that yrbingo is the best match year in master series
            Lu1 = yru1>=yrstart & yru1<=yrend;
            Lw1 = yrw1>=yrstart & yrw1<=yrend;
            yrborn = yrbingo-wind_width+1;
            Lu2 = yru2>=yrborn & yru2<=yrbingo;
            Lw2 = yrw2>=yrborn & yrw2<=yrbingo;
            
            % Raw-data plots
            rr=corrcoef([u1(Lu1) u2(Lu2)]);
            rr=rr(1,2);
            strcor1 = sprintf('r = %6.2f',rr);
            plot(yru2(Lu2),zscore(u1(Lu1)),yru2(Lu2),zscore(u2(Lu2)));
            grid;
            title(['Z-score original series for period of best match: ' strcor1]);
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
            title(['Z-score high-pass series for period of best match: ' strcor2]);
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
                   % Generate nboot bootstrap samples of the test segment
                   % Recall that v1 is the test segment, with years yrv1
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
 