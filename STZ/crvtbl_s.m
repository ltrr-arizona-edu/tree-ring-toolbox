function T=crvtbl_s(S,cmask,nms,yrs,P);
% crvtbl_s:  Summary table of curve fits from crvfit.m
%
%
% See /howto/wordtbl1.txt for followup to put in ms word table

% Excludes from the table (1) any cores not yet fit, as indicated by col 1 entry
% of zero in S, and (2) any masked core, as indicated by 0 in cmask

% Logical pointer to cores that have been fit
L1=S(:,1)~=0;

% Logical pointer to cores not masked
if ~islogical(cmask); error('cmask not logical'); end
L2=cmask;


% Logical pointer to cores that have been fit and are not masked
L3=L1 & L2;


% Cull nms, yrs, S
S=S(L3,:);
yrs=yrs(L3,:);
nms=nms(L3,:);
arord = P(L3,4); % AR order used to create residual core index
varrat = P(L3,5); % ratio of variance of ar residuals to original index

ncores = size(nms,1); % number of selected cores



%--- SHORT VERSION OF TABLE -- FOR ITRDB

% Sequence number
s1 = (1:ncores)';
s1=num2str(s1,'%2d ');

% Core id
s2 = nms;

% First and last year of fit period
s3=S(:,[10 11]);
s3=num2str(s3,'%5.0f %5.0f');

% Flag (*) if fit period differs from full available period of ringwidth in rwl file
Lgo=S(:,10)~= yrs(:,1);
Lsp=S(:,11)~=yrs(:,2);
Lflag = Lgo | Lsp;
sflag = repmat(blanks(2),ncores,1);
if any(Lflag);
   sflag(Lflag,1)='*';
end;
s3=[s3 sflag repmat(blanks(2),ncores,1)];



% ----Curve type

stemp=repmat('XX',ncores,1);
s4=S(:,2);

% negative exponential
L1 = s4==1;
L1sum=sum(L1);
if L1sum>0;
   stemp(L1,:)=repmat('NE',L1sum,1);
end;

% straight line, least squares fit
L2 = s4==2;
L2sum=sum(L2);
if L2sum>0;
   stemp(L2,:)=repmat('SL',L2sum,1);
end;

% horizontal line thru mean
L4 = s4==4;
L4sum=sum(L4);
if L4sum>0;
   stemp(L4,:)=repmat('HL',L4sum,1);
end;

% spline
L9 = s4==9;
L9sum=sum(L9);
if L9sum>0;
   stemp(L9,:)=repmat('CS',L9sum,1);
end;

s4=stemp;  % curve types

%--- If spline, period (yr) at which freq response amplitude is 0.5 (50% wavelength)
stemp = repmat(blanks(9),ncores,1);
if L9sum>0;
   jj = find(L9);
   for j = 1:L9sum;
      % Compute period of 50% freq response
      pspline = S(jj(j),3);
      p50 = spline50(pspline);
      if round(p50)>=1000;
         ss = sprintf('(%4.0f yr)',p50);
      elseif round(p50)>100;
         ss = sprintf('(%3.0f yr)',p50);
      elseif round(p50) >10;
         ss = sprintf('(%2.0f yr)',p50);
      elseif round(p50) >1;
         ss = sprintf('(%1.0f yr)',p50);
      else; 
         error('spline with 50% period shorter than 1 yr');
      end;
      nss = length(ss);
      stemp(jj(j),1:nss) = ss;
   end;
end;
% Lop off any trailing columns of stemp that are all blanks
LL =  all(isspace(stemp));
if any(LL);
   stemp(:,LL)=[];
end;
% Add two blank spaces
s5 = [stemp  repmat(blanks(2),ncores,1)];

%-- AUTOREGRESSIVE MODEL ORDER AND FRACTIONAL PCT VARIANCE DUE TO PERSISTENCE
s6 = num2str(arord,'   AR(%1d)  ');
s7 = num2str(1-varrat,'  %3.2f');

shd=...
   'Sample size in .crn file is the number of trees ';
shd=char(shd,...
   'rather than the number of cores.                ');
shd=char(shd,...
   'Core indices computed as the ratio of ringwidth ');
shd=char(shd,...
   'to fitted trend line.  Core indices averaged to ');
shd=char(shd,...
   'to form tree indices.  Tree indices averaged to ');
shd=char(shd,...
   'to form site index.  Bi-weight mean used if 6 or');
shd=char(shd,...
   'more trees in a given year.  Otherwise chron    ');
shd=char(shd,...
   'computed by arithmetic mean. Residual chron     ');
shd=char(shd,...
   'computed likewise, except from autoregressive   ');
shd=char(shd,...
   'residuals of core indices.  All ringwidth series');
shd=char(shd,...
   'used in chronology are listed in table below,   ');
shd=char(shd,...
   'with sufficient information to allow independent');
shd=char(shd,...
   'researcher to duplicate construction of site    ');
shd=char(shd,...
   'chronology from the ringwidths in the .rwl file ');


stit='N Core    Fit Period    Curve fit    Persistence'
sall=[s1 s2 s3 s4 s5 s6 s7];
sall=[ stit; sall];

sall=char(sall,blanks(3));
sall=char(sall,'Fit Period = period of ring width used in curve fitting');
sall=char(sall,'  * = differs from full available period of ring width data');
sall=char(sall,'Curve Fit = type of model used to detrend ring width');
sall=char(sall,'  NE=modified negative exponential');
sall=char(sall,'  SL=least-square-fit straight line');
sall=char(sall,'  HL=horizontal line at mean of ring width for fit period');
sall=char(sall,'  CS=cubic smoothing spline');
sall=char(sall,'     (Wavelength of 50% frequency reponse in parentheses)');
sall=char(sall,'Persistence = autoregressive prewhitening information');
sall=char(sall,'  AR(q) = order q of autoregressive model fit to core index');
sall=char(sall,'  nn  = decimal proportion of core-index variance due to ');
sall=char(sall,'     persistence; computed as 1-R, where R is the ratio of ');
sall=char(sall,'     variance of residual index to variance of standard index');


shd=char(shd,blanks(3));
sall=char(shd,sall);

[ms,ns]=size(sall);
%-- Write itrdb info text file

[file5,path5]=uiputfile('*.txt','Output .txt file for export to NGDC ');
pf5=[path5 file5];
fid5=fopen(pf5,'w');
for n = 1:ms;
   s = sall(n,:);
   fprintf(fid5,'%s\r\n',s');
end;
fclose(fid5);


%Allocate
T=blanks(100);
T=T(ones(ncores,1),:);


% Convert information on selected cores to string
for n = 1: ncores
   c1type=S(n,2); % curve-type
   nmcore=nms(n,:);
   str1=sprintf('%s\t',nmcore);
   seqno= deblank(int2str(S(n,1))); % sequence number
   str2=sprintf('(%s)\t',seqno);  % sequenc number
   str3=sprintf('%4.0f\t%4.0f\t ', S(n,10),S(n,11)); % start and end year for fit
   str4=sprintf('%s\t',fittype(c1type));
   if c1type==1; % if neg exponential fit
      str5=sprintf('%9.5g\t%0.5g\t%7.2E\t',S(n,3),S(n,4),S(n,5));
   elseif c1type==2; % straight line, no restriction on sign of slope
      str5=sprintf('%9.5g\t%9.5g\t%5s\t',S(n,3),S(n,4),'N/A');   
   elseif c1type==4; % horiz line through mean
      str5=sprintf('%9.5g\t%5s\t%5s\t',S(n,3),'N/A','N/A');
   elseif c1type==9; % spline
      str5=sprintf('%9.5g\t%0.5g\t%8.6g\t',S(n,3),S(n,4),S(n,5));
   end
   str6 = sprintf('%4.0f\t%5.2f',arord(n),varrat(n));
      
   strall=[str1 str2 str3 str4 str5 str6];
   nchars=length(strall);
      
   T(n,1:nchars)=strall;
end


maxlen=0;
for n=1:ncores;
   tt=deblank(T(n,:));
   maxlen=max(maxlen,length(tt));
   sprintf('%s\n',tt);
end

T=T(:,1:maxlen);


% Ascii file
[file1,path1]=uiputfile('*.txt','Ascii table');
pf1=[path1 file1];
fid1=fopen(pf1,'w');

% Header
% Note: set microsoft word font to new courier 8 to show the table lined up right
%
wordinf1 = 'To format file in ms word: (1) open the .txt file as a word doc';
wordinf2 = '(2) change font size to 8; (3) select  the last';
wordinf3 = 'header line and table body and click the table grid';
wordinf4 = '(4)use table formatting to get to look right';
wordinf5= '(5) make sure to make the superscripts super';
wordinf6 = '  ';

hdr1='Period3   Curve4      Parameters5   AR Model6';;

hd1='Core1';  hd2='Seq2';  hd3='Start'; hd4='End'; 
hd5='Type';  hd6='k'; hd7='a';  hd8='b';  hd9='p'; hd10='Ratio';

fprintf(fid1,'%s\n',wordinf1);
fprintf(fid1,'%s\n',wordinf2);
fprintf(fid1,'%s\n',wordinf3);
fprintf(fid1,'%s\n',wordinf4);
fprintf(fid1,'%s\n',wordinf5);

fprintf(fid1,'%s\n',hdr1);
fprintf(fid1,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',hd1,hd2,hd3,hd4,...
   hd5,hd6,hd7,hd8,hd9,hd10);

for n=1:ncores
   fprintf(fid1,'%s\n',T(n,:));
end

% Footer
foot1='1Core: Identifies site, tree and core';
fprintf(fid1,'\n\n%s\n',foot1);
fprintf(fid1,'%s\n','2seq.: Sequence # of core in .rwl ring-width storage file');
fprintf(fid1,'%s\n','3Period: Start and end year of period used to fit detrending curve');
fprintf(fid1,'%s\n','4Type: Curve type used in detrending:');
  fprintf(fid1,'  %s\n','NE - modified negative exponential');
  fprintf(fid1,'  %s\n','CS - cubic smoothing spline');
  fprintf(fid1,'  %s\n','HL - horizontal line through sample mean');
  
  fprintf(fid1,'%s\n','5Parameters: parameters of detrending curve; definition depends on curve type');
fprintf(fid1,'  %s\n','NE: equation y=k + a*exp(-b*(t-start))');
fprintf(fid1,'      %s\n','where y is smooth curve value in 100ths of mm');
fprintf(fid1,'      %s\n','  k, a, and b are parameters fit by least squares');
fprintf(fid1,'      %s\n','  t is the year, and');
fprintf(fid1,'      %s\n','  start is the first year of the fit period (column 3 this table)');


fprintf(fid1,'  %s\n','SL: equation y=a + b*(t-tstart)');
fprintf(fid1,'      %s\n','where y is smooth curve value in 100ths of mm');
fprintf(fid1,'      %s\n','  a and b are parameters fit by least squares');
fprintf(fid1,'      %s\n','  t is the year, and');
fprintf(fid1,'      %s\n','  start is the first year of the fit period (column 3 this table)');

fprintf(fid1,'  %s\n','CS: cubic spline parameter and properties');
fprintf(fid1,'      %s\n','k = spline parameter');
fprintf(fid1,'      %s\n','a = period (years) at which amplitude of frequency response of spline equals b');
fprintf(fid1,'      %s\n','b = amplitude of frequency response at period "a" years');

fprintf(fid1,'  %s\n','HL: Horizontal line information');
fprintf(fid1,'      %s\n','k = mean ring width for the period listed in columns 3 and 4');

fprintf(fid1,'%s\n','6Information on autoregressive model used to generate residual core index');
fprintf(fid1,'  %s\n','p = order of the model');
fprintf(fid1,'  %s\n','Ratio = ratio of variance of residual index to variance of standard index');
fprintf(fid1,'      %s\n','(high persistence in standard index corresponds to a low ratio');


fclose(fid1);

%SUBFUNCTION
function sc=fittype(c);
% returns a string identifier for a curve-fit type number
if c==0;  % No second curve-type
   sc='No';
elseif c==1; % negative exp
   sc='NE';
elseif c==2; % straight line, fit by least squares
   sc='SL';
elseif c==4; % horizontal line thru mean
   sc='HL';
elseif c==9; % spline
   sc='CS';
else
end
   
   
   