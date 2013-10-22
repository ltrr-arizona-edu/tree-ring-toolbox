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
   
   
   