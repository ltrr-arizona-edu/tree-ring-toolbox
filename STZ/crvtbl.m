function crvtbl
% crvtbl:  Summary table of chronology development for ITRDB submission and more
% crvtbl;
% Last Revised 12-19-00
%
% Two ascii tables summarizing chronology development by the crvfit, corei,... 
% sequence of Matlab Functions. One of the .txt files produced is in a form 
% that can be copied/pasted from wordpad into the online ITRDB submission form. 
% The other, more detailed .txt file, can be inserted into MS Word. 
%
%*** INPUT 
%
% Prompted to click on template text description file on strategy for detrending
% Prompted to click on input .mat file with stored indices, ringwidths, etc
% Prompted for various types of site information (e.g., long, lat, elev).
% Prompted to click on input and output filenames
%
%*** OUTPUT
%
% Two ascii files of data summarizing transformation, curve fits, 
% prewhitening, chronology make-up. By default, cfita.txt and cfitb.txt are
% the filenames.  cfita.txt is for export to NGDC. cfitb.txt is for reformatting
% as MS word table.
%
%*** REFERENCES -- none
%
%*** UW FUNCTIONS CALLED 
%
% spline50.m
%
%*** TOOLBOXES NEEDED -- spline
%
%
%*** NOTES
%
% See /howto/wordtbl1.txt for followup to put in ms word table

% String and numeric versions of common powers
pstr={'  -2','  -1', '-2/3', '-1/2','-1/3',  '   0', ' 1/3', ' 1/2', ' 2/3', '   1', '   2', '   3'};
pnum=[   -2     -1   -2/3     -1/2   -1/3        0     1/3     1/2     2/3       1       2       3];

%--- Prompt for template of file-fit info
ktext=menu('Choose strategy used for inital automatic detrending',...
    'Ratio index; negative slope combo (non-inc spline, HL, SL+), followed by manual over-ride',...
    'Ratio index, Specific-wavelength spline (e.g., 200 yr)',...
    'Difference Index; negative slope combo ...');
    
if ktext==1;
    C = textread('c:\mlb\stz\cfitv1.txt','%s','delimiter','\n','whitespace','');
    C =char(C);
    C=C(1:33,:);
elseif ktext==2; % specific-wavelength spline
    C = textread('c:\mlb\stz\cfitv3.txt','%s','delimiter','\n','whitespace','');
    C =char(C);
    C=C(1:34,:);
    prompt={'Enter wavelength (yr):'};
    def={'200'};
    dlgTitle='Wavelength at which amplitude 0.5 for default detrending';
    lineNo=1;
    answer=inputdlg(prompt,dlgTitle,lineNo,def);
    wavelen = answer{1};
    if length(wavelen)>3;
        error('Wavelength more than 3 digits -- must alter program crvtbl');
    elseif length(wavelen)==1;
        wavelen=[wavelen '  '];
    elseif length(wavelen)==2;
        wavelen=[wavelen ' '];
    end;
    for j=1:7;
        str200=C(j,:);
        i = findstr(str200,'200');
        if ~isempty(i);
            strnew = strrep(str200,'200',wavelen);
            C(j,:)=strnew;
        end;
    end;
        
    
elseif ktext==3;
    C = textread('c:\mlb\stz\cfitv2.txt','%s','delimiter','\n','whitespace','');
    C =char(C);
    C=C(1:53,:);
end;
  
%--- Prompt for input file
[file2,path2]=uigetfile('*.mat','Input .mat file with stored indices, ringwidths, etc');
pf2=[path2 file2];
eval(['load ' pf2 ' S P nms yrs cmask Fwhen A ;']);



fstem = strtok(file2,'.');

kmode=3; % specifies that want both the sparse and detailed ascii summary files

% Excludes from the table (1) any cores not yet fit, as indicated by col 1 entry
% of zero in S, and (2) any masked core, as indicated by 0 in cmask

% Logical pointer to cores that have been fit
L1=S(:,1)~=0;

% Logical pointer to cores not masked
if ~islogical(cmask); error('cmask not logical'); end
L2=cmask;

% Logical pointer to cores that have been fit and are not masked
L3=L1 & L2;


% Cull nms, yrs, S for cores used in chronology
A=A(L3,:);
S=S(L3,:);
yrs=yrs(L3,:);
nms=nms(L3,:);

ncores = size(nms,1); % number of selected cores
PNUM=repmat(pnum,ncores,1);
p=repmat(NaN,ncores,1);
cshift=p;
acoef=p;
bcoef=p;
sfact=p;
arord = P(L3,4); % AR order used to create residual core index
varrat = P(L3,5); % ratio of variance of ar residuals to original index
if ktext==3; % diff index
    for jj=1:ncores;
        p(jj) = A{jj,2}; % power of transformation
        cshift(jj)=A{jj,1};
        acoef(jj)=A{jj,3};
        bcoef(jj)=A{jj,4};
        
    end;
    sfact=S(:,7); % scale factor
else; % ratio index
    p=ones(ncores,1); %  no power transformationp;
    cshift=[];
    acoef=[];
    bcoef=[];
    sfact=[];
end;

% Get string versions of power transformation parameters
Spower = repmat(blanks(4),ncores,1); % power
pmtx = repmat(p,1,12);
Lmatch=pmtx==PNUM;
sum1= sum(Lmatch');
if ~all(sum1==1);
    error('Not exactly one match between pst and tranformation');
end;
for j3 = 1:ncores;
    jthis = Lmatch(j3,:);
    ithis=find(jthis);
    Spower(j3,:)=pstr{ithis};
end


%--- SPARSE VERSION OF TABLE -- FOR ITRDB

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
s3=[s3 sflag ];


%--- Power of transform
s3a=[Spower  repmat(blanks(2),ncores,1)];


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
      if jj(j)==59;
         disp('this');
      end;
      
      % Compute period of 50% freq response
      pspline = S(jj(j),3); % spline parameter
      if S(jj(j),5)>=0.4999 & S(jj(j),5)<=0.50001; % Appropriate wavelength stored in S
         p50=S(jj(j),4);
      else;
         p50 = spline50(pspline);
         error('Spline fit, but column 5 of S is not 0.50, and I do not trust spline50.m yet');
         
      end;
      if abs(round(p50))>=1000;
         ss = sprintf('(%4.0f yr)',abs(p50));
      elseif abs(round(p50))>100;
         ss = sprintf('(%3.0f yr)',abs(p50));
      elseif abs(round(p50)) >10;
         ss = sprintf('(%2.0f yr)',abs(p50));
      elseif abs(round(p50)) >1;
         ss = sprintf('(%1.0f yr)',abs(p50));
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
s5 = [ stemp  repmat(blanks(2),ncores,1)];

%-- AUTOREGRESSIVE MODEL ORDER AND FRACTIONAL PCT VARIANCE DUE TO PERSISTENCE
s6 = num2str(arord,'   AR(%1d)  ');
s7 = num2str(1-varrat,'  %3.2f');


shd=C;

if ktext==3;
    stit=' N Core    Fit Period  Power  Curve fit    Persistence';
    sall=[s1 s2 s3 s3a  s4 s5 s6 s7];
else;
    stit=' N Core    Fit Period   Curve fit  Persistence';
    sall=[s1 s2 s3 s4 s5 s6 s7];
end;


sallchar=size(sall,2);
lentit=length(stit);
if lentit<sallchar;
   stit=[stit blanks(sallchar-lentit)];
elseif lentit>sallchar;
   sall=[sall repmat(blanks(lentit-sallchar),size(sall,1),1)];
end;
sall=[stit; sall];

sall=char(sall,blanks(3));
sall=char(sall,'Fit Period = period of ring width used in curve fitting');
sall=char(sall,'  * = subset of full period of ring width in .rwl file');
if ktext==3;
    sall=char(sall,'p = power of transformation');
    sall=char(sall,'  p>0 ring width is raised to this power');
    sall=char(sall,'  p<0 ring width is raised to this power, then sign is changed');
    sall=char(sall,'  p=0 log0 transformation');
    sall=char(sall,'  p=1 equivalent to no transformation');
else;
end;
sall=char(sall,'Curve Fit = type of model used to detrend ring width');
sall=char(sall,'  NE=modified negative exponential');
sall=char(sall,'  SL=least-square-fit straight line');
sall=char(sall,'  HL=horizontal line at mean of ring width for fit period');
sall=char(sall,'  CS=cubic smoothing spline');
sall=char(sall,'     (wavelength of 50% frequency reponse in parentheses)');
sall=char(sall,'Persistence = autoregressive prewhitening information');
sall=char(sall,'  AR(q) = order q of autoregressive model fit to core index');
sall=char(sall,'  nn  = decimal proportion of core-index variance due to ');
sall=char(sall,'     persistence; computed as 1-R, where R is the ratio of ');
sall=char(sall,'     variance of residual index to variance of standard index');
sall=char(sall,'References: ');    
sall=char(sall,'Cook, E.R., Briffa, K., Shiyatov, S., Mazepa, V. 1990.  Tree-');
sall=char(sall,'  ring standardization and growth-trend estimation. In: Methods ');
sall=char(sall,'  of dendrochronology, applications in the environmental  ');
sall=char(sall,'  sciences.  Kluwer Academic Press, Boston, 104-123.     ');
sall=char(sall,'Cook, E.R., Peters, K., 1981.  The smoothing spline: a new ');
sall=char(sall,'  approach to standardizing forest interior tree-ring width ');
sall=char(sall,'  series for dendroclimatic studies.  Tree Ring Bulletin, 41, ');
sall=char(sall,'  p. 45-53  ');
if ktext==3;
    
    sall=char(sall,'Hoaglin, D.C., Mosteller, F., Tukey, J.W.,1983. ');
    sall=char(sall,'  Understanding robust and exploratory analysis. John Wiley & ');
    sall=char(sall,'  Sons, Inc., New York. ');
else;
end;

shd=char(shd,blanks(3));
sall=char(shd,sall);

[ms,ns]=size(sall);
%-- Write itrdb info text file

if kmode==1 | kmode==3;
    
   deffn5=[fstem '.txt'];
   
   [file5,path5]=uiputfile(deffn5,'Output .txt file for export to ITRDB/NGDC ');
   pf5=[path5 file5];
   fid5=fopen(pf5,'w');
   for n = 1:ms;
      s = sall(n,:);
      fprintf(fid5,'%s\r\n',s');
   end;
   fclose(fid5);
end;


% MORE DETAILED TABLE FOR INSERTING INTO MS WORD


% Arrays of string blanks
for j2= 1:5;
    eval(['sbk' int2str(j2) ' =  repmat(blanks(' int2str(j2) '),' int2str(ncores) ',1);']);
end;

% Core id -- 
coreid = nms; % 8 char width

% Sequence no -- 
seqno = num2str(S(:,1),'%3.0f');
seqno = padher(seqno,3); % right justify, padding on left with blanks if needed


% Fit period; -- 
% s3 already holds this, including flag * if not entire period


%** TRANSFORMATION PARAMETERS

if ktext==2;
    
    % shift of ringwidt - 3
    rwshift = num2str(cshift,3);
    rwshift = padher(rwshift,3);
    
    % transformation power -- Spower, as created before
    
    % match coef for shifting -- acoef
    acoef=strjust(num2str(acoef,'%7G'));
    acoef=padher(acoef,8);
    
    % match coef for scaling -- bcoef
    bcoef=strjust(num2str(bcoef,'%7G'));
    bcoef=padher(bcoef,8);
else;
end;


%*** DETRENDING PARAMETERS

% First age-curve parameter (e.g., spline parameter  -- 13
p1 = num2str(S(:,3),'%8.4G'); % first detrending param (e.g., spline parameter if a spline fit)
p1 = strjust(padher(p1,10));

% Second age-curve parameter; 
p2 = S(:,4); 
% If spline, this is sline "length"; round to nearest whole number of years
L9=S(:,2)==9; % these were fit with spline
if any(L9);
    ss2=p2(L9);
    n9 = length(ss2);
    for j6=1:n9;
        ss2(j6)=round(ss2(j6));
    end;
    p2(L9)=ss2;
    clear  ss2 n9 j6;
end;
clear L9;
   

p2 = num2str(p2,'%8.4G'); % first detrending param (e.g., spline parameter if a spline fit)
p2 = strjust(padher(p2,10));

% Third age-curve parameter
p3 = num2str(S(:,5),'%8.4G'); % first detrending param (e.g., spline parameter if a spline fit)
p3 = padher(p3,10);

% Fourth age-curve parameter -- Scale factor for departures; one if ratio-index (no scaling).  Departures from 
% age curve are multiplied by this factor then added to 1.0 to ensure that difference index has mean 1.0 and 
% standard deviation equal to the ratio index derived from the same ring width (transformed) and age curve
p4=num2str(S(:,7),'%7.3g');
p4=padher(p4,9);


if ktext==3;
    H1=[coreid sbk1 seqno sbk2  s3 sbk1 rwshift sbk2 Spower sbk2 acoef sbk2 bcoef sbk3  s4];
    H2=[sbk1 p1 sbk1 p2 sbk1 p3  sbk1 p4 sbk2 s6 sbk1 s7];
    H=([H1 H2]);
else;
    H=[coreid sbk1 seqno sbk2  s3 sbk1 sbk3  s4 sbk1 p1 sbk1 p2 sbk1 p3  sbk1 p4 sbk2 s6 sbk1 s7];
end;
    
    
    
% Ascii file
deffn1=[fstem 'w.txt'];
[file1,path1]=uiputfile(deffn1,'MS-Word Insertable Ascii Info');
pf1=[path1 file1];
fid1=fopen(pf1,'w');

[mH,nH]=size(H);
for m = 1:mH;
    h=H(m,:);
    fprintf(fid1,'%s\n',h);
end;

    
% Header
% Note: set microsoft word font to new courier 8 to show the table lined up right
%
wordinf1 = 'To format file in ms word: (1) open the .txt file as a word doc';
wordinf2 = '(2) change font size to 8; (3) select  the last';
wordinf3 = 'header line and table body and click the table grid';
wordinf4 = '(4)use table formatting to get to look right';
wordinf5= '(5) make sure to make the superscripts super';
wordinf6 = '  ';



fclose(fid1);

% Update history
ctime=clock;
ctime=num2str(ctime(4:5));
dtime=date;
Fwhen{8,1}='crvtbl';
Fwhen{8,3}=file2;
Fwhen{8,4}=file5;
Fwhen{8,2}=[dtime ', ' ctime];
eval(['save ' pf2 ' Fwhen -append;']);


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
   
function dout = padher(din,nchar);
[m,n]=size(din);
d=max([(nchar-n) 0]);
if d==0;
    dout=din;
elseif d>0;
    paddy = repmat(blanks(d),m,1);
    dout = [paddy din];
else;
    error('specified too small a field');
end;

    


   