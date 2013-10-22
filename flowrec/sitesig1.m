function sitesig1
% sitesig1:  ascii summary table of sites used in model
% CALL: sitesig1;
%
% Meko 8-3-98
%
%********** IN
%
% No input args.
% User prompted to select several files:
%   -sitenms.mat: string matrix, 50 chars wide, Snames holds info on:
%       seq number of site, crn flname, sitename, state, species 
%   -sitelst?.txt ascii of site information.  From this will lat,long,el
%   -gosp?.dat  start and end year of standard chron
%   -pty?.mat  ptyr holds "ass-year", adjusted further of any loss due to AR startup
%   -climwhch.mat
%       Lomit (? x 1)L  1 if chron to be omitted from recon modeling, 0 otherwise
%       iclim(? x 1)i  which water-balance variable represented by single-site filter
%              1=pcp, 4=Z index
%   -resar?.mat: output stats from single-site filtering
%      C1 (? x 16)r R-sqd
%      C3 (? x 16)r  p-value for overall F
%         Note.  the 16 cols refer to different water-balance variables. For example,
%         col 1 is pcp, col 4 is Z index
%
%
%******** OUTPUT
%
% Ascii table with the following:
%
% cols 1-50 sitenms info (see above)
% col  51 blank
% cols 52-56 lat, deg and min
% cols 59-65 long
% col  67 -- Lomit
% col  69 -- iclim
% col  72-81    R-sq and p value of overall F:    0.65(1.3E-5)


%******************   Get sequence number and initial identifying info
[file1,path1]=uigetfile('sitenms.mat','.mat file with Snames info');
pf1 = [path1 file1];
eval(['load ' pf1 ';']);
%   now have Snames
str1 = Snames;


nsites = size(Snames,1);

parleft=repmat('(',nsites,1);
parright = repmat(')',nsites,1);
blnkcol = repmat(' ', nsites,1);


%******************   Get long, lat, el
[file2,path2]=uigetfile('sitelst?.txt','.txt file with site info');
pf2 = [path2 file2];
A = caseread(pf2);

latdeg = A(:,65:66);
latmin = A(:,72:73);
londeg = A(:,78:81);
lonmin = A(:,88:89);
elm = A(:,91:94);

str2 = [latdeg  blnkcol latmin  blnkcol londeg lonmin blnkcol elm blnkcol];

clear A;


%******************   Get start and end year of standard chron
[file3,path3]=uigetfile('gosp?.dat','aosp?.dat file with chron start and end yr');
pf3 = [path3 file3];
eval(['load ' pf3 ';']);

str3 = char(num2str(Gospa));

%************* Get pointer year

[file4,path4]=uigetfile('ptyr?.mat','.mat file with pointer years of sites');
pf4=[path4 file4];
eval(['load ' pf4 ';']);

str4 = char(num2str(ptyr));


%************  Get Lomit and iclim
[file5,path5]=uigetfile('climwhch.mat','.mat file with Lomit and iclim');
pf5=[path5 file5];
eval(['load ' pf5 ';']);

str5= char(num2str(Lomit));
str6 = char(num2str(iclim));

%********* Single-site filtering stats
[file1,path1]=uigetfile('resar?.mat','stats from single-site filtering');
pf1 = [path1 file1];
eval(['load ' pf1 ';']);

siggy = repmat(blanks(3),nsites,1);
c4 = repmat(NaN,nsites,1);

CC1 = repmat(NaN,nsites,1);
for n = 1:nsites;
   c1 = C1(n,iclim(n));
   c4(n)=c1;
   c3 = C3(n,iclim(n));
   if c3<0.001;
      siggy(n,:)='***';
   elseif c3< 0.01;
      siggy(n,:) = '** ';
   elseif c3 < 0.05;
      siggy(n,:) = '*  ';
   else
   end
end

  
str7 = [char(num2str(c4))  parleft siggy parright];
  

S = [str1 str2 str3 blnkcol str4 blnkcol str5 blnkcol str6 blnkcol str7];


%******** OUTPUT

[file2,path2]=uiputfile('sslist?.txt','.txt file to hold output');
pf2 = [path2 file2];
casewrite(S,pf2);