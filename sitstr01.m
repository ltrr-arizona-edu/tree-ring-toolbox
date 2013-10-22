function S=sitstr01(Lcull,infile,nlast)
% sitstr01: make string matrix from ascii matrix of site info
% CALL: sitstr01;
%
% Used in Sacramento River reconstruction to get string matrix with
% sequential number, filename, site name, state, and species
%
%************  IN ************
% 
% Lcull -- logical vector to desired rows in text file infile
% infile -- asci text file with site information
% nlast -- greatest number of chars needed from any row of infile
%
%********** OUT *********
%
% S (? x ?)s string matrix with form:
%  1-3 sequential site number
%  5-9 original sequence number (in parens)
%  11-22 file name
%  24-37 first 14 chars of site name
%  39-42 species code

clc

% Get the input file
fid1 = fopen(infile,'r');


% ****************** Count rows of input file
k1=1;
nsites=0;
while k1
   c=fgetl(fid1);
   if ~(feof(fid1) & length(c)<5);
      nsites=nsites+1;
   else
      k1=0;
   end
end

if length(Lcull)~=nsites;
   fclose all
   error('Lcull should have length equal to nsites');
end


disp(['Total of ' int2str(nsites) ' counted in ' infile]);

%****************** Size storage matrix
blnks = blanks(nlast);
Sin=repmat(blnks,nsites,1);


%*************** Fill storage matrix

frewind(fid1);
for n = 1:nsites
   c=fgetl(fid1);
   Sin(n,:)=c(1:nlast);
end

% Cull lines of storage matrix
Sin=Sin(logical(Lcull),:);

ngood =sum(Lcull);  % number of culled sites
S = repmat(blanks(50),ngood,1);

for n = 1:ngood;
   sthis = Sin(n,:);
   str1 = sprintf('%2.0f ',n);
   ctemp = sthis(3:31);
   str2 = sprintf('%s ',ctemp); % file name and 14 chars of site name
   ctemp = sthis(41:48);
   str3 = sprintf('%2s-',ctemp); % state
   ctemp =sthis(98:101);
   str4 = sprintf('%4s ',ctemp);  % species
   strall = [str1 str2 str3 str4];
   S (n,1:47)= strall;
end

fclose all;