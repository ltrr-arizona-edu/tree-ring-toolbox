function sitlst01(Lcull,infile,nlast)
%
% Lcull -- logical vector to desired rows in text file infile
% infile -- asci text file with site information
% nlast -- greatest number of chars needed from any row of infile
clc

% Get the input file
fid1 = fopen(infile,'r');


% ******************8Count rows of input file
k1=1;
nsites=0;
while k1
   c=fgetl(fid1);
   clen = length(c);
   if ~feof(fid1) | clen>10;
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
S=repmat(blnks,nsites,1);



%*************** Fill storage matrix

frewind(fid1);
for n = 1:nsites
   c=fgetl(fid1);
   S(n,:)=c(1:(nlast));
end



% Cull lines of storage matrix
S=S(Lcull,:);


ngood =sum(Lcull);  % number of culled sites


%*******************  ASCII MATRIX OF SITE INFO, CULLED SITES

[file2,path2]=uiputfile('*.txt','Output site list');
pf2=[path2 file2];

for n = 1:ngood
   str1 = sprintf('%3.0f ',n);
   str2 = S(n,:);
   str3=[str1 str2];
   fprintf(pf2,'%s\n',str3);
end
fclose all
   

disp(' ')
   