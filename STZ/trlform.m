function trlform(X,yrs,n,id)
% USAGE : trlform(X,yrs,n,id);
%   Writes a data file in Tree Ring Lab format. The output file has 
%   the following format : 1st 6 columns of the matrix contain
%   Core ID strings. Next 4 columns contain years in decades. Column 
%   11 and up hold the data for individual years and their corresponding
%   sample size. 
%
% INPUTS  :  X (? x 1) - Strung out vector containing tree ring info
%	     yrs (? x 1) - Year vector corresponding to X
%	     n (? x 1) - Sample size vector corresponding to X
%	     id (1 x 6) - Core ID for the standard tree
%
% NO OUTPUTS
%
% USER WRITTEN FUNCTIONS NEEDED
%------------------------------
% UDISP.M : A user interactive graphic display function
% JDISP.M : A non-interactive graphic display function
%
% FILES PRODUCED
% --------------
% ?.MAT		A .mat file with user given name. Contains 
%		a data matrix, a data ID matrix, and a text
%		matrix explaining different columns of the 
%		data ID matrix.
%________________________________________________________
    

% Prompt for data file name
fln=uiputfile('*.crn','OUTPUT filename');

% If Cancel button is pressed in the input window, flag is 
% set to -1 and control is returned to calling function
if fln==0,
  fl=-1;
  return;
end 

% Display text window
jh=jdisp(['Writing into the file: ',fln,'. Please wait ... ']); 
pause(1);

% Open the data file for writing
fid=fopen(fln,'wt');

% Rewind the file to the top
frewind(fid);

% Write the string data at the top of the file
smtxt=[id,' : Run ITRDB program HEADER to add header lines'];
fprintf(fid,'%s   \n',smtxt);

% Convert all the NaN's in X into 9990's
lmx=isnan(X);
X(lmx)=ones(sum(lmx),1)*9.990;
n(lmx)=zeros(sum(lmx),1);

bcount=1;
while lmx(bcount)==1,
  bcount=bcount+1;
end
bcount=bcount-1;
rcount=1;
k=length(lmx);
while lmx(k)==1,
  k=k-1;
  rcount=rcount+1;
end
rcount=rcount-1;

% Append 9990's in the beginning of the data if necessary
if rem(yrs(1),10)~=0,
  fapnd=ones(rem(yrs(1),10),1)*9990;
end
bfcount=bcount+length(fapnd);

% Append 9990's at the end of the data if necessary
if rem(yrs(length(yrs)),10)~=0,
  rapnd=ones(9-rem(yrs(length(yrs)),10),1)*9990;
end
rfcount=rcount+length(rapnd);

% Append 9990's in the beginning of the data if necessary
brem=rem(bfcount,10);
if brem~=0,
  bwhln=fix(bfcount/10);
  if bwhln==0,
    fprintf(fid,'%s',id);
    fprintf(fid,'%4d',yrs(1));
    for i=1:length(fapnd),
      fprintf(fid,'%4d%3d',fapnd(i),0);
    end
  else
    fprintf(fid,'%s',id);
    fprintf(fid,'%4d',yrs(bwhln*10+brem-length(fapnd)+1));
    for i=1:brem,
      fprintf(fid,'%4d%3d',9990,0);
    end
  end
end

disp('Writing the data. Please wait...');

% Main loop : write the actual data X
for i=bcount+1:length(X)-rcount,
  if rem(yrs(i),10)==0,
    fprintf(fid,'\n%s%4d',id,yrs(i));
  end
  fprintf(fid,'%4d%3d',round(X(i)*1000),n(i));
end

% Append 9990's at the end of the data if necessary
rrem=rem(rfcount,10);
if rrem~=0,
  for i=1:rrem,
    fprintf(fid,'%4d%3d',9990,0);
  end
end

% Close the display window
close(jh);

fclose(fid);  	% Close the input file


% End of file
