function matdat01(nfiles)
% matdat01: ascii .dat files of PDSI for steve leavitt from watbalm1.m output
% CALL: matdat01(nfiles);
%
% Meko 10-3-98
%
%************* IN 
%
% nfiles (1 x 1)i number of watbalm1.m output files to get data from. These
%    assumed to be sequentially number wbout1.mat, wbout2.mat, ... and to be int
%    the current directory
%*********** OUT
%
% Ascii files named pdsi1.dat, pdsi2.dat, etc
%
%********** NOTES
%
% Input files assumed to be generated from watbalm1.m.  The PDSI is datout{2} in 
% the output storage 

fmt1='%4.0f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n';

for n = 1:nfiles;
   flnm1 = ['wbout' int2str(n)];
   eval([' load ' flnm1 ';']);
   Z = datout{2};
   
   [m1,n1]=size(Z);
   
   flnm2 = ['pdsi' int2str(n) '.dat'];
   fid1=fopen(flnm2,'w');
   fprintf(fid1,fmt1,Z');
   fclose(fid1);
end


