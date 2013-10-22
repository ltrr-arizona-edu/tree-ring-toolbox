function [Ireg,S]=eisch1(latrng,lonrng);
% eisch1: get sequential station index and string info for Eischeid daily P stations within a lat-lon box
% [istn,S]=eisch1(latrng,lonrng);
% Last revised 2-13-01
%
% Need a station pointer file for later use in pulling a time series matrix of daily stn precip.  
% Later will use day2fle1.m to make the time series matrix.  Day2fle1.m requires a list of sequential
% station numbers.  That list is built and stored by eisch1.m
%
%*** INPUT
%
% latrng (1 x 2)r  lower and upper latitude boundaries for box. E.g., [32.1  35.5]
% lonrng (1 x 2)r lower and upper longitude box .  E.g., [-120.1  -112.5]
% 
%
%*** OUTPUT
%
% Ireg (? x 1)i sequential station index into site information database for Eischeid daily P
% S (? x ?)i   station information for selected stations
%
% Prompted for name of ouput file "regstn?.mat" to store list of station indices as Ireg in
%
%*** REFERENCES -- NONE
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES
%
% Can use xyedge.m beforehand to find the latrng and lonrng suitable for a cluster of tree ring sites


%-- Get the Eischeid site info file

[file1,path1]=uigetfile('c:\data\eischeid\dlypcp\psurvey','Infile with Eischeid station list info');
pf1=[path1 file1];
eval(['load ' pf1 ' C;']);

lat1=str2num(C(:,35:40));
lon1=str2num(C(:,41:48));

L1 = lat1>=latrng(1) & lat1 <=latrng(2) & lon1>=lonrng(1) & lon1<=lonrng(2);
if any(L1);
    nsum=sum(L1);
    Ireg=find(L1);
    S=[num2str((1:nsum)')   repmat(blanks(2),nsum,1)   C(L1,:)];
   
else;
    error('No stns in range');
end;





