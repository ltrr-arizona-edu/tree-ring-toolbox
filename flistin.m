function [S,len]=flistin()

% Read in a list of file names, for use by rwlook
%
% D. Meko 12-22-93
%
%
%************ INPUT ARGS
%
%***********  OUTPUT ARGS
%
% S (? x 12) list of file names
% len (? x 1) number of valid chars in each row of S
%
%
%***************  NOTES ***********************************
%
% Used with rwlook to avoid having to type in the entire rw file name
% each time you want to read it in.
%
% File of file names must have 'eof' as last line

ff=input('Name of file holding names of sub-files: ','s');
fid=fopen(ff,'r');

k=1;
while k<100; % default in case forget eof line
	eval(['s=fgetl(',int2str(fid),');']);
	if s(1:3)=='eof';, break, end;
	lens=length(s);
	len(k)=lens;
	need=12-lens;
	if need>0
		b=' ';
		b=b(:,ones(need,1));
		S(k,:)=[s b];
	else 
		S(k,:)=s(1:12);
	end

	
	k=k+1;
end



	
