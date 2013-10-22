function flag=skipln1(fid,n)
% skipln1: skip lines in reading a file
% CALL: flag=skipln1(fid,n);
%
% Meko 5-15-97
%
%************** IN *************
%
% fid -- file identifier of file being read
% n  -- number of lines to skip
%
%
%************** OUT *************
%
% flag -- error status flag
%		==0 no problem
%		==1 eof reached in skipping


if n>0;
	for nn = 1:n;
		if feof(fid);
			flag=1;
			return
		else
			flag=0;
			c=fgetl(fid);
		end
	end
else
	return;
end




