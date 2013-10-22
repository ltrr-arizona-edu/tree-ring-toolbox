function rvout=dlgi001(tit,prompt,lineno,rvin)
% dlgi001: input dialog to edit elements of row vector
% CALL: rvout=dlgi001(tit,prompt,lineno,rvin);
%
% Meko 5-3-97
%
%************************ IN **************
%
% tit (1 x ?)s  title for dialog box
% prompt{}   cell of character prompt
% lineno (1 x 1)i  number of lines in prompt, and in response
%   *** function not generalized to any except 1 line ****
% rvin (1 x ?)r    row vector to be changed
%
%
%****************************** OUT 
%
% rvout (1 x ?)r  altered row vectro
%
%**************************************



if ~ischar(tit)
	error('Dialog title must be character');
end
if ~iscell(prompt);
	error('Dialog prompt must be cell');
end
if ~isnumeric(lineno);
		error('Dialong line numbers variable must be numeric')
end


if ~isnumeric(rvin);
	error('rvin must be numeric');
end

[m1,n1]=size(rvin);
if m1~=1;
	error('rvin must be rv');
end

% change default answer from numeric rv to cell of characters
deflt = {num2str(rvin)};


ans1=inputdlg(prompt,tit,lineno,deflt); % class is cell
ans1=(cat(1,ans1{:})); % class is char
rvout=str2num(char(ans1)); % class is numeric

