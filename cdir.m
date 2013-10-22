function cdir
% cdir: change current directory to one of listed directories
% CALL: cdir
%
% Meko 9-16-97
%
%*********  Notes 
%
% Revise txt1,2,etc, and txtm as needed


txt1 = 'cd c:\projs\ac4\mapsmlb';
txt2 = 'cd c:\data\ghcn\tmpv2_0\dolt';
txt3 = 'cd c:\projs\ac4\tmp';
txt4 = 'cd c:\projs\ac4\pdsi';
txt5= 'cd c:\projs\af1\sacrflow';

txtm = '1 -- cd c:\projs\ac4\mapsmlb';
txtm = strvcat(txtm,'2 -- cd c:\data\ghcn\tmpv2_0\dolt');
txtm = strvcat(txtm,'3 -- cd c:\projs\ac4\tmp');
txtm = strvcat(txtm,'4 -- cd c:\projs\ac4\pdsi');
txtm = strvcat(txtm,'5 -- cd c:\projs\af1\sacrflow');
disp(blanks(5));
disp(txtm)

k = input('Change to directory of which of above: ');

switch k
case 1
   eval(txt1);
case 2
   eval(txt2);
case 3
   eval(txt3);
case 4
   eval(txt4);
case 5 
   eval(txt5);

otherwise
end



