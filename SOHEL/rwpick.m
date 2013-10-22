function [x,yx,nmx]=rwpick(zname,z,zy,zind,svr)
%
% USAGE : [x,yx]=rwpick(zname,z,zy,zind,ncol)
%   Picks a data series x out of a strung out vector z and
%   also returns the corresponding year vector yx
%
%
% INPUTS
%-------
% zname 	String matrix of Core IDs.
% z (? x 1)	The data array in a strung-out column vector.
% zy (? x 3)	Year index
% zind (? x 1)	Index matrix to the strung out vector z
% svr 		String variable 'X' or 'Y'.
%     
%   
% OUTPUTS
%--------
% x (? x 1)	Time series vector x
% yx (? x 1)	Year vector corresponding to x
% nmx  		Core ID of chosen series x
%
%
% USER WRITTEN FUNCTIONS NEEDED 
%------------------------------
% SVMENU.M	A modified menu function
%__________________________________________________________

% Prompt for the data series .
mntl=['Choose series : ',svr];
n = dvmenu(mntl,zname);
   
% Display the chosen series Core ID
    disp(['You have chosen this series : ', zname(n,:)]);

    nmx=zname(n,:);
    yx=(zy(n,2):zy(n,3))';
    x=z(zind(n):zind(n+1)-1);

% End of file
