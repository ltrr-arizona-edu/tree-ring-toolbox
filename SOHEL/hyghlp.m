function [m,N,ps,n95,n99] = hyghlp(my,Lx,Ly,klq)
%
% USAGE : [m,N,ps,n95,n99] = hyghlp(my,Lx,Ly,klq)
%   Computes the hypergeometric statistics as described in HYGEOM1.M
%
%
% INPUTS
%-------
% my (1 x 1)	Segment length for x time series
% Lx (1 x ny)	Logical Pointer for x matrix
% Ly (1 x ny)	Logical Pointer for y matrix
% klq (1 x 1)	If klq=1, the critical values n95 and n99 are computed
%		If klq~=1, zeros are returned for n95 and n99.
%
%
% OUTPUTS
%--------
% N (nx x 1) 	The number of "samples", defined as above as the 
%		number of x values below the single-sided threshold 
%		or outside the double-sided threshold
% m (ny x 1) 	The number of successes.  A success is a sample value
%  		for which y also is below its single-sided threshold
%		or outside the double-sided threshold
% ps (1 x ny) 	The number of y's outside qu and/or ql percentile
% n95 (1 x ny)	95 % significance level
% n99 (1 x ny) 	99 % significance level
%
%
% NO USER WRITTEN FUNCTIONS NEEDED
%__________________________________________________________________

 m = sum(Lx & Ly);
 ps = sum(Ly)/my;
 N = sum(Lx);
 if klq,
   n95=hygeinv(0.95,my,ps*my,N(1));
   n99=hygeinv(0.99,my,ps*my,N(1));
 else
   n95=zeros(1,length(ps));
   n99=zeros(1,length(ps));
 end

% End of file
