function  p = respoint(A,w,Ac)
% respoint:  period at which a filter has desired amplitude of frequency response
% CALL: p = respoint(A,w,Ac);
%
% D Meko 11-5-92
% Revised  2-24-93.
%
%*********  INPUT ARGS
% 
% A (mA x 1)  :  Magnitude of frequency response function
% w (mA x 1)  :  Vector of frequencies corresponding to A
% Ac (1 x 1)  :  A particular magnitude of frequency response for which
%		you want the corresponding wavelength
%
%*********  OUTPUT ARG
%
% p (1 x 1) wavelength (yr) corresp to the Ac value of freq response
%
%**********  METHOD
%
% Interpolation.  Given a frequency response function.  Look down the
% vector of mangnitudes -- starting at the lowest frequency -- until
% the first magnitude less than or equal to Ac.  Interpolate between this
% frequency and the next frequency in the vector.  Finds the first 
% satisfactory point from the left (low-frequency) side.
%
%*************  CAUTION
%
% Written for low-pass filters.  Assumes that first occurence of Ac will
% be on a falling limb of the frequency response funtion.

i=find(A>Ac);  % Find index of all amplitudes greater than the critical
id=diff(i);  % First difference;  Will equal 1 as long as monotonic 
L=id~=1;  %  If not monotonic decreasing, some elements of L will not  
%		equal 1

if (~any(L));   % All elements of id are 1;  all elements of L are 0.
	%  In other words, the freq response is monotonic decreasing.
	i=length(L)+1;
	j=i+1;
else  %   A has multiple peaks -- A is not monotonic decreasing. 
% Use the first leftmost occuring A interval including Ac to find
% the wavelength for Ac
	Li=find(L==1)
	i=min(i);
	j=i+1;
end

fract=(Ac-A(j))/(A(i)-A(j));
wc=w(j)  -  fract * (w(j) - w(i));
p=1.0 / wc;  % Convert frequency to wavelength
