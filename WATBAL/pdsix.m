function [x,xm]=pdsix(z)
% pdsix: Compute PDSI from z index, as in Table 12, Palmer, 1965
% CALL: [x,xm]=pdsix(z);
%
%******************************  IN **********************************
%
% z (mz x 1)r   Palmer 'Z index' (see notes)
%
%******************************  OUT *********************************
%
% x PDSI, same length as z
% xm modified PDSI (see in-text comments)
%
%******************************** NOTES ***************************
%
% z is assumed to be a multiple of 12; the data is assumed to begin
% with jan of the first year and end with dec of the last year
%
% See Palmer, Table 12, p. 31 for info on variables


a=NaN;

% Check Input
[mz,nz]=size(z);
if nz~=1;
	error('z must be col vect');
end
if rem(mz,12)~=0;
	error('length of z must be multiple of 12');
end

%**************************** Size variables

Uw=repmat(a,mz,1);
Ud=repmat(a,mz,1);
V=zeros(mz,1);
Pe=zeros(mz,1);
Q=repmat(a,mz,1);
Ze=repmat(a,mz,1);
Q=repmat(a,mz,1);
x1=zeros(mz,1);
x2=zeros(mz,1);
x3=zeros(mz,1);


%************************ Easy vector computations

z3=z/3;


%************************ Allocate

LLd=repmat(logical(0),mz,1); % 1 means established drought 
LLw=repmat(logical(0),mz,1); % 1 means established wet period
LLn=repmat(logical(1),mz,1); % 1 means neither established dry or wet period occurring
pullx1=repmat(logical(0),mz,1);
pullx2=repmat(logical(0),mz,1);
pullnorm=repmat(logical(0),mz,1);


%*****************************************************************
% Initialize values for first month; assumes x1,x2, x3==0 in previous month

Uw(1)=NaN; % effective wetness; applies only with pre-existing drought
Ud(1)=NaN; % effective dryness; applies only with pre-existing wet period
V(1)=0; % numerator term for prob(end of drought or wetness)
Ze(1)=NaN; % z-value needed to end drought or wet period in one month
Q(1)=NaN; % denominator term ...
Pe(1)=0; % prob that drought or wet period has ended this month

x1(1)=max(0,z3(1));
x2(1)=min(0,z3(1));
x3(1)=z3(1);

if z3(1)<=-1.0;
	status='dry';
elseif z3(1)>=1.0;
	status='wet';
else
	status='normal';
end

switch status
case 'wet';
	LLw(1)=1; % Logical marker for a wet period month
	newwet=1; % Have started new wet period
	firstdry=0; % Have not yet encountered first effective dry month since start of wet per
	nump=0; % Number of active previous months in the probability computation
case 'dry'
	LLd(1)=1;
	newdry=1;
	firstwet=0;
	nump=0;
case 'normal';
	LLn(1)=1;
	x3(1)=0;
	newnorm=1;
	pullnorm(1)=1;
end

for i =2:mz;
	
	switch status

%***************
	case 'normal'
		x1(i)=max(0,0.897*x1(i-1)+z3(i));
		x2(i)=min(0,0.897*x2(i-1)+z3(i));
		x3(i)=0;
		if x1(i)<1.0 & x2(i)>-1.0; % no new drought or wet perod
			LLn(i)=1;
			pullnorm(i)=1;
		elseif x1(i)>=1.0; % new wet period
         status='wet'; newwet=1; firstdry=0;  nump=0;
         LLw(i)=1;
			x3(i)=x1(i);
		elseif x2(i)<=-1.0;
         status='dry'; newdry=1; firstwet=0;  nump=0;
         LLd(i)=1;
         x3(i)=x2(i);
		end
		

%*************

	case 'wet'
		x1(i)=max(0,0.897*x1(i-1)+z3(i));
		x2(i)=min(0,0.897*x2(i-1)+z3(i));
		x3(i)=0.897*x3(i-1)+z3(i);
		Ud(i)=z(i)-0.15;
		if newwet==1;
			newwet=0;
			x1(i)=0;
		end
		if firstdry==0; % have not previously hit an effective dry Ud
			if Ud(i)>=0;
				Ud(i)=NaN;
				Pe(i)=0;
			else
				firstdry=1;  % flag that now have hit first effect dry Ud
				Ze(i)=-2.691 * x3(i-1) + 1.50;  % z-value needed to end wet
      		% Compute prob(drought has ended)
				V(i)=Ud(i);
         	Q(i)=Ze(i);
				Pe(i)=V(i)*100/Q(i);
			end
		else; % firstdry==1
			Ze(i)=-2.691 * x3(i-1) + 1.50;  % z-value needed to end wet
      	% Compute prob(drought has ended)
			V(i)=V(i-1)+Ud(i);
        Q(i)=Ze(i)+V(i-1);
			Pe(i)=V(i)*100/Q(i);
		end; % if firstdry
     if Pe(i)<0;
        Pe(i)=0;
     elseif Pe(i)>=100;
        Pe(i)=100;
     end
      
		if Pe(i)==0;  % zero prob that wet period over
         LLw(i)=1; nump=0;
         x1(i)=0;
			if Pe(i-1)>0; % if fizzled out ending, zero x1 and x2
				x1(i)=0;
            x2(i)=0;
            firstdry=0; Ze(i)=NaN; V(i)=0; Q(i)=NaN; nump=0;
			end
		elseif Pe(i)==100;  % wet period ends
			if x2(i)<=-1.0;  % entering a new drought
				status='dry';  LLd(i)=1;  newdry=1;  firstwet=0;
            x3(i)=x2(i);
            pullx2(i:-1:(i-nump))=1;
         elseif x1(i)>=1.0;
				%error('Ending wet period and entering new wet period');
			else; % x1 and x2 do not indicate starting drought or wet period
				status='normal'; LLn(i)=1;  
				x3(i)=0;
				pullnorm(i:-1:(i-nump))=1;
         end
         nump=0;
		else; % Pe(i) is greater than 0 but less than 1; wet period continues
			nump=nump+1;
		end


% *************

	case 'dry'
		x1(i)=max(0,0.897*x1(i-1)+z3(i));
		x2(i)=min(0,0.897*x2(i-1)+z3(i));
		x3(i)=0.897*x3(i-1)+z3(i);
		Uw(i)=z(i)+0.15;
		if newdry==1;
			newdry=0;
			x2(i)=0;
		end
		if firstwet==0; % have not previously hit an effective wet Uw
			if Uw(i)<=0;
				Uw(i)=NaN;
				Pe(i)=0;
			else
				firstwet=1;  % flag that now have hit first effect wet Uw
				Ze(i)=-2.691 * x3(i-1) - 1.50;  % z-value needed to end dry
      		% Compute prob(dry period has ended)
				V(i)=Uw(i);
         	Q(i)=Ze(i);
				Pe(i)=V(i)*100/Q(i);
			end
		else; % firstwet==1
			Ze(i)=-2.691 * x3(i-1) - 1.50;  % z-value needed to end dry
      	% Compute prob(drought has ended)
			V(i)=V(i-1)+Uw(i);
         Q(i)=Ze(i)+V(i-1);
			Pe(i)=V(i)*100/Q(i);
		end; % if firstwet
     if Pe(i)<0;
        Pe(i)=0;
     elseif Pe(i)>=100;
        Pe(i)=100;
     end
      
		if Pe(i)==0;  % zero prob that dry period over
         LLw(i)=1; nump=0;
         x2(i)=0;
			if Pe(i-1)>0; % if fizzled out ending, zero x1 and x2
				x1(i)=0;
            x2(i)=0;
            firstwet=0; nump=0; V(i)=0; Q(i)=NaN; Ze(i)=NaN;
			end
		elseif Pe(i)==100;  % drought ends
			if x1(i)>=1.0  % entering a new wet period
				status='wet';  LLw(i)=1;  newwet=1;  firstdry=0;
				x3(i)=x1(i);
            pullx1(i:-1:(i-nump))=1;
         elseif x2(i)<=-1.0;
				%error('Ending drought and entering new drought!');
			else; % x1 and x2 do not indicate starting drought or wet period
				status='normal'; LLn(i)=1;  
				x3(i)=0;
				pullnorm(i:-1:(i-nump))=1;
         end
         nump=0;
		else; % Pe(i) is greater than 0 but less than 1; drought continues
			nump=nump+1;
		end

	end; % case 'normal', 'wet', or 'dry'

end ; %  for i = 2:mz



%******************************************************************

% column vectors of monthly data x1,x2, and x3 have now been calculated.  Now must
% substitute values from series x1 and x2 for those in x3 when appropriate

x=x3; % initialize storage vector of final pdsi
if any(pullx1);  x(pullx1)=x1(pullx1); end; % Handle background wet periods
if any(pullx2);  x(pullx2)=x2(pullx2); end;  % Handle background droughts
% Handle normal periods.  Substitute into x3 the whichever of x1 and x2 has
% the max absolute value
Lwet = pullnorm & abs(x1)>=abs(x2);
Ldry = pullnorm & abs(x2)>abs (x1);
if any(Lwet);   x(Lwet)=x1(Lwet);end;
if any(Ldry);  x(Ldry)=x2(Ldry);  end;

% Modified Palmer Index.
% The series xm is a modified index that may differ from x when the 
% probability of ending a drought or wet period is greater than zero but
% less than 100%.  For such months, the modified index is a weighted average
% of x3 and either x1 or x2, depending on whether the existing period as
% indicated by x3 is a drought (x3<0) or a wet period (x3>0).  The weights are
% the probabilities (Pe) of ending the drought or wet period. If this probability
% is  prob=100 * Pe, the weighted value for dry period is:
%
%  x = prob*x2 + (1-prob)* x3;
%
% For a wet period, it is:
%
%  x = prob*x1 + (1-prob)*x3

prob= Pe/100; % prob of ending a drought or wet period
L7 = prob>0 & prob<1;
L7d = L7 & x3<0;
L7w = L7 & x3>0;

xm=x;
if any(L7d);
   xm(L7d)= prob(L7d) .* x1(L7d) + (1-prob(L7d)) .* x3(L7d);
end

if any(L7w);
   xm(L7w)=prob(L7w) .* x2(L7w) + (1-prob(L7w)) .* x3(L7w);
   
end

