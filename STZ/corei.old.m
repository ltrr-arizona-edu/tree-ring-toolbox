function corei

% DESCRIPTION : corei
% Apply curve-fits saved from crvfit.m to get core indices; store
% core indices in strung-out vector IX
% Saves  IX along with previously generated yrs, nms, and S in a file
%
% INPUTS   :  NONE
%
% OUTPUTS  :  NONE
%___________________________________________________________

clear all

% Prompt for the .mat file name

% Get the .mat filename interactively; load file
[flmat,path1]=uigetfile('*.mat','Input .MAT filename ?');
pf1=[path1 flmat];
flold=flmat; % Will need this file name later
eval(['load ',pf1]);
clear flmat;


[ns,dum1]=size(nms);  % Store the size of nms string matrix

if ~exist('S'),
  error('No curves have been test-fitted yet');
end
if S(:,1)==zeros(ns,1),
  error('S is zero matrix');
  
end

disp('Please wait ...');

IX=ones(length(X),1)*NaN; % Indices will be stored in same-size 
  % vector as ring widths

nblk=5; % Maximum possible number of blocked intervals
scid=1;
while scid<=ns,
  if S(scid,1)~=0,
  disp(['Fitting Ring-width Series ',nms(S(scid,1),:)]);

  % Get the Full-length ring-width series and its years vector 
  yrv=(yrs(scid,1):yrs(scid,2))';
  xv=X(yrs(scid,3):yrs(scid,3)+yrs(scid,2)-yrs(scid,1));

  %Plot the rw data with scid core ID
  %hf0=figure('Color','w');
  %plot(yrv,xv,'b');

  % Get row index within yrv (and xv) of beginning and ending years
  % of fit-segment as marked in cfit1.m and stored in S
  xind1=find(yrv==S(scid,10));
  xind2=find(yrv==S(scid,11));

  % Get the coresponding years and ring-width data 
  yrvn=S(scid,10):S(scid,11);
  xvn=xv(xind1:xind2);

% VECTORIZE?
  xindb=zeros(2*nblk,1); 
  for i=1:2*nblk,
   if S(scid,12+i)~=0,
    if S(scid,12+i)<yrvn(1),
      S(scid,12+i)=yrvn(1);
    elseif S(scid,12+i)>yrvn(length(yrvn)),
      S(scid,12+i)=yrvn(length(yrvn));
    end
    xindb(i)=find(yrvn==S(scid,12+i));
   end
  end
  xindb(find(xindb==0))=[]; % If no blocking, makes xindb empty mtx


  if ~isempty(xindb);  % skip this section if no blocking
  for i=1:length(xindb)/2,
    lbdm=length(xindb(2*i-1):xindb(2*i));
    yrvn(xindb(2*i-1):xindb(2*i))=ones(lbdm,1)*NaN;
    xvn(xindb(2*i-1):xindb(2*i))=ones(lbdm,1)*NaN;
  end
  end
  
  % COMMENT OUT
  % Plot the ring-width series with any blocked segments blank
  %hf1=figure('Color','w');
  %plot(yrv(xind1:xind2),xvn,'b');


  
%******  APPLY CURVE FITS TO GENERATE SMOOTHED 'AGE' TREND *************

% Done by retrieving the curve-fit parameters from S(?,3:5) and
% applying appropriate curve type
 
   nfit1=S(scid,2); % curve-type for detrending
   if nfit1==9; % Spline -- in this case, wil re-fit 
     cvx1 = (cfspl(S(scid,3),yrv(xind1:xind2),yrvn,xvn))';
	elseif nfit1==1; % negative exponential -- generate using stored parameters 
	  kne=S(scid,3);
     ane=S(scid,4);
     bne=S(scid,5);
	  t = yrvn-yrvn(1) + 1;
     cvx1 = kne + ane * exp(-bne*t);
	  cvx1=cvx1'; % needs to be col vector
	  % Fits data  to the equation
	  % g(t) = k + a * exp(-b*t),
	  %   where t is the shifted time variable t = yrvn-yrvn(1)+1
	  %   In other words, t is same length as yrvn after 
	  %   dropping NaNs

   elseif nfit1==4; % Horizontal line through mean
      cvx1 = repmat(S(scid,3),length(yrvn),1); %cfmean(yrv(xind1:xind2),yrvn,xvn);
   elseif nfit1==2; % straight line, any slope
      t=yrvn-yrvn(1)+1;
      Xtemp = [ones(length(t),1) t'];
      aconst = S(scid,3);
      aregr = S(scid,4);
      cvx1  = Xtemp * [aconst aregr]';
   else;
      error('Unnacceptable value for nfit1');
      
   end
   
   % COMMENT OUT
   %figure(hf1);hold on;
   %plot(yrv(xind1:xind2),cvx1,'r');
   %pause
   %************

% Compute  index -- by ratio
rwin1=xvn./cvx1;

% Adjust ring-width index by shifting to mean 1.0
mnx = nanmean(rwin1);
rwin2 = rwin1 -mnx + 1.0;

% Store index
IX(S(scid,12):S(scid,12)+S(scid,11)-S(scid,10))=rwin2;

end
 scid=scid+1;
end		% End of scid while loop


%******************** AR MODELING AND GENERATION OF RESIDUAL CORE INDICES

nnmx=input('Maximum AR order to consider in prewhitening model: ');
[EX,P]=armod(IX,S,yrs,nms,nnmx);

%******************* STATISTICS OF CORE INDICES *********************

[CE1,CE2,CE3,CS1,CS2,CS3]=stats1(IX,EX,yrs,nms);


%*********************************


% Save the vectors in a .mat file
setnew=' IX EX P CE1 CE2 CE3 CS1 CS2 CS3 '; % New variables to be saved
nsv=menu('Save variables?','Add New Variables to original file?',...
   'Store new variables only in a new file?','QUIT');
if nsv==1,
  eval(['save ' pf1 setnew ' -append']);
elseif nsv==2,
   [ofmat,path2]=uiputfile('*.mat','New .MAT filename to store new stuff in: ');
   pf2=[path2 ofmat];
  eval(['save ',pf2, setnew]);
end


% End of file
