function cfit2

% DESCRIPTION : cfit2
% Apply curve-fits saved from cfit1.m to get core indices; store
% core indices in strung-out vector IX
% Saves  IX along with previously generated yrs, nms, and S in a file
%
% INPUTS   :  NONE
%
% OUTPUTS  :  NONE
%___________________________________________________________



% Prompt for the .mat file name

% Get the .mat filename interactively
flmat=uigetfile('*.mat','.MAT filename ?');

% Check if the mat file exists in the workspace
if ~exist('X'),
  eval(['load ',flmat]);
end

[ns,dum1]=size(nms);  % Store the size of nms string matrix

if ~exist('S'),
  udisp('No curves have been test-fitted yet');
  return;
end

if S(:,1)==zeros(ns,1),
  udisp('S is zero matrix');
  return;
end

disp('Please wait ...');

IX=ones(length(X),1)*NaN; % Indices will be stored in same-size 
  % vector as ring widths

nblk=5; % Maximum possible number of blocked intervals
scid=1;
while scid<=ns,
 if S(scid,1)~=0,
  disp(['Fitting RingWidth Series ',nms(S(scid,1),:)]);

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

  % Plot the ring-width series with any blocked segments blank
  %hf1=figure('Color','w');
  %plot(yrv(xind1:xind2),xvn,'b');


  
%******  FIRST CURVE FITTING *************

   nfit1=S(scid,2); % curve-type for first detrending
   nfit2=S(scid,6); % curve-type for second detrending
   if nfit1==9; % Spline
     cvx1 = cfspl(S(scid,3),yrv(xind1:xind2),yrvn,xvn);
   else; % Horizontal line through mean
     cvx1 = cfmean(yrv(xind1:xind2),yrvn,xvn);
   end

%  figure(hf1);hold on;
%  plot(yrv(xind1:xind2),cvx1,'r');

% Compute first index -- by ratio
   rwin1=xvn./cvx1;

%  hf2=figure('Color','w');
%  plot(yrv(xind1:xind2),rwin1,'r');
%  ht1=text('Position',[0.5,0.9],'Color','k','String',...
%  [nms(scid,:),'- Index after 1st fit'],'Units','Normalized');



%*************  SECOND DETRENDING  ******************************

   if nfit2==9; % Spline
      cvx2 = cfspl(S(scid,7),yrv(xind1:xind2),yrvn,xvn);
      rwin2 = rwin1 ./ cvx2;
   else; % Other than spline as second curve-fit
      rwin2=rwin1;
   end

%  figure(hf2);hold on;
%  plot(yrv(xind1:xind2),cvx2,'r');

% Accept the current 2nd detrending ?
%  hf3=figure('Color','w');

% Adjust second index to mean 1.0
  rwin2 = rwin2 + 1.0 - mean(rwin2(~isnan(rwin2)));


%  plot(yrv(xind1:xind2),rwin2,'r');
%  text('Position',[0.5,0.9],'Color','k','String',...
%  [nms(scid,:),'- Index after 2nd fit'],'Units','Normalized');



% Store the second index
  IX(S(scid,12):S(scid,12)+S(scid,11)-S(scid,10))=rwin2;

%  pause;
%  close all;
end
 scid=scid+1;
end		% End of scid while loop

nsvf=usinp('Save index data in file ?');
if nsvf,
  flsn=uiputfile('*.mat','Save as ?');
  eval(['save ' flsn ' X ' ' IX ' ' yrs ' ' nms ' ' S ']);
end


% End of file
