function [R,porder,torder,pratio,tratio] =...
corr3(C1,D,E,L1,endmo,arord,k,kw,nhi)

% Returns matrix of corr coefs in R.  Elements are product-moment 
% correls between different types of series, depending on k.
%
% k=5  ppt vs temp
% k=6  tree.ppt  vs temp.ppt,  where y.x means residuals of linear
%		regression of y against x.  
% k=7  tree vs whitened ppt
% k=8  tree vs whitened temp
% k=9  whitened ppt vs whitened temp  (whitened version of k=5)
% k=10  tree.(whitened ppt)  vs   (whitened temp).(whitened ppt)

% C1=tree has been previously prewhitened--the residuals from fitting the
% tree-ring index to a arord model. 

% kw specifies whether ppt and temp series should be prewhitened 
% before correlation analysis.  kw=1: whiten    kw=0: do not whiten.

%****************  ARGUMENTS ********************************

% D - monthly ppt array. Each row a year t.  Each col a month.
% 	36 cols represent three years.   Col 1 is the year t.
% 		col 2   jan t-2
% 		col 3	feb t-2
%		col 13	dec t-2
%		col 14  jan t-1
%		col 26  jan t
%     col 37  dec t

% E - monthly temperature array

% C1 (cv) the tree-ring series, possibly prewhitened, depending

% L1 (cv) logical, specifying rows (years) of D  and E for analysis.
%	 entries 1 for years to be included in the analysis,
%	 and 0 for years not to be included.  If whitening (kw=1), the
%   first n elements will be lopped off, where n is the maximum 
%   ar order of the model fit to the ppt and temp series for the 
%   particular ending month and number of months.

% endmo - ending month of tree-ring growth year.  1=jan, 12=dec.
%  	So if want to gear analysis to growth year ending in august of
%	   year t, set endmo=8.

% arord - order of AR model previously fit to tree-rings.  If zero, the
%   corr analysis is done using the non-prewhitened tree-ring index.  

% kw:  1 prewhiten ppt and temp series
%      0 do not prewhiten temp and ppt series

% nhi:  max ar order to consider fitting to ppt, temp series
% 
%********************  CAUTIONS ******************************

% D and E should have been built so that cols 26-37 correspond to jan-dec of
%   year t.  
% D and E should not have missing values over the year range specified by
%   L1

%********************  OUTPUT   *********************************

% R matrix of partial correlations  tree.ppt vs temp.ppt.

% 	The ending month of the monthly grouping varies down the rows.
% 	The number of months in the grouping varies across the cols.

% 	Example:  assume endmo=8.

% 	Col 1 contains r between tree and one-month averages
%		col 1, row 1 is for August of year t
%		col 1, row 2 is for July   of year t
%     col 1,row 12 is for Sept   of year t-1
% 	Col 2 contains r between tree and two-month averages
%		col 2, row 1 is for July-Aug, year t
%		col 2, row 2 is for June-July, year t
% 		col 2,row 12 is for Aug t-1 thu Sept t-1
% 	Col 12 contains r between tree and 12-month averages
%		col 12, row  1 is for Sept t-1 thru Aug t
%		col 12, row  2 is for Aug  t1  thru July t
%    	col 12, row 12 is for Oct t-2 thru Sept t-1

% porder (12 x 12)  order of ar model fit to ppt series.  Monthly
%    groupings as in R.
% torder (12 x 12)  order of ar model fit to temp series.

% pratio (12 x 12) ratio of residual to original variance  for ar
%   models of ppt

% tratio (12 x 12) ratio of residual to original variance for ar 
%   models ot temp
%
% nn   sample size (valid # yr) for corrl mtx



%*********  PREALLOCATE *********************************************

R = zeros (12);  
porder =zeros(12,12);
torder=zeros(12,12);
pratio=ones(12,12);
tratio=ones(12,12);


C2 = C1;

S5=zeros(sum(L1),12);
S5T=zeros(sum(L1),12);

lastmo=37- (12-endmo);  % tie ending month to column of D

for k1=1:12;  % loop for each number of months in seasonal grouping
	top=lastmo-(k1-1):lastmo;  % rv.  Say [32 33] for k1 = 2,endmo=8
	S1 = top(ones(12,1),:); % dupe rows of top
	S2 = (0:11)';
	S3 = S2 (:,ones(k1,1)); % dupe cols of S2
	S4 = S1-S3;

	for k2 = 1:12;  % sum the clim variable over months
		if k1==1; % special case for vector
			S5(:,k2) = D(L1,S4(k2,:));
			S5T(:,k2)=E(L1,S4(k2,:));
		else 
			S5(:,k2) = (sum((D(L1,S4(k2,:)))'))';
			S5T(:,k2)= (sum((E(L1,S4(k2,:)))'))';
		end ;  % of if loop
	end;  % of for loop k2

	S5=dtrend(S5);
	S5T=dtrend(S5T);

	for k3=1:12;  % loop for ending month for given number of months
		
		if kw==0;  % do not prewhiten climate data
		
			W=S5(:,k3);  % a ppt time series
			V=S5T(:,k3); % a temperature time series
	  		if k==6;  % want partial corrs of tree.ppt vs temp.ppt
				x=[ones(length(C2),1)  W]\C2;
				y=[ones(length(C2),1)  W]\V;

				W1=[ones(length(C2),1) W] * x -C2;
				V1=[ones(length(C2),1)  W] * y -V;
	  		elseif k==5;  % want corrs between ppt and temp
				W1=W;
				V1=V;
			else
				error('CALL TO CORR3 WITH kw=0 REQUIRES K=5 OR K=6');
	  		end
			nn=length(C2);  % num of obs correls based on.			

		else;  % kw=1, whiten the climate data before analysis.

	   	[W,porder(k3,k1),pratio(k3,k1),arcs]=whit1(S5(:,k3),nhi,1);
			[V,porder(k3,k1),pratio(k3,k1),arcs]=whit1(S5T(:,k3),nhi,1);
			po=porder(k3,k1);  % ar order for this ppt series
			to= torder(k3,k1); % ar order for this temp series
			ord = max([po to]); %  maximum ar order for this ppt-temp pair
			nn=length(C2)-ord;


	  		if k==10 ;  % want partial corrs of tree.ppt vs temp.ppt
				x=[ones(nn,1)  W(ord+1:length(C2))]\C2(ord+1:length(C2));
				y=[ones(nn,1)  W(ord+1:length(C2))]\V(ord+1:length(C2));

				W1=[ones(nn,1) W(ord+1:length(C2))] * x -C2(ord+1:length(C2));
				V1=[ones(nn,1) W(ord+1:length(C2))] * y -V(ord+1:length(C2));

	  		elseif (k==9 | k==7 | k==8);  % want corrs between 
%				whitened ppt and whitened temp, between tree and
%				whitened ppt, or between tree and whitened temp.

				W1=W(ord+1:length(C2));
				V1=V(ord+1:length(C2));

			else
				error('CALL TO CORR3 WITH kw=1 REQUIRES K=7,8,9 or 10');
	  		end
      end;   % of for-loop for prewhitening (kw==1)

		if k~= 7 & k~= 8;  %  If correlations are between climate data only 
			r=corrcoef(W1,V1);
		elseif k==7;  % tree vs whitened ppt
			r=corrcoef(C2(ord+1:length(C2)),W1);
		elseif k==8;  % tree vs whitened temp
			r=corrcoef(C2(ord+1:length(C2)),V1);
		end
		R(k3,k1)=r(1,2);

	end;  % k3 for loop
	disp(k1)

end; % of loop for k1
 
