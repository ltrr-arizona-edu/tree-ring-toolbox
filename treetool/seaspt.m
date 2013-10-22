function F = seaspt(A,begmo,endmo,k)
% seaspt:  seasonalize a monthly precipitation or temperature series
% F = seaspt(A,begmo,endmo,k);
% Last revised: 4-30-99
%
% Converts monthly climatic data to seasonal data.  Seasonal series can be
% summed over months (e.g., PPT), or averaged (e.g., PDSI). Seasons may cross the
% calendar-year boundary.  For example, a cool-season-total PPT might include
% the months from the previous Novemver through the current April.<P>
%
% The input monthly data is assumed to be stored in a Matlab matrix (13 col, with 
% year as col 1).  Be sure that missing data is coded as NaN in this matrix. That way,
% any seasonalized value relying on the missing data for any month will 
% result in a missing value (NaN) for the seasonal total
%
%*** INPUT 
%
% A(? x 13)r  monthly ppt or tmp; year in col 1; Jan data in col 2 ... Dec in col 13
% begmo (1 x 1)i Begining month of season. Jan=1, Dec=12.
% endmo (1 x 1)i Ending month of season.   Jan=1, Dec=12.
% k (1 x 1)i option for type of variable
%   k=1, A is ppt, runoff, effective ppt, or any variable summed over months  
%   k=2, A is temperature, pdsi, Z-index, or soil moisture
%   Only difference is
%		seasonalized ppt, etc,  is computed as sum over months, 
%     temperature, etc., is computed as average over months.
%
%***  OUTPUT ARGUMENT
% F (? x 2)  Seasonalized ppt or temp.  First col is year of ending month
%		of seasonal grouping.  Say begmo=11, endmo=1, and yrs=[1930 1950].
%		Then beginning year of seasonalized series (Nov-Jan) is 1931, which
%		is based on Nov and Dec of 1930 and Jan of 1931.
%
%
%*** REFERENCES -- none
%*** UW FUNCTIONS CALLED--none
%*** TOOLBOXES NEEDED -- none
%
%*** NOTES
%
% Maximum allowable number of months in season is 12.  So cannot 
% use seaspt.m to make, say, a 14-month season from June of last
% year through July of this year.  If you set begmo==6 and endmo==7,
% seaspt.m assumes you want a 2-month season, June+July for each year
% not a 14-month season.
%
% Revision 4-30-99 changed the input arguments to delete yrs, which
% specified the years of monthly climate data to use in forming the
% seasonal series.  Now all years of monthly data in A are used. The
% user can trim this monthly data before passing it to seaspt.m if he
% wants to limit the years of data used.
%
% The 4-30-99 revision also makes user assure that A is a 13-col matrix.
% Before, seaspt.m trimmed first 13 cols and user could have gotten away
% with passing trailing numeric columns -- such as an annual total, or a 
% numeric site id, in columns beyond column 13.
%
% The units of the seasonal data are the same as the units of the monthly 
% data. No conversion or scaling is done within seaspt.m

%-- Check sizes
[mA,nA]=size(A);
if nA~=13 | mA<2;
   error('A must be a 13-col matrix with at least 2 rows');
end

C=A; m1=mA; n1=nA; % renaming for consistency with earlier version


%****  Adjust beginning and ending years for seasonalized period to
%	    account for possible loss of initial year.  

if endmo < begmo;  % Ending month less than beginning month: must cross year boundary
	begmo = begmo-12; % For example, Nov (11) become (-1)
	m2go=1;   m2stop=m1-1;
	m3go=2;  m3stop=m1;
else
	m2go=1;  m2stop=m1;
	m3go=1;  m3stop=m1;
 end
nmos = endmo-begmo+1;  % number of months in "season"

D=[C(m2go:m2stop,:)  C(m3go:m3stop,2:13)];  % Has 25 cols
E=D(:,begmo+13:endmo+13);  % Subset of needed climatic data, not
%  including the year column.

% Compute the sum over months
if nmos ==1;  % special case to handle vector instead of matrix
	F(:,2)=E;
else;
   F(:,2) = (sum(E'))';  % Sum variable over months in season.
end

% If an "averaged" variable, like temperature, convert sum to average
if k==2;  
	F(:,2) = F(:,2) / nmos;
end

% Put on the year column
F(:,1) = C(m3go:m3stop,1);  % Appropriate years for seasonalized series.


