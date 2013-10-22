function infout=hydacc1(P,T,PE,awcs,awcu,snowinf,kopt,ssi,su1)
% hydacc1:  hydrologic accounting to be used in PDSI and other indices
% CALL: infout=hydacc1(P,T,PE,awcs,awcu,snowinf,ssi,su1);
%
% Meko 5-24-97
%
%******************************  IN ***************************
%
% P (mP x 12)r monthly ppt (in) ; columns represent jan,feb,...,dec
% T (mT x 12)r monthly mean T (oF)
% PE (mPE x 12) monthly mean PE (in)
% awcs (1 x 1)r available water capacity in surface soil layer (in)
% awcu (1 x 1)r likewise for underlying layer
%
% snowinf{} input relating to snow storage model
%		Tcrit (1 x 1)r critical monthly mean temperature (oF) below which
%			ppt is assumed to be snow.  Typically 30.2 oF.  May set to []
%			if kopt(1)==1.
%		mons1 (1 x 1)i start, end months of 'snow period', typically [10 5], 
%			meaning Oct-May.  This is the group of months snow is allowed to
%			fall in.  Then any snow that accumulated must be totally melted
%			by the start of the next snow year.  For the [10 5] example,
%			snow must be melted by the end of september. User flagged with
%			error if snow falls later than mons(2) month. May set to []
%			if kopt(1)==1;
%		melttbl (? x ?)r  lookup table to be developed relating number of
%			months for the melt to total snow accumulation in cold preceding
%			months and T of the first month with T>=Tcrit. May set to [] if
%			kopt(1)==1;
%		snowgo (1 x 1)r Jan 1 initial starting value of snow accumulation.
%			May set to [] if kopt(1)==1; or if kopt(1)==2 and this a pass-1 run;
%
% kopt -- options
%		kopt(1) handling of snow vs rain
%			==1 ignore snow, treat all ppt as rain
%			==2 treat snow by redistributing monthly ppt input to soil moisture
%
% <ssi> (1 x 1)r initial soil moisture (in) for Jan 1 of first year for
%		surface soil layer
% <ssu> (1 x 1)r likewise for underlying soil layer
%
%
%***********************  OUT ***********************************
%
% infout{} output data from accounting -- same size tsm's as T and P. Much
% of this info follows the notation of Palmer.
%
%		{1} DELSS change in soil moisture for surface layer
%		{2} DELSU change in soil moisture for undelying layer
%		{3} DELS change in soil moisture in both layers
%		{4} SSP starting ('prime', in notation of Palmer) soil moisture, surface
%		{5} SUP likewise for underlying
%		{6} SP  likewise for both layers combined
%		{7} SS end of month soil moisture in surface layer
%		{8} SU end of month soil moisture in underlying layer
%		{9} S end of month soil moisture in both layers
%		{10} HP start of month snow storage
%		{11} H end of month snow storage
%		{12} PP revised ppt, adjusted for snow storage addition or subtraction
%		{13} R recharge
%		{14} PR potential Recharge
%		{15} RO runoff
%		{16} L loss
%		{17} PL potential loss
%		{18} ET evapotranspiration
%
%
%
%******************* NOTES *************************************
%
% ssi and ssu are optional -- either real data or empty.  If empty, 
% you want a 'pass-1" run which will come up with reasonable initial
% values by running a balance for a 35 year synthetic 'time series'
% of long-term monthly normals
%
% Snow storage model redistributes monthly ppt for those years with
% snow accumulation.  Snow is assumed if the monthly T <Tcrit. Accumulated
% snow is assumed to melt completely in the first  1-several months with
% T>=Tcrit.  The number of months (nmelt) is interpolated from a 2-d
% lookup table based on the total snow accumulation for the year and 
% the temperature for the month following the last month with 
% T<Tcrit.  So far, this is a constant lookup table, because I 
% have not developed an empirical relation ship. I assume a 2-month melting
% period.
%
% The percentage of accumulated snow assigned as redistributed ppt to the
% first warm months is proportional to T-Tcrit for those months.  For
% example, if the first two warm months have T = 32.2 and 40.2, the
% weights are 2/10 and 8/10
%
% Stepwise strategy. 
% (1) In calling function, compute monthly normals of P, T, PE; and
%		make synthetic 35-year time series of those normals repeated 35 rows 
% (2) call hydacc1.m with these sythetic series and with ssi sui and snowgo as [].
%		These [] settings tell the function to automatically assume the soil
%		is saturated and snow accum is zero at the start of Jan of the first
%		year.  By the end of the 35-year accounting, more realistic values
%		for these parameters will be arrived at.
% (3) Call hydacc1.m with the true historical monthly time series data
%		of T, P, PE, and with the newly arived at values for ssi,sui, snowgo
%
%
%*************************  FUNCTIONS CALLED
%
% highrow.m -- utility function to find highest row with nonzero logical element
%		in each col of a matrix
% snowppt.m -- utility function to redistribute monthly ppt over months by
%		adjusting for storage of ppt in a "snow buffer". Required inputs are
%		same-size matrices T and P.  snowppt.m returns the revised P.  Function
%		snowppt.m is not yet written.  Is needed only if a monthly T is below
%		Tcrit, and if kopt(1)==2.
% 

% Thats where i'm at so far.
