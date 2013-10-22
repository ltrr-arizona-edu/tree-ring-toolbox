for i=1:nlags
	if i==1  ;   %  sum term zero for first-order r
		sum1 = 0;  
	else
		rsub=r(1:i-1); 
		rsubsq=rsub .* rsub;
		sum1=2*sum(rsubsq);
	end

	var=(1/m)  *  (1.0 + sum1);  % variance of ac coef
	se2(i) = 2.0 * sqrt(var);  % two standard errors
end

	
