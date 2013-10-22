
k=0;
kill=0;   % will exit following while loop when kill=1

while (kill==0),
	mm=input('Span of Daniell Filter? ')	
	k=k+1;  % keep count of number of filters applied
	k1(k)=mm;   % store length of each successive filter
	
	%compute the new Daniell filter

		atemp=[.5 .5]';
	        wnew= conv(atemp,ones(mm-1,1)/(mm-1));
		% plot(w)
		% title('daniell filter')
		%pause

	% convolute Previous filter with new one, if not first in loop

		if  k==1,
			w=wnew;
		%	igo = n+1+(mm-1)/2;
		%	istop=2*n+(mm-1)/2;
			wsave=w;
		else
			w=wnew;
			shift=(mm-1)/2;
		%	igo=igo+shift ;
		%	istop=istop+shift ;
			wsave=conv(wsave,w);
		end
	
		plot(wsave)
		title('Daniel Filter Weights')
	kill=input('key 1 to kill loop, 0 otherwise : ');
end



