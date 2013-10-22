for i=1:250
	plot(D(i,1:20),R2(i,:),'*')
	
	title(['SITE ',int2str(i)])
	pause
end
