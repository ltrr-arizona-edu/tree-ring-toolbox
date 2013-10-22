% code for plotting factor scores of factors 1 and 2

clg
subplot(211)
yr=(1700:1979)';

m1=length(yr);
ze=zeros(m1,1);

plot(yr,F(:,1),yr,ze);
text(1800,0.5,'Factor 1');
text(1705,2.0,'74 %');

plot(yr,F(:,2),yr,ze);
text(1800,0.5,'Factor 2');
text(1705,2.0,'10 %');
