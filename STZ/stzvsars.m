function stzvsars
% stzvsars:  stz versus arsstan: comparison of standard chrons by programs stz and arstan
% stzvsars;
% Last revised 1-29-01
%
% User has made a standardized site chronology using the sequence of function in the
% STZ suite, and wants to compare time series features with those of a chronology by ratio method in 
% ARSTAN with curve option 2:  neg exp, SL neg slope or HL
%
%*** INPUT
%
% No arguments
% User prompted to point to the file <arstans.dat> with the arstan chronology, and the .mat storage 
% file with the STZ version of the chronology
%
%
%*** OUTPUT
%
% No args
% Time series plot, box plot
%
%*** REFERENCES -- NONE
%*** UW FFUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED ---NONE
%
%*** NOTES

[file1,path1]=uigetfile('*.mat','Infile storing STZ-version standard chronology (e.g., frywt1.mat)');
pf1=[path1 file1];

[file2,path2]=uigetfile('arstans.dat','Infile with standard chron from pgm ARSTAN');
pf2=[path2 file2]


strtit=['Comparison, standard chrons in ' file1 ' and ' file2];

eval(['load ' pf1 ' yrZI ZI;']);
yr1=yrZI;
x1=ZI(:,1);

eval(['load ' pf2 ';']);
yr2=arstans(:,1);
x2=arstans(:,2);

figure (1);
plot(yr1,x1,yr2,x2);
xlabel('Year');
ylabel('Index');
legend('STZ','ARSTAN');
title(strtit);
grid;
zoom xon;

figure (2);
plot(yr1,zscore(x1),yr2,zscore(x2));
xlabel('Year');
ylabel('zscore of index');
legend('STZ','ARSTAN');
title(['Z-score ' strtit]);
grid;
zoom xon;

figure(3);
[m1,n1]=size(x1);
[m2,n2]=size(x2);
if m1<m2;
    d=m2-m1;
    y1=[x1 ; repmat(NaN,d,1)];
    y2=x2;
    
elseif m2<m1;
    d=m1-m2;
    y2=[x2 ; repmat(NaN,d,1)];
else;
    y1=x1;
    y2=x2;
end;
T=[y1 y2];
boxplot(T);
set(gca,'XTickLabel',{'STZ','ARS'});
title('Boxplot of indices');
    




