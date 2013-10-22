load c:\projs\ac3\w1.dat
load c:\projs\ac3\w2.dat
load c:\projs\ac3\w3.dat
load c:\projs\ac3\w4.dat
load c:\projs\ac3\w5.dat
load c:\projs\ac3\w6.dat

[m,p]=size(w1);
a= -1.0;
w1(:,2) = w1(:,2) .*  a(ones(m,1),:);

[m,p]=size(w2);
clear a
a= -1.0;
w2(:,2) = w2(:,2) .*  a(ones(m,1),:);


[m,p]=size(w3);
clear a
a= -1.0;
w3(:,2) = w3(:,2) .*  a(ones(m,1),:);


[m,p]=size(w4);
clear a
a= -1.0;
w4(:,2) = w4(:,2) .*  a(ones(m,1),:);


[m,p]=size(w5);
clear a
a= -1.0;
w5(:,2) = w5(:,2) .*  a(ones(m,1),:);


[m,p]=size(w6);
clear a
a= -1.0;
w6(:,2) = w6(:,2) .*  a(ones(m,1),:);


plot(w4(:,1),w4(:,2),':',w1(:,1),w1(:,2),'--',w3(:,1),w3(:,2),'-')

text(32,-35,' --  D-13-14-25 CAD')
text(32,-40,' -   D-13-15-33 CBB')
text(32,-45,'...  D-14-15-02 DDA1')
xlabel('YEAR')
ylabel('DEPTH TO WATER (FT)')
