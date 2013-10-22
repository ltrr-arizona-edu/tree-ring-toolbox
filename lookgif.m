function lookgif
[file1,path1]=uigetfile('*.gif','Get file');
pf1=[path1 file1];
f=['''' file1 ''''];
figure(1);
eval(['p=imread('  f ');']);
imshow(p);
imzoom;
