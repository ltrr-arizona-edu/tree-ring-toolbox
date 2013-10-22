function lookjpg
[file1,path1]=uigetfile('*.jpg','Get file');
pf1=[path1 file1];
f=['''' file1 ''''];
figure(1);
eval(['p=imread('  f ');']);
imshow(p);
imzoom;