% ascppt1.m  : script file to get ascii files of
% original and estimated monthly ppt, san pedro basin




% Make files of the original data
fidxx = fopen ('namesx.txt','r'); % file names of input mat file
fidx = fopen('tryx.txt','r'); % corresponding full names
for n = 1:46
	disp(['Starting original series ',int2str(n)])
	fnm = fgetl(fidxx);  % file name for input
	eval(['load ',fnm]); % X will now hold the data
	txt = fgetl(fidx); % text to be centered
	digs = 2;  % digits to convert data to integer
	fn = [fnm '.dat'];
	fpfmon1(X,fn,txt,digs)
	eval(['clear '  fn])
end
fclose('all')

% Make files of the original data, with estimates for missing
fid11 = fopen ('names1.txt','r'); % file names of input mat file
fid1 = fopen('try1.txt','r'); % corresponding full names
for n = 1:35
	disp(['Starting estimate series ',int2str(n)])
	fnm = fgetl(fid11);  % file name for input
	eval(['load ',fnm]); % Y1 will now hold the data
	txt = fgetl(fid1); % text to be centered
	digs = 2;  % digits to convert data to integer
	fn = [fnm '.dat'];
	fpfmon1(Y1,fn,txt,digs)
	eval(['clear '  fn])
end
fclose('all')
