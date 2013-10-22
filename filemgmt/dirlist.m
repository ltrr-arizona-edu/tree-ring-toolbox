function dirlist(suffix)
% dirlist(suffix)
%
% Make listing of files in directory.  If no input argument, list files with any
% suffix.  Otherwise only those with specified suffix.

eval(['dir *.' suffix ' >filnms.out'])