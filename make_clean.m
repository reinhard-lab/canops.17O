%
% #########################################################################
% ### careful, now...this will delete all existing output/figures #########
% #########################################################################
%
% navigate to home directory
cd ~/Documents/MATLAB/canops.17O/;
%
% remove any data in workspace
clear all;
%
% remove any existing output
if exist('~/Documents/MATLAB/canops.17O/_output/','dir')
    rmdir('~/Documents/MATLAB/canops.17O/_output/','s')
end
%
% remove any existing figures
if exist('~/Documents/MATLAB/canops.17O/_figures/','dir')
    rmdir('~/Documents/MATLAB/canops.17O/_figures/','s')
end
%
% -------------------------------------------------------------------------