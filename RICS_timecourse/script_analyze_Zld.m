% script_analyze_Zld
%
% This script loads data from excel file given in the "filename" variable
% just below. You will need to specify the "basepth", which points to the
% location of the excel file on your computer. You will also have to
% specify which sheet you wish to analyze (Sheet1 or "tiles-additional").
%
% On a separate note, to use the excel file, you will have to edit Column A
% on either sheet to point to the folder where you stored the image (czi)
% files. Be sure not to edit anything else about the excel file (unless you
% know what you're doing)...the function "timecourse_parser" expects the
% sheet to be in a certain structure.


clear
close all
yesplot = true;


%
% Run the time course parser so we have the attributes of each file and
% folder that we want to run (stages, filenames, start and end indices)
%
basepth = 'G:\Shared drives\Reeveslab\Trainees\Sadia Siddika Dima\Confocal\sfGFP-Zld-H2ARFP\'; % needs to be edited
filename = [basepth 'Tiles.xlsx'];
% sheetname = 'Sheet1';
sheetname = 'tiles-additional';
genotype = 'Zld-sfGFP'; % optional
outfilename = 'Zld_sfGFP'; % optional
if ~isempty(outfilename) && ~strcmp(outfilename(end),'_')
	outfilename(end+1) = '_';
end

%
% Run the time course parser so we have the attributes of each file and
% folder that we want to run (stages, filenames, start and end indices)
%
[filenames,stages,istarts,iends,embryo_indices,basenames] = ...
	timecourse_parser(filename,sheetname);


%
% Geometric params
%
w0 = 0.25;
w0r = 0.3;
chi = 3;
wz = chi*w0;

%
% RICS hyperparams
%
subtr_bg = true;% subtract background within ftn_imread (calls ftn_calcbg)
stabilize = true; % how many frames are avg'd together in each group? floor(H/ngroups)
usemask = 10;
sws = 5;
ntmax = 100; % comfortably max number of time points

NC_list = 10:14;
clk1 = clock;
clk = [num2str(clk1(1)),'-',num2strDU(clk1(2),2),'-',num2strDU(clk1(3),2),'_',...
	num2strDU(clk1(4),2),'-',num2strDU(clk1(5),2),'-',num2strDU(round(clk1(6)),2)];

%
% Run the for loop
%
Soln = [];
n_embryos = max(embryo_indices);
Data = cell(n_embryos,1);
LE = cell(n_embryos,1);

%%
for iii = 1:n_embryos
	try
		v = find(embryo_indices == iii);

		data = run_analyze_RICS_t(filenames(v),stages(v),istarts(v,:),iends(v,:),...
			ones(size(v)),basenames(v),[],genotype,w0,wz,sws,[],usemask,outfilename,[],[],yesplot);
		Data{iii} = data;
	catch lastE
		LE{iii} = lastE;
		fprintf('%s in %s\n',lastE.message,filename)
	end
end

save(['Mat/',clk,'_',outfilename(1:end-1)],'Data','LE')

for iii = 1:n_embryos
	Soln = [Soln;Data{iii}];
end

save(['Mat/',clk,'_',outfilename(1:end-1)],'Soln','LE')












