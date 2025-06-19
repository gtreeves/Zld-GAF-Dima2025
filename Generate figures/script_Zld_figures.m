% script_Zld_figures
%
% This script extracts the data from the analyzed Zld files, stored in raw
% form in mat files, and averages the data across multiple embryos.
%
% The second half of the script runs the brightness analysis

clear
close all

% {

outfilename = '';
if ~isempty(outfilename) && ~strcmp(outfilename(end),'_')
	outfilename(end+1) = '_';
end

%
% Load timecourses
%
load Mat\2023-10-15_00-49-57_Zld_timecourse
Soln1 = Soln;
load Mat/2024-05-12_08-45-38_Zld_sfGFP_tiles
Soln1 = [Soln1;Soln];
load Mat/2024-05-12_01-21-23_Zld_sfGFP_tiles-additional
Soln = [Soln1;Soln];

%
% Filter because some embryos went out of focus
%
badembryos = [10,11,14,18,19];
Soln(badembryos) = [];

%
% Extract and average the data, which also creates a plot
%
data_zld = extract_RICS_timecourse(Soln,true);


%}

%% ========================================================================
% Corrections for clusters
% =========================================================================
% {
%
struct2vars(data_zld)

%
% Calculate the S-factor
%
[S_factor,b,stats] = totallsq(Inuc(:),shotnoise(:),Inuc_e(:),shotnoise_e(:));
disp(['S = ',num2str(S_factor)])
disp(['b = ',num2str(b)])
disp(stats)

figure
errorbar(Inuc,shotnoise,shotnoise_e,shotnoise_e,Inuc_e,Inuc_e,'o-','linewidth',2)
hold on
plot([0,500],S_factor*[0 500]+b,'k :','linewidth',2)


%
% Plot apparent brightness, molecular brightness
%
figure
errorbar(t+t0,B,B_e,'o-','linewidth',2)
hold on
plot([-90 0],S_factor*[1 1],'k:')
set(gca,'fontsize',24,'xtick',-80:20:0)
xlim([-90 0])
ylim([0 500])
xlabel('time until gast. [min]')
ylabel('apparent brightness')

figure
Q = B/S_factor - 1;
Q_e = B_e/S_factor;
errorbar(t+t0,Q,Q_e,'o-','linewidth',2)
set(gca,'fontsize',24,'xtick',-80:20:0)
xlim([-90 0])
ylim([0 0.2])
xlabel('time until gast. [min]')
ylabel('molecular brightness')



