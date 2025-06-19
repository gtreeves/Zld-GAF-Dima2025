% script_GAF_figures
%

clear
close all

outfilename = '';
if ~isempty(outfilename) && ~strcmp(outfilename(end),'_')
	outfilename(end+1) = '_';
end

% {

%
% Load timecourses
%
load Mat\2024-01-27_14-09-13_GAF_timecourse
Soln1 = Soln;
load Mat/2024-03-08_00-00-20_GAF_sfGFP_tiles
Soln1 = [Soln1;Soln];
load Mat/2024-05-11_13-58-41_GAF_sfGFP_tiles-additional
Soln = [Soln1;Soln];


%
% Run the plotting function
%
data_GAF = extract_RICS_timecourse(Soln,true);

%}


%% ========================================================================
% Corrections for clusters
% =========================================================================
% {
struct2vars(data_GAF)

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
plot([0,300],S_factor*[0 300]+b,'k :','linewidth',2)


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
Qmean = mean(Q(:,end),'omitmissing');
Qmax = max(Q(:,end),[],'omitmissing');
errorbar(t+t0,Q,Q_e,'o-','linewidth',2)
hold on
plot(xlim,Qmean*[1 1],'k:')
plot(xlim,Qmax*[1 1],'k:')
set(gca,'fontsize',24,'xtick',-80:20:0)
xlim([-90 0])
ylim([0 0.2])
xlabel('time until gast. [min]')
ylabel('molecular brightness')



%}







%% ========================================================================
% Recalculate everything based on knowing that B and vareps are constant
% throughout 
% =========================================================================
% {



%
% Normalize brightness by some brightness value. We could eyeball this as
% the average brightness during nc 14, so we are not trying to correct for
% every small variation in Q. Really the only variation in Q we're worried
% about is the spike at the end of nc 13. So if the Q is less than the
% average of nc 14, we don't do anything.
%
Qmin = [max(Q(:,end),[],'omitmissing'),mean(Q(:,end),'omitmissing')];
Anuc_old = Anuc;
phinuc_old = phinuc;
Cnuc_old = Cnuc;

ii = 2; % ii = 1 means the normalization is to the max of nc14, 
% while ii = 2 means that the normalization is to the mean of nc14
q = Q/Qmin(ii); q(q < 1) = 1;

%
% Recalc variables
%
Anuc = Anuc./q;
phinuc = 1 - q.*(1 - phinuc_old);
Cnuc = 1./(V_eff*Anuc)/ 0.602;
Cfree = Cnuc.*(1 - phinuc);
CDNAbound = Cnuc.*phicc;
Cuncorrel = Cnuc.*(phinuc - phicc);
 
%}







%% ========================================================================
% Make plot of average curves
% =========================================================================
% {


%
% Make cell variables for each of the plots (X,Y,S)
%
NAN = NaN(size(Cfree));
NAN(idx_de) = 1;

X = [repmat({t+t0},1,4),Cfree.*NAN
	repmat({t+t0},1,4),Cfree.*NAN]; % the "1" is a placeholder since the final plot is different from the others
Y = {Cnuc, Cfree, CDNAbound, Cuncorrel, CDNAbound.*NAN
	Cnuc_ch, 1-phinuc, phicc, phinuc-phicc, Cuncorrel.*NAN};
Sx = {0, 0, 0, 0, Cfree_e.*NAN
	0, 0, 0, 0, Cfree_e.*NAN};
Sy = {Cnuc_e, Cfree_e, CDNAbound_e, Cuncorrel_e, CDNAbound_e.*NAN
	Cnuc_ch_e, phinuc_e, phicc_e, phiuncorrel_e, Cuncorrel_e.*NAN};

%
% Make graph-labeling variables
%
Xlabel = [repmat({'time until gast. [min]'},1,4),{'Free nuc conc'}
	repmat({'time until gast. [min]'},1,4),{'Free nuc conc'}];
Ylabel = {'Total nuc conc. [nM]', 'Free nuc conc. [nM]', 'Correlated conc. [nM]', 'Uncorrelated conc. [nM]', 'Correlated conc'
	'Histone conc [nM]', 'free (1-\phi_{2c,nuc})', '\phi_{cc}', 'uncorrel (\phi_{2c,nuc}-\phi_{cc})', 'Uncorrelated conc'};
YLIM2 = [inf(1,5)
	inf,1,1,1,inf];

%
% Make a plot of 10 tiles. 
%
figure('pos',[9         102        1899         592])
TL = tiledlayout(2,5,'TileSpacing','Compact','Padding','Compact');
C = colormap('lines');
for i = 1:2
	for j = 1:5
	nexttile % 1,1

	if j < 5
		errorbar(X{i,j},Y{i,j},Sy{i,j},'o-','linewidth',2)
	else
		errorbar(X{i,j},Y{i,j},Sy{i,j},Sy{i,j},Sx{i,j},Sx{i,j},'o-','linewidth',2)
		xlim([0 inf])
	end

	if j == 5
		Vmax = [35 9]; KD = [5 0.4];
		hold on
		x = linspace(0,4)';
		z = x/KD(i);
		plot(x,Vmax(i)*(z./(1 + z)))
	end

	xlabel(Xlabel{i,j})
	ylabel(Ylabel{i,j})
	ylim([0 YLIM2(i,j)])
	set(gca,'fontsize',12)

	end
end




	