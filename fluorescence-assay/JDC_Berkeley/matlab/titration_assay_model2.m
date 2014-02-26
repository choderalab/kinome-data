% Model of plate-based fluorescence titration binding assay to help optimize assay conditions for a range of Kd values.
% 24 Apr 2012
% John D. Chodera <jchodera@berkeley.edu>

clear;

ligand_Kds = [10e-12, 100e-12, 1e-9, 10e-9, 100e-9, 1e-6, 10e-6, 100e-6, 1e-3]; % dissociation constants (M)
naffinities = length(ligand_Kds); % number of dissociation constants to plot binding curves for

protein_concentration = (1e-6)/2.5^3; % protein concentration in plate (M)
protein_volume = 100e-6; % protein volume in plate (L)

% Detection limit based on BioTek H4 specs and Nick's estimate of bosutinib fluorescence efficiency.
fluorescein_detection_limit = 15e-12 * 200e-6; % detection limit (mol)
detection_limit = 100 * fluorescein_detection_limit;

% Detection limit based on Nick's bosutinib:Abl paper.
% 1 nM in 1.5 mL cuvette.
%detection_limit = 1e-9 * 1.5e-3; % detection limit (mol)

titrant_concentrations = [10, 10, 10, 10, 10, 10, 10, 10, 10, 10] * 1e-6; % fluorescent ligand concentrations (M)
titrant_volumes        = [10, 10, 10, 10, 10, 10, 10, 10, 10, 10] * 1e-6; % fluorescent ligand titrant volumes (L)

%titrant_concentrations = [10, 10, 10, 10, 10, 100, 100, 100, 100, 100, 500, 500, 500, 500, 500] * 1e-6; % fluorescent ligand concentrations (M)
%titrant_volumes        = [10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10] * 1e-6; % fluorescent ligand titrant volumes (L)

%titrant_concentrations = [2, 2, 2, 2, 2, 20, 20, 20, 20, 20, 200, 200, 200, 200, 200] * 1e-6; % fluorescent ligand concentrations (M)
%titrant_volumes        = [10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10] * 1e-6; % fluorescent ligand titrant volumes (L)

%titrant_concentrations = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096] * 1e-6; % fluorescent ligand concentrations (M)
%titrant_volumes        = [10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10] * 1e-6; % fluorescent ligand titrant volumes (L)

%titrant_concentrations = 2.^(-(11:-1:0)) * 10e-6; % fluorescent ligand concentrations (M)
%titrant_volumes        = 10 * ones(size(titrant_concentrations)) * 1e-6; % fluorescent ligand titrant volumes (L)

titrant_concentrations = 2.^(-(7:-1:0)) * 10e-6; % fluorescent ligand concentrations (M)
titrant_volumes        = 10 * ones(size(titrant_concentrations)) * 1e-6; % fluorescent ligand titrant volumes (L)

%N = 1;
%titrant_concentrations = [ 1*ones(1,N) 10*ones(1,N) 100*ones(1,N) 1000*ones(1,N) 10000*ones(1,N) 100000*ones(1,N) 1000000*ones(1,N)] * 10e-9; % fluorescent ligand concentrations (M)
%titrant_volumes        = [10*ones(1,N) 10*ones(1,N)  10*ones(1,N)   10*ones(1,N)    10*ones(1,N)     10*ones(1,N)      10*ones(1,N)] * 1e-6; % fluorescent ligand titrant volumes (L)

%titrant_concentrations = [1e-9, 10e-9, 100e-9, 1e-6, 10e-6, 100e-6, 1e-3, 10e-3, 100e-3]; % fluorescent ligand concentrations (M)
%titrant_volumes        = 10 * 1e-6 * ones(size(titrant_concentrations)); % fluorescent ligand titrant volumes (L)

% 10x dilution series [1 nM to 100 uM] followed by titration with 100 uM; 10 uL injection volumes
%titrant_concentrations = [1e-9, 10e-9, 100e-9, 1e-6, 10e-6, 100e-6, 100e-6, 100e-6, 100e-6, 100e-6]; % fluorescent ligand concentrations (M)
%titrant_volumes        = 10 * 1e-6 * ones(size(titrant_concentrations)); % fluorescent ligand titrant volumes (L)

% 10x dilution series [1 nM to 10 uM] followed by titration with 10 uM; 10 uL injection volumes
%titrant_concentrations = [100e-9, 100e-9, 100e-9, 100e-9, 100e-9, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 10e-6, 10e-6, 10e-6, 10e-6, 10e-6]; % fluorescent ligand concentrations (M)
%titrant_volumes        = 10 * 1e-6 * ones(size(titrant_concentrations)); % fluorescent ligand titrant volumes (L)

%titrant_concentrations = [1*ones(1,6) 100*ones(1,6) 10000*ones(1,6) 1000000*ones(1,6)] * 1e-9; % fluorescent ligand concentrations (M)
%titrant_volumes        = [1, 2, 4, 8, 16, 32, 1, 2, 4, 8, 16, 32, 1, 2, 4, 8, 16, 32, 1, 2, 4, 8, 16, 32] * 1e-6; % fluorescent ligand titrant volumes (L)

%titrant_concentrations = [1*ones(1,4) 100*ones(1,4) 10000*ones(1,4) 1000000*ones(1,4)] * 1e-9; % fluorescent ligand concentrations (M)
%titrant_volumes        = [1, 5, 10, 20, 1, 5, 10, 20, 1, 5, 10, 20, 1, 5, 10, 20] * 1e-6; % fluorescent ligand titrant volumes (L)

%titrant_concentrations = [10*ones(1,3) 100*ones(1,3) 1000*ones(1,3) 10000*ones(1,3)] * 1e-9; % fluorescent ligand concentrations (M)
%titrant_volumes        = [5, 10, 20, 5, 10, 20, 5, 10, 20, 5, 10, 20] * 1e-6; % fluorescent ligand titrant volumes (L)

%titrant_concentrations = [10*ones(1,3) 100*ones(1,3) 1000*ones(1,3) 10000*ones(1,3)] * 1e-9; % fluorescent ligand concentrations (M)
%titrant_volumes        = [5, 10, 20, 5, 10, 20, 5, 10, 20, 5, 10, 20] * 1e-6; % fluorescent ligand titrant volumes (L)

%titrant_concentrations = [10*ones(1,2) 100*ones(1,2) 1000*ones(1,2) 10000*ones(1,2) 100000*ones(1,2)] * 1e-9; % fluorescent ligand concentrations (M)
%titrant_volumes        = [10, 20, 10, 20, 10, 20, 10, 20, 10, 20] * 1e-6; % fluorescent ligand titrant volumes (L)

%titrant_concentrations = [10*ones(1,2) 100*ones(1,2) 1000*ones(1,2) 10000*ones(1,2) 100000*ones(1,2), 1000000*ones(1,2)] * 1e-9; % fluorescent ligand concentrations (M)
%titrant_volumes        = [10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10] * 1e-6; % fluorescent ligand titrant volumes (L)

%titrant_concentrations = [0.2, 0.2, 0.2, 0.2, 0.2, 2, 2, 2, 2, 2, 20, 20, 20, 20, 20, 200, 200, 200, 200, 200] * 1e-6; % fluorescent ligand concentrations (M)
%titrant_volumes        = [10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10] * 1e-6; % fluorescent ligand titrant volumes (L)

%titrant_concentrations = [10, 10, 10, 10, 5, 5, 5, 5, 5, 100, 100, 100, 100, 100, 100, 1000, 1000, 1000, 1000, 1000] * 1e-6; % fluorescent ligand concentrations (M)
%titrant_volumes        = [10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10] * 1e-6; % fluorescent ligand titrant volumes (L)

%titrant_concentrations = [45, 1, 1, 1, 1, 1, 10, 10, 10, 100, 100, 100, 1000, 1000, 1000] * 1e-6; % fluorescent ligand concentrations (M)
%titrant_volumes        = [10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10] * 1e-6; % fluorescent ligand titrant volumes (L)

threshhold = 1e-12 * 100; % threshhold for detection of fluorescent species (M)

ninjections = length(titrant_volumes); % number of injections

% Solve for complex concentration given total ligand and protein concentrations.
% Kd = P*L / PL
% P + PL = Pt
% L + PL = Lt
% Kd = (Pt - PL)*(Lt - PL) / PL
% PL Kd = Pt Lt - PL Pt - PL Lt + PL^2
% PL^2 - (Pt + Lt + Kd) PL + Pt Lt = 0
% PL = ((Pt + Lt + Kd) +- sqrt((Pt + Lt + Kd)^2 - 4 Pt Lt)) / 2

% Pt and Lt are in M, result is in M.
complex_concentration = @(Pt,Lt,Kd) ((Pt + Lt + Kd) - sqrt((Pt + Lt + Kd).^2 - 4*Pt*Lt)) / 2;

% Compute cumulative volumes for different injections.
cumulative_volumes = [protein_volume];
for injection = 1:ninjections
  cumulative_volumes = [cumulative_volumes, cumulative_volumes(end)+titrant_volumes(injection)];
end

% Compute total quantities.
total_protein_quantities = [protein_volume*protein_concentration];
total_ligand_quantities = [0];
for injection = 1:ninjections
  total_protein_quantities = [total_protein_quantities, total_protein_quantities(end)]; % stays the same
  total_ligand_quantities = [total_ligand_quantities, total_ligand_quantities(end)+titrant_volumes(injection)*titrant_concentrations(injection)];
end
molar_ratios = total_ligand_quantities ./ total_protein_quantities; % molar ratio of ligand:protein

% Compute concentrations as a function of injection.
total_ligand_concentrations = total_ligand_quantities ./ cumulative_volumes;
total_protein_concentrations = total_protein_quantities ./ cumulative_volumes;

% Compute concentrations of fluorescent species.
nrows = naffinities;
ncols = ninjections+1;
complex_concentrations = zeros(nrows, ncols); % total concentration of fluorescent complex (M)
complex_quantities = zeros(nrows, ncols); % total quantity of fluorescent complex (mol)
for i = 1:nrows
  Kd = ligand_Kds(i);
  for j = 1:ncols
    complex_concentrations(i,j) = complex_concentration(total_protein_concentrations(j), total_ligand_concentrations(j), Kd);
    complex_quantities(i,j) = complex_concentrations(i,j) * cumulative_volumes(j);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(4,1,1);

titrant_volumes_in_uL = titrant_volumes / 1e-6; % uL
h = bar(1:ninjections, titrant_volumes_in_uL);
axis([0 ninjections 0 1.1*max(titrant_volumes_in_uL)]);

xlabel('injection number');
ylabel('injection volume (uL)');

title(sprintf('initial well contents: %.0f uL of %.3f nM protein', protein_volume / 1e-6, protein_concentration / 1e-9));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(4,1,2);

%h = bar(1:ninjections, titrant_concentrations);
%set(gca, 'yscale', 'log');
h = semilogy(1:ninjections, titrant_concentrations, 'b.');
axis([0 ninjections min(titrant_concentrations)*0.4 max(titrant_concentrations)*2]);
set(gca, 'ytick', [1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e-0]);
set(h, 'MarkerSize', 20);

xlabel('injection number');
ylabel('injection concentrations (M)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%subplot(5,1,3);
%
%h = plot([0:ninjections], 1e6*cumulative_volumes, '.-');
%set(h, 'MarkerSize', 20);
%
%xlabel('injection number');
%ylabel('well volume (uL)');
%hold on;
%max_safe_volume = 250e-6;
%max_volume = 360e-6;
%plot([0 ninjections], 1e6*max_safe_volume*[1 1], 'g-');
%plot([0 ninjections], 1e6*max_volume*[1 1], 'r-');
%axis([0 ninjections 0 1e6*max_volume]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(4,1,3);

h = plot(0:ninjections, complex_quantities', '.-');

xlabel('injection number');
ylabel('fluorescent species (mol)');
hold on;
set(h, 'MarkerSize', 20);
plot([0 ninjections], detection_limit * [1 1], 'k--'); % detection limit
text(5, detection_limit, 'detection limit')

axis([0 ninjections 0 max(max(complex_quantities))]);

legend_text = {'1 mM', '100 uM', '10 uM', '1 uM', '100 nM', '10 nM', '1 nM', '100 pM', '10 pM'};
%legend(fliplr(legend_text), 'location', 'best');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(4,1,4);

%h = plot(molar_ratios, complex_quantities', '.-');
h = semilogx(molar_ratios, complex_quantities', '.-');

xlabel('molar ratio (ligand/protein)');
ylabel('fluorescent species (mol)');
hold on;
set(h, 'MarkerSize', 20);
plot([molar_ratios(2) molar_ratios(end)], detection_limit * [1 1], 'k--'); % detection limit
oldaxis = axis;
axis([molar_ratios(2) molar_ratios(end) 0 oldaxis(4)]);
text(molar_ratios(5), detection_limit, 'detection limit')

%legend_text = {'100 uM', '10 uM', '1 uM', '100 nM', '10 nM', '1 nM', '100 pM', '10 pM', '1 pM'};
%legend(fliplr(legend_text), 'location', 'northwest');

legend_text = {'1 mM', '100 uM', '10 uM', '1 uM', '100 nM', '10 nM', '1 nM', '100 pM', '10 pM'};
legend(fliplr(legend_text), 'location', 'best');

%filename = 'titration-model-standard.eps';
%filename = 'titration-powers-of-two.eps';
%filename = 'titration-powers-of-ten.eps';
%filename = 'titration-powers-of-ten-solubility-limited.eps';
filename = 'titration-test.eps';
exportfig(gcf, filename, 'width', 10, 'height', 10, 'color', 'cmyk');
system(sprintf('epstopdf %s', filename));

return

C
imagesc(C);
colorbar;

% Plot assay.
clf;
subplot(2,1,1);
hold on;
Cmax = max(max(C));
npoints = 100;
theta = linspace(0,2*pi,npoints);
r = 0.3;
xc = r.*cos(theta);
yc = r.*sin(theta);
for i = 1:nrows
  for j = 1:ncols
    x0 = j;
    y0 = i;
    color = [1 1 1] - (C(i,j)/Cmax)*[0 1 1];
    h = fill(x0+xc, y0+yc, color);
    set(h, 'edgecolor', [1 1 1]);
  end
end
set(gca,'YDir', 'reverse');
axis equal;
axis([0.5 ncols-0.5 0.5 nrows-0.5]);
set(gca,'XTick', 1:ncols);
set(gca,'YTick', 1:nrows);
xlabel('ligand dilution');
ylabel('kinase dilution');
title('K_d = 250 nM');

subplot(2,2,3);
plot(ligand_concentrations / threshhold, 'k-', 'linewidth', 2);
hold on;
plot(C' / threshhold);
ylabel('relative fluorescence over threshhold');
xlabel('ligand dilution');
axis([0 ncols+1 0 max(max(C / threshhold))]);
legend('max signal', '1 uM kinase', '100 nM kinase', '10 nM kinase', '1 nM kinase');
title('assuming 100 pM detection threshhold');

subplot(2,2,4);
plot(log10(ligand_concentrations / threshhold), 'k-', 'linewidth', 2);
hold on;
plot(log10(C' / threshhold));
ylabel('log10 rel. fluorescence over threshhold');
xlabel('ligand dilution');
axis([0 ncols+1 0 max(max(log10(ligand_concentrations / threshhold)))]);
legend('max signal', '1 uM kinase', '100 nM kinase', '10 nM kinase', '1 nM kinase');
title('assuming 100 pM detection threshhold');

%exportfig(gcf, '250nM.eps', 'width', 10, 'height', 7.5, 'color', 'cmyk');
