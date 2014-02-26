% Model of plate-based fluorescence titration binding assay to help optimize assay conditions for a range of Kd values.
% This script models digital titration by the Tecan HP D300 Digital Dispenser.
% 26 May 2012
% John D. Chodera <jchodera@berkeley.edu>

clear;

ligand_Kds = [10e-12, 100e-12, 1e-9, 10e-9, 100e-9, 1e-6, 10e-6, 100e-6, 1e-3]; % dissociation constants (M)
naffinities = length(ligand_Kds); % number of dissociation constants to plot binding curves for

%added_protein_concentration = 1e-6 / 2.5^2; % protein concentration added to plate (M)
added_protein_concentration = 150e-9; % protein concentration added to plate (M)
protein_volume = 130e-6; % protein volume added to plate (L)

%added_ligand_concentration = 20e-3; % 20 mM : ligand in DMSO concentration added to plate (M)
added_ligand_concentration = 2e-3; % 2 mM : ligand in DMSO concentration added to plate (M)
%added_ligand_volumes = 1.3e-6 * 10.^(-(5:-1:0)); % ligand volumes added to plate (L); note 13 pL minimum volume, 5 uL max volume, <= 8% CV for > 100 pL volume
nwells = 12;
added_ligand_volumes = logspace(log10(13e-12), log10(1.3e-6), nwells);

solubility_limit = 50e-6; % solubility limit of compound in final assay conditions(M)

nwells = length(added_ligand_volumes);
well_volumes = protein_volume + added_ligand_volumes;

% Compute total protein and ligand concentrations.
total_protein_concentrations = added_protein_concentration * (protein_volume ./ well_volumes);
total_ligand_concentrations = added_ligand_concentration * (added_ligand_volumes ./ well_volumes);

% Detection limit based on BioTek H4 specs and Nick's estimate of bosutinib fluorescence efficiency.
fluorescein_detection_limit = 15e-12 * 200e-6; % detection limit (mol)
detection_limit = 100 * fluorescein_detection_limit;

% Detection limit based on Nick's bosutinib:Abl paper.
% 1 nM in 1.5 mL cuvette.
detection_limit = 1e-9 * 1.5e-3; % detection limit (mol)

threshhold = 1e-12 * 100; % threshhold for detection of fluorescent species (M)

symbols = 'sdv^<>ph';
nsymbols = length(symbols);
markersize = 7;
%cmap = winter(naffinities); % for Nick
cmap = jet(naffinities);

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

% Compute total quantities.
molar_ratios = total_ligand_concentrations ./ total_protein_concentrations; % molar ratio of ligand:protein

% Compute concentrations of fluorescent species.
nrows = naffinities;
ncols = nwells;
complex_concentrations = zeros(nrows, ncols); % total concentration of fluorescent complex (M)
complex_quantities = zeros(nrows, ncols); % total quantity of fluorescent complex (mol)
for i = 1:nrows
  Kd = ligand_Kds(i);
  for j = 1:ncols
    complex_concentrations(i,j) = complex_concentration(total_protein_concentrations(j), total_ligand_concentrations(j), Kd);
    complex_quantities(i,j) = complex_concentrations(i,j) * well_volumes(j);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,1,1);

h = plot(1:nwells, complex_quantities', '.-');

% Use unique markers.
for i = 1:length(h)
  set(h(i), 'marker', symbols(mod(i,nsymbols)+1));
  set(h(i), 'markerfacecolor', cmap(i,:));
  set(h(i), 'markeredgecolor', cmap(i,:));
  set(h(i), 'color', cmap(i,:));
end

set(gca, 'XTick', 1:nwells);
xlabel('well number');
ylabel('complex quantity (mol)');
hold on;
set(h, 'MarkerSize', markersize);
plot([1 nwells], detection_limit * [1 1], 'k--'); % detection limit
text(5, detection_limit, 'detection limit')

axis([1 nwells 0 max(max(complex_quantities))]);

legend_text = {'1 mM', '100 uM', '10 uM', '1 uM', '100 nM', '10 nM', '1 nM', '100 pM', '10 pM'};
%legend(fliplr(legend_text), 'location', 'best');

title(sprintf('%.0f uL of %.3f nM protein, digital titration of 1.3 pL to 1.3 uL of %.0f mM ligand in DMSO', protein_volume / 1e-6, added_protein_concentration / 1e-9, added_ligand_concentration / 1e-3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,1,2);

h = semilogx(total_ligand_concentrations, complex_quantities', '.-');

% Use unique markers.
for i = 1:length(h)
  set(h(i), 'marker', symbols(mod(i,nsymbols)+1));
  set(h(i), 'markerfacecolor', cmap(i,:));
  set(h(i), 'markeredgecolor', cmap(i,:));
  set(h(i), 'color', cmap(i,:));
end

xlabel('total ligand concentration (M)');
ylabel('complex quantity (mol)');
hold on;
set(h, 'MarkerSize', markersize);
plot([total_ligand_concentrations(1) total_ligand_concentrations(end)], detection_limit * [1 1], 'k--'); % detection limit
oldaxis = axis;
axis([total_ligand_concentrations(1) total_ligand_concentrations(end) 0 oldaxis(4)]);
text(total_ligand_concentrations(5), detection_limit, 'detection limit')

% Mark solubility limit of compound.
plot(solubility_limit * [1 1], [0 oldaxis(4)], 'r-');

%legend_text = {'100 uM', '10 uM', '1 uM', '100 nM', '10 nM', '1 nM', '100 pM', '10 pM', '1 pM'};
%legend(fliplr(legend_text), 'location', 'northwest');

legend_text = {'1 mM', '100 uM', '10 uM', '1 uM', '100 nM', '10 nM', '1 nM', '100 pM', '10 pM'};
legend(fliplr(legend_text), 'location', 'west');

filename = 'dilution-series-model-hpd300.eps';
exportfig(gcf, filename, 'width', 10, 'height', 10, 'color', 'cmyk');
system(sprintf('epstopdf %s', filename));

