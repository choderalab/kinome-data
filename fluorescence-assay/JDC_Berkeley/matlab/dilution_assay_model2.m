% Model of plate A8 plate-based fluorescence titration binding assay to help optimize assay conditions for a range of Kd values.
% 24 Apr 2012
% John D. Chodera <jchodera@berkeley.edu>

clear;

ligand_Kds = [10e-12, 100e-12, 1e-9, 10e-9, 100e-9, 1e-6, 10e-6, 100e-6, 1e-3]; % dissociation constants (M)
naffinities = length(ligand_Kds); % number of dissociation constants to plot binding curves for

added_protein_concentration = 160e-9; % protein concentration added to plate (M)
protein_volume = 100e-6; % protein volume added to plate (L)

%added_ligand_concentrations = 10e-6 * 2.^(-(11:-1:0)); % ligand concentrations added to plate (M)
added_ligand_concentrations = 20e-6 * 2.^(-(7:-1:0)); % ligand concentrations added to plate (M)
ligand_volume = 100e-6; % ligand volume added to plate (L)

nwells = length(added_ligand_concentrations);
well_volume = protein_volume + ligand_volume;

% Compute total protein and ligand concentrations.
total_protein_concentrations = added_protein_concentration * (protein_volume / well_volume) * ones(1,nwells);
total_ligand_concentrations = added_ligand_concentrations * (ligand_volume / well_volume);

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
cmap = winter(naffinities);

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
    complex_quantities(i,j) = complex_concentrations(i,j) * well_volume;
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

xlabel('well number');
ylabel('complex quantity (mol)');
hold on;
set(h, 'MarkerSize', markersize);
plot([1 nwells], detection_limit * [1 1], 'k--'); % detection limit
text(5, detection_limit, 'reliable detection limit')

axis([1 nwells 0 max(max(complex_quantities))]);

legend_text = {'1 mM', '100 uM', '10 uM', '1 uM', '100 nM', '10 nM', '1 nM', '100 pM', '10 pM'};
%legend(fliplr(legend_text), 'location', 'best');

title(sprintf('%.0f uL of %.3f uM protein + %.0f uL ligand at variable concentration', protein_volume / 1e-6, added_protein_concentration / 1e-6, ligand_volume / 1e-6));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data from plate A8 (280,480) gain 125

observed_mixture_fluorescence = [28823   27300   27387   26747   27165   26085   21672   10045   5581    3279    2505    1873];
observed_ligand_fluorescence = [3410    3031    1807    1305    1161    979     1135    973     1042    1020    1063    993];
observed_buffer_fluorescence = [1010    974     1042    1002    1059    989     1017    961     1067    1040    1007    1001];
observed_plate_fluorescence = [705     748     773     758     743     868     781     884     887     825     836     791];

hold on;
x = observed_mixture_fluorescence - observed_ligand_fluorescence - observed_buffer_fluorescence;
x = (x - min(x)) / (max(x) - min(x));
shift = 0.06e-11;
x = x * (max(max(complex_quantities)) - shift) + shift;
plot(1:12, fliplr(x), 'k.-', 'markersize', markersize*4);

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
text(total_ligand_concentrations(5), detection_limit, 'reliable detection limit')

%legend_text = {'100 uM', '10 uM', '1 uM', '100 nM', '10 nM', '1 nM', '100 pM', '10 pM', '1 pM'};
%legend(fliplr(legend_text), 'location', 'northwest');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data from plate A8 (280,480) gain 125

observed_mixture_fluorescence = [28823   27300   27387   26747   27165   26085   21672   10045   5581    3279    2505    1873];
observed_ligand_fluorescence = [3410    3031    1807    1305    1161    979     1135    973     1042    1020    1063    993];
observed_buffer_fluorescence = [1010    974     1042    1002    1059    989     1017    961     1067    1040    1007    1001];
observed_plate_fluorescence = [705     748     773     758     743     868     781     884     887     825     836     791];

hold on;
x = observed_mixture_fluorescence - observed_ligand_fluorescence - observed_buffer_fluorescence;
x = (x - min(x)) / (max(x) - min(x));
shift = 0.06e-11;
x = x * (max(max(complex_quantities)) - shift) + shift;
plot(total_ligand_concentrations, fliplr(x), 'k.-', 'markersize', markersize*4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

legend_text = {'observed', 'limit', '1 mM', '100 uM', '10 uM', '1 uM', '100 nM', '10 nM', '1 nM', '100 pM', '10 pM'};
legend(fliplr(legend_text), 'location', 'west');

filename = 'dilution-series-model.eps';
exportfig(gcf, filename, 'width', 10, 'height', 10, 'color', 'cmyk');
system(sprintf('epstopdf %s', filename));


