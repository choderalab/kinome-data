% Model of plate-based fluorescence binding assay.

Kd = 1e-6; % 1 uM (worst case)
Kd = 250e-9; % 250 nM
%Kd = 1e-9; % 1 nM

protein_concentrations = 10.^(-[6:10]);
%protein_concentrations = 2.^(-[0:7]);

threshhold = 1e-12 * 100;

dilution_factor = 2;
ndilutions = 7;
ligand_concentrations = 10e-6 * dilution_factor.^(-[0:ndilutions]);
%ligand_concentrations = 2e-6:-1e-7:0;

% Solve for complex concentration given total ligand and macromolecule.
% Kd = P*L / PL
% P + PL = Pt
% L + PL = Lt
% Kd = (Pt - PL)*(Lt - PL) / PL
% PL Kd = Pt Lt - PL Pt - PL Lt + PL^2
% PL^2 - (Pt + Lt + Kd) PL + Pt Lt = 0
% PL = ((Pt + Lt + Kd) +- sqrt((Pt + Lt + Kd)^2 - 4 Pt Lt)) / 2

complex_concentration = @(Pt,Lt) ((Pt + Lt + Kd) - sqrt((Pt + Lt + Kd).^2 - 4*Pt*Lt)) / 2;

nrows = length(protein_concentrations);
ncols = length(ligand_concentrations);

C = zeros(nrows, ncols);
for i = 1:nrows
  for j = 1:ncols
    C(i,j) = complex_concentration(protein_concentrations(i), ligand_concentrations(j));
  end
end

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
