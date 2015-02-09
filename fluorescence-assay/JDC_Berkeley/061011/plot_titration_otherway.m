% Plot data.
clf;

clear;


ligand_concentration = 20 * 2.^(-(0:11));


wavelengths = { '280-350', '280-480', '350-480', '280-350-480-FRET', '290-420', '300-480', '270-420' };
nwavelengths = length(wavelengths);

nx = 4;
ny = nwavelengths;

for i = 1:nwavelengths
  wavelength = wavelengths{i};
  filename = sprintf('quinazoline-titration-%s.txt', wavelength);
  x = textread(filename);
  x = x(:,1:12);
  
  subplot(ny,nx,nx*(i-1)+1);
  plot(ligand_concentration', x(1:2,:)', '.-');
  title('erlotinib');
  if i == 1
    legend('1 uM p38', 'buffer');
  end
  xlabel('ligand (uM)');
  ylabel(wavelength);
  %set(gca,'XTick', 1:12);
  oldaxis = axis;
  %axis([0 13 0  oldaxis(4)]);
  %set(gca, 'XTickLabel', {'20 uM', '10 uM', '5 uM', '2.5 uM', '1.25 uM', '625 nM', '313 nM', '156 nM', '78 nM', '39 nM', '20 nM', '10 nM'})
  
  subplot(ny,nx,nx*(i-1)+2);
  plot(ligand_concentration', x(3:4,:)', '.-');
  %legend('1 uM p38', 'buffer');
  title('bosutinib');
  xlabel('ligand (uM)');
  ylabel(wavelength);
  %set(gca,'XTick', 1:12);
  oldaxis = axis;
  %axis([0 13 0  oldaxis(4)]);

  subplot(ny,nx,nx*(i-1)+3);
  plot(ligand_concentration', x(5:6,:)', '.-');
  %legend('1 uM p38', 'buffer');
  title('gefitinib');
  xlabel('ligand (uM)');
  ylabel(wavelength);
  %set(gca,'XTick', 1:12);
  oldaxis = axis;
  %axis([0 13 0  oldaxis(4)]);
  xlabel('ligand (uM)');

  subplot(ny,nx,nx*(i-1)+4);
  plot(ligand_concentration', x(7:8,:)', '.-');
  %legend('1 uM p38', 'buffer');
  title('quinazoline');
  xlabel('ligand (uM)');
  ylabel(wavelength);
  %set(gca,'XTick', 1:12);
  oldaxis = axis;
  %axis([0 13 0  oldaxis(4)]);
end

filename = 'quinazoline-titration-concentrations.eps';
%print('-depsc', filename);
scale = 1.5;
exportfig(gcf, filename, 'width', scale*7.5, 'height', scale*10, 'color', 'cmyk');
system(sprintf('epstopdf %s', filename));

