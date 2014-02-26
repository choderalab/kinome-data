% Plot data.
clf;

clear;


ligand_concentration = 20 * 2.^(-(0:11));


wavelengths = { '280-350', '280-480', '350-480' };
nwavelengths = length(wavelengths);

nx = 2;
ny = nwavelengths;

for i = 1:nwavelengths
  wavelength = wavelengths{i};
  filename = sprintf('quinazoline-displacement-%s.txt', wavelength);
  x = textread(filename);
  x = x(1:8,1:12);
  
  subplot(ny,nx,nx*(i-1)+1);
  plot(x(1:4,:)', '.-');
  legend('erlotinib', 'bosutinib', 'gefitinib', 'buffer');
  title('imatinib titration');
  xlabel('ligand 1:2 dilution');
  ylabel(wavelength);
  set(gca,'XTick', 1:12);
  oldaxis = axis;
  axis([0 13 0  oldaxis(4)]);
  
  subplot(ny,nx,nx*(i-1)+2);
  plot(x(5:8,:)', '.-');
  legend('erlotinib', 'bosutinib', 'gefitinib', 'buffer');
  title('fasudil titration');
  xlabel('ligand 1:2 dilution');
  ylabel(wavelength);
  set(gca,'XTick', 1:12);
  oldaxis = axis;
  axis([0 13 0  oldaxis(4)]);
end

filename = 'quinazoline-displacement-dilution.eps';
%print('-depsc', filename);
scale = 0.7;
exportfig(gcf, filename, 'width', scale*7.5, 'height', scale*10, 'color', 'cmyk');
system(sprintf('epstopdf %s', filename));


