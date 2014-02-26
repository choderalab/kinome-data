% Plot data.
clf;

clear;

wavelengths = { '280-480' };
nwavelengths = length(wavelengths);

nreads = 4;

nx = 2;
ny = nwavelengths;

for i = 1:nwavelengths
  wavelength = wavelengths{i};

  % read data
  x = zeros(nreads,8,4);
  for j = 1:nreads
    filename = sprintf('quinazoline-injection-%d-%s.txt', j, wavelength);
    y = textread(filename);
    x(j,1:8,1:4) = y(1:8,1:4) ./ mean(y(:,5)); % normalize to quinazoline
  end
  
  subplot(ny,nx,nx*(i-1)+1);
  plot(0:3, squeeze(x(:,1:4,1)), '.-');
  hold on
  plot(0:3, squeeze(x(:,1:4,2)), '.-');
  plot(0:3, squeeze(x(:,1:4,3)), '.:');
  plot(0:3, squeeze(x(:,1:4,4)), '.:');
  %legend('erlotinib', 'bosutinib', 'gefitinib', 'buffer', 'location', 'northwest');
  title('20 uL injections of 40 uM');
  xlabel('injection');
  ylabel(sprintf('%s normalized to QZA', wavelength));
  set(gca,'XTick', 0:nreads);
  oldaxis = axis;
  axis([-0.5 (nreads-0.5) 0  oldaxis(4)]);
  
  subplot(ny,nx,nx*(i-1)+2);
  plot(0:3, squeeze(x(:,5:8,1)), '.-');
  hold on
  plot(0:3, squeeze(x(:,5:8,2)), '.-');
  plot(0:3, squeeze(x(:,5:8,3)), '.:');
  plot(0:3, squeeze(x(:,5:8,4)), '.:');
  legend('erlotinib', 'bosutinib', 'gefitinib', 'buffer', 'location', 'northwest');
  title('20 uL injections of 4 uM');
  xlabel('injection');
  ylabel(sprintf('%s normalized to QZA', wavelength));
  set(gca,'XTick', 0:nreads);
  oldaxis = axis;
  axis([-0.5 (nreads-0.5) 0  oldaxis(4)]);
end

filename = 'quinazoline-injections.eps';
%print('-depsc', filename);
scale = 0.6;
exportfig(gcf, filename, 'width', scale*7.5, 'height', scale*4.5, 'color', 'cmyk');
system(sprintf('epstopdf %s', filename));


