% Check dilution series accuracy by using both absorbance and fluorescence.

clear;

% from "2012-05-14 plate B1 quick absorbance and fluorescence.txt", row B

absorbance = [0.152   0.101   0.081   0.067   0.062   0.057   0.057   0.056   0.055   0.055   0.055   0.053]; % absorbance at 280 nm
fluorescence = [3348    2223    1134    606     422     349     316     298     264     240     268     271]; % excite 280 nm / emission 480 nm (top)

wells = 1:12;

% Fit absorbance and fluorescence.
well_volume = 75e-6; % L
well_area = 0.1586 / 100 / 100; % m^2
path_length = (well_volume * 1000 / (100^3)) / well_area ; % m

%fluorescence_model = @(concentration,epsilon,I0Q,path_length,buffer_fluorescence) I0Q*epsilon.*path_length.*concentration + buffer_fluorescence; % linear
fluorescence_model = @(concentration,epsilon,I0Q,path_length,buffer_fluorescence) I0Q*(1-exp(-epsilon.*path_length.*concentration)) + buffer_fluorescence; % primary inner filter effect 
%fluorescence_model = @(concentration,epsilon,I0Q,path_length,buffer_fluorescence) I0Q*(1-10.^(-path_length.*absorbance)) ./ (path_length*absorbance).*concentration + buffer_fluorescence; % primary inner filter effect (alternative correction)
%fluorescence_model = @(concentration,epsilon,I0Q,path_length,buffer_fluorescence) I0Q*10.^(-absorbance).* (1-exp(-epsilon.*path_length.*concentration)) + buffer_fluorescence; % primary inner filter effect (alternative correction)
%fluorescence_model = @(concentration,epsilon,I0Q,path_length,buffer_fluorescence) I0Q*10.^(-absorbance).*(epsilon.*path_length.*concentration) + buffer_fluorescence; % primary inner filter effect (alternative correction)

absorbance_model = @(concentration,epsilon,epsilon_buffer) concentration*epsilon*path_length + epsilon_buffer*path_length;

concentration = 20.0e-6 * 2.^(0:-1:-11); % M
%objective = @(epsilon,I0Q,buffer_fluorescence,epsilon_buffer) sum( ((fluorescence - fluorescence_model(concentration,epsilon,I0Q,path_length,buffer_fluorescence)) ./ fluorescence).^2 ) + sum( ((absorbance - absorbance_model(concentration,epsilon,epsilon_buffer)) ./ absorbance).^2 ); % relative to individual measurement values
objective = @(epsilon,I0Q,buffer_fluorescence,epsilon_buffer) sum( ((fluorescence - fluorescence_model(concentration,epsilon,I0Q,path_length,buffer_fluorescence)) ./ fluorescence(1) ).^2 ) + sum( ((absorbance - absorbance_model(concentration,epsilon,epsilon_buffer)) ./ absorbance(1)).^2 ); % relative to largest measurement value

buffer_fluorescence = fluorescence(end);
epsilon_buffer = absorbance(end) / path_length;
epsilon = (absorbance(1) - absorbance(end)) / (path_length * concentration(1));
I0Q = (fluorescence(1) - fluorescence(end)) / (epsilon * path_length * concentration(1)) / 10.^(-absorbance(1));

x = log([epsilon, I0Q, buffer_fluorescence, epsilon_buffer]);
x = fminsearch(@(x) objective(exp(x(1)), exp(x(2)), exp(x(3)), exp(x(4))), x);
epsilon = exp(x(1));
I0Q = exp(x(2));
buffer_fluorescence = exp(x(3));
epsilon_buffer = exp(x(4));

expected_fluorescence = fluorescence_model(concentration, epsilon, I0Q, path_length, buffer_fluorescence);
expected_absorbance = absorbance_model(concentration, epsilon, epsilon_buffer);

clf;

subplot(2,1,1);

plot(wells, fluorescence, 'o', wells, expected_fluorescence, 'x-');
legend('observed', 'expected');
xlabel('dilution well');
ylabel('280 excite / 480 emission (top)');

subplot(2,1,2);

plot(wells, absorbance, 'o', wells, expected_absorbance, 'x-');
legend('observed', 'expected');
xlabel('dilution well');
ylabel('280 absorbance');

% Write
print -depsc dilution-series-check.eps
system('epstopdf dilution-series-check.eps');
