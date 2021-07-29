clear
clc
filename = '130mm_copy.txt';
data = dlmread(filename, '', 1, 0);
MagCenter = 47; 
x = (data(:,1)-MagCenter)';
v = (data(:,2))';
F = griddedInterpolant(x,v);
xq = linspace(0,300,300);
depliney = (F(xq))';

% hold on
% plot(x,v,'ro')
% hold on
% plot(xq,vq,'.')
% legend('Sample Points','Interpolated Values')

plot(xq, depliney)
profile_x_len = length(depliney);

save('depline_exp_130mm_copy.mat', 'depliney', 'profile_x_len')
