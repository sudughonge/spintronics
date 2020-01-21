clc;
clear all;
filename = input("Enter the name of the file\n", "s");
A = dlmread(filename);
T = A(:, 2);
T = T.*3.483;
[f, x] = hist(T, 100);
figure(1)
bar(x, f/trapz(x, f));
ax = gca();
set(ax, 'fontsize', 20);
hx = xlabel('t (picoseconds)');
hy = ylabel('P(t)');
set(hx, 'fontsize', 20);
set(hy, 'fontsize', 20);
print(strcat('plot_', filename, '.png'), "-dpng");

