%--- DOUBLE ROTATION

clear;
clc;
clf;

set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',14,'DefaultTextFontName','Times New Roman'); 

%path to dep profile
%sputter_profile = 'depline_Kaufman.mat';
sputter_profile = 'ExpData/depline_exp_130mm.mat';

C = 4.46; %thickness [nm] per minute
filename = 'depz.txt';
RELdeposition_coords_map_z = rot90(dlmread(filename, '', 1, 0),1);
row_dep = max(RELdeposition_coords_map_z(:));

%%%%INPUTS%%%%

%Choose source of get thickness data; 1 - seimtra, 0 - experiment
source = 0; 
val = 1;                    %1, 2, 3 - magnetron position
holder_inner_radius = 20;   % mm
holder_outer_radius = 145;  % radius sampleholder, mm

deposition_offset_x = -145; % mm
deposition_offset_y = -145; % mm
deposition_len_x = 290; % mm
deposition_len_y = 290; % mm
deposition_res_x = 1; % 1/mm
deposition_res_y = 1; % 1/mm

alpha0_sub = 0*pi;
R = 0;                     %radius of the planet orbit, mm               
NR = 5;                    %number of revolutions
omega = 3;                  %speed rev/min
a = 0:2*pi/360:NR*2*pi;     %angular position of the planet in its orbit
deltat = a(2)./(2*pi*omega);
process_time = NR/omega;

%Ns = 160;                  %numbers of teeth on the solar gears
%Np = 20;                   %numbers of teeth on the planet gears
%k = 1 + (Ns/Np);
k = 4.18;                    %if k%2 is integer then the period is 1 revolution. %k = 1 + (Ns/Np);

%%%%%%%

substrate_x_len = 100; % Substrate width, mm
substrate_y_len = 100; % Substrate length, mm

substrate_x_res = 10; % Substrate x resolution, mm
substrate_y_res = 10; % Substrate y resolution, mm

substrate_coords_x = R-substrate_x_len/2:substrate_x_len/substrate_x_res:R+substrate_x_len/2;
substrate_coords_y = -substrate_y_len/2:substrate_y_len/substrate_y_res:substrate_y_len/2;

if (R+sqrt(substrate_x_len^2+substrate_y_len^2)/2>holder_outer_radius)
    error ('Incorrect substate out of holder border.');
end

clf;

%%%% depoition profile meshing
deposition_coords_x = deposition_offset_x:deposition_res_x:deposition_offset_x+deposition_len_x-1;
deposition_coords_y = deposition_offset_y:deposition_res_y:deposition_offset_y+deposition_len_y-1;
[deposition_coords_map_x, deposition_coords_map_y] = meshgrid(deposition_coords_x, deposition_coords_y);

if source == 1
    deposition_coords_map_z = C * (RELdeposition_coords_map_z./row_dep);
    smooth_dep = smooth(deposition_coords_map_z(deposition_len_x/2,:),50, 'moving');
    max_point_y = max(smooth_dep);
    kk = find(smooth_dep == max_point_y);
    depliney = smooth_dep(kk:end);
elseif source == 0
    load(sputter_profile)
    if val == 1
        deposition_coords_map_z = dep_profile(deposition_coords_map_x, deposition_coords_map_y, -80, 59, C, 1);
    elseif val == 2
        deposition_coords_map_z = dep_profile(deposition_coords_map_x, deposition_coords_map_y, -80, -59, C, 1);
    elseif val == 3
        deposition_coords_map_z = dep_profile(deposition_coords_map_x, deposition_coords_map_y, 105.8, 0, C, 1);
    else
        error ('Incorrect magnetron position.');
    end
end

% % RING TARGET SPUTTERING PROFILE
% deposition_coords_map_z = 1*ring_sputter(deposition_coords_map_x, deposition_coords_map_y, 0, 0, 32, 40, 1, 30) + ... %X, Y, center_x, center_y, ring_r, ring_w, ring_i, substr_h
%                           0*ring_sputter(deposition_coords_map_x, deposition_coords_map_y, 78.25, -50, 35, 20, 1, 130) + ...
%                           0*ring_sputter(deposition_coords_map_x, deposition_coords_map_y, -78.25, -50, 35, 20, 1, 130);

                      
%GAUSS TARGET SPUTTERING PROFILE
%deposition_coords_map_z = gauss2d(deposition_coords_map_z, 200, [31  9])  + ...
%                          gauss2d(deposition_coords_map_z, 200, [15 31]) + ...
%                          gauss2d(deposition_coords_map_z, 200, [31 48]);

%%%% plot deposition area and holder circles
ang=0:0.01:2*pi;

holder_circle_inner_x=holder_inner_radius*cos(ang);
holder_circle_inner_y=holder_inner_radius*sin(ang);

holder_circle_outer_x=holder_outer_radius*cos(ang);
holder_circle_outer_y=holder_outer_radius*sin(ang);

%%%%
deposition_rect_x = [deposition_offset_x, deposition_offset_x+deposition_len_x, ...
    deposition_offset_x+deposition_len_x, deposition_offset_x, deposition_offset_x];
deposition_rect_y = [deposition_offset_y, deposition_offset_y, deposition_offset_y+deposition_len_y, ...
    deposition_offset_y+deposition_len_y, deposition_offset_y];

subplot(2,3,1);
hold on;
contourf(deposition_coords_map_x, deposition_coords_map_y, deposition_coords_map_z, 100, 'LineStyle', 'none');
colorbar;

%plot holder circle inner
plot(holder_circle_inner_x, holder_circle_inner_y, 'LineWidth', 2, 'Color', 'black', 'LineStyle', '--');
%plot holder circle outer
plot(holder_circle_outer_x, holder_circle_outer_y, 'LineWidth', 2, 'Color', 'black');
%plot drawing borders
plot(deposition_rect_x, deposition_rect_y, 'LineWidth', 2, 'Color', 'green');
 
xlabel('x, mm');
ylabel('y, mm');

plot1 = sprintf('Holder orientation, R = %d', R);
title(plot1);

%%%% PLOT SUBSTRATE

substrate_rect_x = [min(substrate_coords_x), max(substrate_coords_x), max(substrate_coords_x), min(substrate_coords_x), min(substrate_coords_x)];
substrate_rect_y = [max(substrate_coords_y), max(substrate_coords_y), min(substrate_coords_y), min(substrate_coords_y), max(substrate_coords_y)];
hold on;

plot(substrate_rect_x, substrate_rect_y, 'Color', 'black');

[substrate_coords_map_x, substrate_coords_map_y] = meshgrid(substrate_coords_x, substrate_coords_y);
scatter(reshape(substrate_coords_map_x, [numel(substrate_coords_map_x) 1]), ...
        reshape(substrate_coords_map_y, [numel(substrate_coords_map_y) 1]), 'x');


% The path of the point on the planet may be calculated using
xs = holder_outer_radius * cos(a); 
ys = holder_outer_radius * sin(a);
x = R*cos(a); 
y = R*sin(a);

%%%%%%
rho = zeros(size(substrate_coords_map_x));
alpha0 = zeros(size(substrate_coords_map_x));
xp = zeros([size(substrate_coords_map_x) numel(a)]);
yp = zeros([size(substrate_coords_map_x) numel(a)]);
z = zeros([size(substrate_coords_map_x) numel(a)]);
rel_th = zeros(1);
I = zeros([size(substrate_coords_map_x)]);

%gridded interpolant
F = griddedInterpolant(deposition_coords_map_x', deposition_coords_map_y', deposition_coords_map_z');

for i = 1:size(substrate_coords_map_x, 1)
    for j = 1:size(substrate_coords_map_x, 2)
        if substrate_coords_map_x(1, j) >= R
            rho(i, j) = sqrt(abs(substrate_coords_map_x(i, j)-R).^2 + substrate_coords_map_y(i, j).^2);
        else
            rho(i, j) = -sqrt(abs(substrate_coords_map_x(i, j)-R).^2 + substrate_coords_map_y(i, j).^2);
        end

        if i == round((substrate_y_res+1)/2) && j == round((substrate_x_res+1)/2)
            alpha0(i, j) = 0;
        else
            alpha0(i, j) = asin(substrate_coords_map_y(i, j) / rho(i, j));
        end
        %Drawing path
        for ii = 1:length(a)
            xp(i, j, ii) = R*cos(a(ii)+alpha0_sub)+rho(i, j)*cos(a(ii)*k + alpha0(i, j));
            yp(i, j, ii) = R*sin(a(ii)+alpha0_sub)+rho(i, j)*sin(a(ii)*k + alpha0(i, j));
            z(i, j, ii) =  F(xp(i, j, ii), yp(i, j, ii)); %get value z for xp, yp
            rel_th = 0;
        end
        for jj = 1:length(a)
            if jj > 1
                rel_th = rel_th + (z(i, j, jj-1) + z(i, j, jj))/2;
            else
                rel_th = rel_th + z(i, j, jj);
            end
        end
        I(i,j) = rel_th * deltat;
    end
end

% Drawing sampleholder
grid on
hold on 

h1 = plot(xs,ys,'LineWidth',0.5,'Color','b','LineStyle', '--');
axis([-deposition_len_x/2 deposition_len_x/2 -deposition_len_x/2 deposition_len_x/2]); 

subplot(2,3,2);
hold on
%contourf(deposition_coords_map_x, deposition_coords_map_y, deposition_coords_map_z, 100, 'LineStyle', 'none');
%colorbar;

%[xmin; ymin]
plot(squeeze(xp(substrate_x_res+1, round((substrate_x_res+1)/2), :)),squeeze(yp(substrate_y_res+1, round((substrate_y_res+1)/2), :)), 'LineWidth', 1, 'Color', 'r', 'LineStyle', '-')
%[xmax; ymax]
plot(squeeze(xp(1, round((substrate_x_res+1)/2), :)),squeeze(yp(1, round((substrate_y_res+1)/2), :)), 'LineWidth', 1, 'Color', 'r', 'LineStyle', '--')
%[xmax; ymin]
plot(squeeze(xp(round((substrate_x_res+1)/2), 1, :)),squeeze(yp(round((substrate_y_res+1)/2), 1, :)), 'LineWidth', 1, 'Color', 'g', 'LineStyle', '--')
%[xmin; ymax]
plot(squeeze(xp(round((substrate_x_res+1)/2), substrate_x_res+1, :)),squeeze(yp(round((substrate_y_res+1)/2), substrate_y_res+1, :)), 'LineWidth', 1, 'Color', 'c', 'LineStyle', '--')
%[xmax/2; ymax/2]
plot(squeeze(xp(round(substrate_x_res/2)+1, round(substrate_x_res/2)+1, :)),squeeze(yp(round(substrate_y_res/2)+1, round(substrate_y_res/2)+1, :)), 'LineWidth', 1, 'Color', 'k')
plot(substrate_rect_x, substrate_rect_y, 'Color', 'black');
hold off
grid on
axis([-deposition_len_x/2 deposition_len_x/2 -deposition_len_x/2 deposition_len_x/2]); 
xlabel('x, mm');
ylabel('y, mm');

RotPath = sprintf('Rotation path\nk = %3.2f, NR = %3.2f, omega = %3.2f rev/min time=%3.2f min', k, NR, omega, process_time);
%RotPath = sprintf('Rotation path k = %3.2f NR = %3.2f', k, NR);
title(RotPath);

subplot(2,3,3);
load(sputter_profile)
dep_norm = depliney./max(depliney);
plot(dep_norm, 'LineWidth', 1, 'Marker', 'x')
hold on
if val == 3
    pos3 = deposition_coords_map_z(145, 1:end);
    depsim_norm = pos3./max(pos3);
    plot(flip(depsim_norm), 'LineWidth', 1, 'Marker', 'x')
end
hold off
grid on

xlabel('x, mm');
ylabel('thickness, nm');
%title('Numerical profile. SIMTRA');
title('Sputtering profile. Time 40 min. 250W');

subplot(2,3,4);
contourf(substrate_coords_map_x, substrate_coords_map_y, I./max(I(:)), 100, 'LineStyle', 'none');
c = colorbar;
c.Label.String = 'Physical thickness, nm';

hold on

%plot line along x
plot([substrate_rect_x(1),substrate_rect_x(2)],[0,0], 'LineStyle' , '--', 'LineWidth', 2, 'Color', 'k')
%plot line along y
plot([R,R],[substrate_rect_y(4),substrate_rect_y(5)], 'LineStyle' , '--', 'LineWidth', 2, 'Color', 'r')
hold off

xlabel('x, mm');
ylabel('y, mm');

avg_thickness = mean(I(:));
avg_speed = avg_thickness/process_time;
%title(sprintf ('Relative thickness'));
title(sprintf ('Absolute thickness\n Mean = %3.2f Max = %3.2f Min = %3.2f\nAvg speed = %3.2f nm/min', avg_thickness, max(I(:)), min(I(:)), avg_speed));

subplot(2,3,5);
plot(substrate_coords_y,I(round(numel(substrate_coords_y)/2),:)./max(I(round(numel(substrate_coords_y)/2),:)), 'LineWidth', 2, 'Color', 'k', 'Marker', 'x')
grid on
xlabel('x, mm');
ylabel('y, rel unit');
title('Relative thickness along black line');

subplot(2,3,6);
plot(substrate_coords_y,I(:,round(numel(substrate_coords_y)/2))./max(I(:,round(numel(substrate_coords_y)/2))), 'LineWidth', 2, 'Color', 'r', 'Marker', 'x')
grid on
xlabel('x, mm');
ylabel('y, rel unit');
title('Relative thickness along red line');

figure
subplot(1,2,1)
G = griddedInterpolant(deposition_coords_map_x', deposition_coords_map_y', deposition_coords_map_z', 'linear', 'none');

rot_line_x = 0:holder_circle_outer_x/100:holder_circle_outer_x*(1-1/100);
rot_line_y = zeros(size(rot_line_x));
rot_line_z = zeros(size(rot_line_x));


for i = 1:numel(ang)
    % find nearest deposition point for current position of substrate point 
    rotated_x_map = rot_line_x.*cos(ang(i))-rot_line_y.*sin(ang(i));
    rotated_y_map = rot_line_x.*cos(ang(i))+rot_line_y.*sin(ang(i));
    
    rotated_z_map = G(rotated_x_map, rotated_y_map);
    
    if ~isempty(rotated_z_map(isnan(rotated_z_map)==1))
        error('Substrate point out of deposition area.');
    end
    
    rot_line_z = [rot_line_z + rotated_z_map];
end

ntt = 50;                                  % define how many angular division for the plot
theta = linspace(0,2*pi,ntt);              % create all the angular divisions
[rr,tt] = meshgrid(rot_line_x, theta);     % generate a grid 
z = repmat( rot_line_z, ntt , 1 );         % replicate our rot_line_z vector 50 times
[xx,yy,zz] = pol2cart(tt,rr,z);            % convert polar to cartesian coordinates
zzz = zz./max(zz(:));
contourf(xx,yy,zzz, 100, 'LineStyle', 'none');

colorbar;

hold on;

%plot(substrate_rect_x, substrate_rect_y, 'Color', 'black');
xlabel('x, mm');
ylabel('y, mm');
title('Single Rotation');

subplot(1,2,2)
len_y = size(zzz());
plot(rot_line_x, zzz(round(len_y(1)/2),:), 'LineWidth', 2)
grid on

name1 = 'Sample Thickness: ';
var1 = [name1, num2str(avg_thickness)];
disp(var1)

thick_witness = NR * deposition_coords_map_z(145, 145);
name2 = 'Witness Thickness: ';
var2 = [name2, num2str(thick_witness)];
disp(var2)

tooling_factor = thick_witness / avg_thickness;
name3 = 'Tooling Factor: ';
var3 = [name3, num2str(tooling_factor)];
disp(var3)
