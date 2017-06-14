% Analytical solution to steady flow over a Wendand bump

clear all;
close all;

g = 9.81;

% Flow per unit width
q = 1.41;

% Bump maximum height
zb_max = 0.5;

% Bump half width
b = 150;

% Location of the bump maximum
x_max = 1000;

% Downstream water depth
h_B = .92;

% Upstream water depth
h_A = 1.32;

dx = 0.01;
x = (-500:dx:2000);
n = max(size(x));

for i = 1:n
    r = abs(x(i) - x_max)/b;
    zb(i) = 0;
    if (r >= 0) && (r<= 1)
            zb(i) = zb_max*(1 - r)^6*(1 + 6*r + (35/3.0)*r^2);
    else
        zb(i) = 0.0;
    end       
end



%subplot(2,1,1)
plot(x,zb,'-b')
axis([500 1500 0 2])
hold on

% subplot(2,1,2)
% plot(x,zb,'-b')
% axis([0 50 0 1])
hold on

% Calculate h_0 so that critical depth occurs over bump
h_max = 100;
h_min = 0.;
for i = 1:50
    h_0 = (h_max + h_min)/2.;
    u_0 = q/h_0;
    Fr_0 = u_0/sqrt(g*h_0);
    H_c = Fr_0*Fr_0/2. + 1. - 3./2.*Fr_0^(2./3.);
    if H_c*h_0 < zb_max
        h_min = h_0;
    else
        h_max = h_0;
    end
end
h_0
error = max(abs(h_0 - h_max),abs(h_0 - h_min))
u_0 = q/h_0
Fr_0 = u_0/sqrt(g*h_0)
E_T_0 = u_0*u_0/2./g + h_0

h_A
u_A = q/h_A
Fr_A = u_A/sqrt(g*h_A)
E_T_A = u_A*u_A/2./g + h_A

h_B
u_B = q/h_B
Fr_B = u_B/sqrt(g*h_B)
E_T_B = u_B*u_B/2./g + h_B

h_c = (q*q/g)^(1./3.)
u_c = q/h_c
Fr_c = u_c/sqrt(g*h_c)
E_T_c = u_c*u_c/2./g + h_c + zb_max

if h_A < h_0
% Calculate the downstream hydraulic jump location, toe velocity and depth
% and downstream conditions

    'critical flow is established based on obstacle height'
    
% Hydraulic jump only when downstream boundary condition is not
% supercritical
% calculate subsequent depth

    h_max = h_c;
    h_min = 0.;
    for i = 1:50
        h_s = (h_max + h_min)/2.;
        u_s = q/h_s;
        if u_s*u_s/2./g + h_s < E_T_0
            h_max = h_s;
        else
            h_min = h_s;
        end
    end
    h_s
    u_s
    Fr_s = u_s/sqrt(g*h_s)
    E_T_s = u_s*u_s/2./g + h_s
    error = max(abs(h_s - h_max),abs(h_s - h_min))
    
    if h_s > h_B
        
% Supercritical flow downstream of the obstacle (no jump)
        for j = 1:n
            if x(j) <= x_max
                h_max = 100.;
                h_min = zb(j);
                for i = 1:50
                    h_x = (h_max + h_min)/2.;
                    u_x = q/h_x;
                    if u_0*u_0*h_0*h_0/2./g/h_x/h_x + h_x + zb(j) < E_T_0
                        h_min = h_x;
                    else
                        h_max = h_x;
                    end 
                end
            else
                h_max = h_c;
                h_min = 0;
                zb(j)
                for i = 1:50
                    h_x = (h_max + h_min)/2.;
                    u_x = q/h_x;
                    if u_0*u_0*h_0*h_0/2./g/h_x/h_x + h_x + zb(j) >= E_T_0
                        h_min = h_x;
                    else
                        h_max = h_x;
                    end 
                end    
            end
%           h_x
%           error = max(abs(h_x - h_max),abs(h_x - h_min))

            h(j) = h_x + zb(j);
            u(j) = q/h_x;
        end
        
    else
        
% Subcritical flow downstream boundary (jump formed)
        h_max = 100;
        h_min = 0.;
        for i = 1:50
            h_toe = (h_max + h_min)/2.;
            u_toe = q/h_toe;
            zb_toe = E_T_c - (u_toe*u_toe/2./g + h_toe);
            h_j = (-h_toe + sqrt(h_toe*h_toe + 8./g*h_toe*u_toe*u_toe))/2.;
            u_j = q/h_j;
            E_T_j = u_j*u_j/2./g + h_j + zb_toe;
            if u_B*u_B/2./g + h_B < E_T_j
                h_max = h_toe;
            else
                h_min = h_toe;
            end
        end
        h_toe
        zb_toe
        u_toe
        Fr_toe = u_toe/sqrt(g*h_toe)
        x_toe = x_max + b*sqrt(1. - zb_toe/zb_max)
        error = max(abs(h_toe - h_max),abs(h_toe - h_min))
        
        h_j
        u_j
        Fr_j = u_j/sqrt(g*h_j)
        E_T_j = u_j*u_j/2./g + h_j
    
        for j = 1:n
            if x(j) <= x_max
                h_max = 100.;
                h_min = zb(j);
                for i = 1:50
                    h_x = (h_max + h_min)/2.;
                    u_x = q/h_x;
                    if u_x*u_x/2./g + h_x + zb(j) < E_T_0
                        h_min = h_x;
                    else
                        h_max = h_x;
                    end 
                end
            elseif x(j) < x_toe
                h_max = h_c;
                h_min = 0;
                for i = 1:50
                    h_x = (h_max + h_min)/2.;
                    u_x = q/h_x;
%                    if u_c*u_c*h_c*h_c/2./g/h_x/h_x + h_x + zb(j) < E_T_c
                    if u_x*u_x/2./g + h_x + zb(j) > E_T_0
                        h_min = h_x;
                    else
                        h_max = h_x;
                    end 
                end
            else
                u_j = q/h_j;
                E_T_j = u_j*u_j/2./g + h_j + zb_toe;
                h_max = 100.;
                h_min = zb(j);
                for i = 1:50
                    h_x = (h_max + h_min)/2.;
                    u_x = q/h_x;
                    if u_j*u_j*h_j*h_j/2./g/h_x/h_x + h_x + zb(j) < E_T_j
                        h_min = h_x;
                    else
                        h_max = h_x;
                    end 
                end
            
            end
%           hh
%           error = max(abs(h_x - h_max),abs(h_x - h_min))

            h(j) = h_x + zb(j);
            u(j) = q/h_x;
        end    
    end
   
   % subplot(2,1,1)
    plot(x,h,'-r')
    
   % subplot(2,1,2)
   % plot(x,u,'-r')
   % axis([0 50 0 4])


else

% Subcritical flow
    'sub-critical flow over obstacle'
    
    for j = 1:n
      h_max = h_0*2;
      h_min = zb(j);
      for i = 1:50
         h_x = (h_max + h_min)/2.;
         if u_0*u_0*h_0*h_0/2./g/h_x/h_x + h_x + zb(j) < E_T_A
             h_min = h_x;
         else
             h_max = h_x;
         end 
      end
      if zb(j) == zb_max
         h_x
         error = max(abs(h_x - h_max),abs(h_x - h_min))
         u_x = q/h_x
         Fr_x = u_x/sqrt(g*h_x)
      end
      h(j) = h_x + zb(j);
      u(j) = q/h_x;
    end
        
    %subplot(2,1,1)
    plot(x,h,'-r')
    
    %subplot(2,1,2)
    %plot(x,u,'-r')
    %axis([0 50 0 3])
    
    
end
plot(x,u,'-g')

fileID = fopen('xhub.txt','w');
A = [x;h;u;zb];
%fprintf(fileID,'%32s %32s %32s %32s\n','x','h','u','b');
formatSpec = '%4.28f,%4.28f,%4.28f,%4.28f\n';
fprintf(fileID,formatSpec,A);
fclose(fileID);