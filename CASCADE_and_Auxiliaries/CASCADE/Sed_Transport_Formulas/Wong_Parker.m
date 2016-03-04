function [Qs] = Wong_Parker(h,d50,v,Wac,S,g,taucrit,I)
%gravel_velocity calculates the velocity scaling according to Cuba &
% Foufoulas (2014), Eq.: 25 for sediemnt larger 2mm

%%% Inputs:
    %h: water depth 
    %d50: mean grain size 
    %Wac: Active channel width
    %v: velocity 
    %S: slope
    %g: gravitational acceleration
    %alpha, beta: from the MPM equation, see Wong & Parker 2006
    %Lai: active transport layer

%%% Output 
    % Qs: Sediment flux [kg/day]

% parameters from Wong&Parker 2006, EQ. 24
alpha=3.97;
beta=1.5; 

rho_w=1000; % water density [kg/m^3]
rho_s=2600; % sediment densit [kg/m^3]
R=1.6; % Relative sediment density []

% Bed shear stress, tau
    tau=rho_w*g*h*S; 
% Dimensionless shear stress (see Flussbauskript, Bezzola, 2011, p. 5-8)
    tau_star=tau/(g*(rho_s-rho_w)*d50);

%%% Calculate sediment transport per unit width (m^2/s)
qs=(R*g*d50).^0.5*d50*alpha*(tau_star-taucrit)^beta;     
    
%%% Calculate volumetric sediment transport(m^3/s), convert than to [kg/s]
    Qs=qs*Wac*rho_s;  
    
% Convert from  [kg/s] to per year [kg/d]
    Qs=Qs*60*60*24; 


if isreal(Qs)==0
    Qs=inf;
end

end

