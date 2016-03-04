function [Qs] = Engelund_Hansen(h,d50,v,Wac,S,g)
% Calculates teh sediemnt transport according to the Engelund and Hansen Formula 
% Compare Czuba and Foufoulas (2014) and Engelund and Hansen (1967),see: http://repository.tudelft.nl/assets/uuid:81101b08-04b5-4082-9121-861949c336c9/Engelund_Hansen1967.pdf)
%%% INPUT: 
    %h: water depth 
    %d50: mean grain size 
    %Wac: Active channel width
    %v: velocity 
    %S: slope
    %g: gravitational acceleration
    %theta: scale factor: in which fraction of the water colum is the sediment
    %transport taking place?
    
%%% Output 
    % Qs: Sediment flux [kg/day]
       
rho_w=1000; % water density [kg/m^3]
rho_s=2600; % sediment densit [kg/m^3]
R=1.6; % Relative sediment density []

% Calculation of the Engelund/Hansen friction factor(Engelund/Hansen, p.28,3.1.3 ) 
    Cf=(2*g*S*h)./v.^2;
% Bed shear stress, tau
    tau=rho_w*g*h*S; 
% Dimensionless shear stress (see Flussbauskript, Bezzola, 2011, p. 5-8)
    tau_star=tau/(g*(rho_s-rho_w)*d50);

%%% Calculate sediment transport per unit width (m^2/s)
qs=(R*g*d50)^0.5*d50*0.05./Cf*tau_star^2.5;

%%% Calculate volumetric sediment transport(m^3/s), convert than to [kg/s]
    Qs=qs*Wac*rho_s;  
% Transport per day 
    Qs=Qs*60*60*24; 
  


end

