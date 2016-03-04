function [dQ] = hydraulic_solver(h)
%HYDRAULIC SOLVER: Solves open channel flow equations and calculates grain
% size given a bankfull discharge. This implementation makes some changes in the way the d50 is calculated 

global kk wac q15 slp v kst_analytic rho s taucrit fminout d90
kk=kk+1;

% 1.: h is a given input: calculate the resulting hydraulic radius 
Rh=(h*wac)/(2*h+wac);

% 2.: Estimate what grain size would result from flow depth and
% local morphology, using a fixed value for taucrit. 

d90=(rho*h*slp)/((s-rho)*taucrit); 

%3.: Calculate an analytic value for taucrit 
if h/d90<10  % for low water levels in comparison to the grain size, use the formula proposed by Suszka L (1991) Modification of transport rate formula for steep channels. Fluvial Hydraulics of Mountain Regions, Lecture Notes in Earth Sciences., eds Armanini PA, Silvio PGD (Springer Berlin Heidelberg), pp 59–70. Available at: http://link.springer.com/chapter/10.1007/BFb0011182 [Accessed March 2, 2015].
    taucrit=0.0851*(h/d90)^(-0.0261);
else  
    Re_d50=(1.6*9.81*d90)^0.5*d90/10^(-6); % particle Reynolds number
    taucrit=(0.22*Re_d50^(-0.6)+0.06*10^(-7.7*Re_d50^(-0.6)));
end  

if taucrit>1
asd=1;
end

% use the analytical value of taucrit to calculate a new equilibrium grain
% size. 
d90=(rho*h*slp)/((s-rho)*taucrit); 
        
%3.: Calculate an analytic value for the friction factor based on the grain size
kst_analytic=21.6/(d90^(1/6)); 

% 4.: Calculate a flow velocity
v=kst_analytic*Rh.^(2/3)*slp^0.5;

if v>20 
asd=0;
end 

% 5.: Estimate resulting discharge
Q_calc=v*wac*h;

% Objective function
dQ=abs(Q_calc-q15);

% fminout(kk,1)=h;
% fminout(kk,2)=v1; 
% fminout(kk,3)=d50;
% fminout(kk,4)=dQ;

end

% 
% syms rho h s tcrit w h S Q kst
% 
% solve(Q==21/(rho*h*S/(1.6*tcrit))^(1/6)*(w*h/(w+2*h))^(2/3)*sqrt(S)*w*h,h)
% 
% 
% solve(Q/w/h==kst*(h*w/(2*h+w))^(2/3)*S^0.5,h)