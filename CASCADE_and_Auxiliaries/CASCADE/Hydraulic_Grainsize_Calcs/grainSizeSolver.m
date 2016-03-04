function [dQ] = hydraulic_solver(h)
%HYDRAULIC SOLVER: Solves open channel flow equations and calculates grain
% size given a bankfull discharge. This implementation makes some changes in the way the d50 is calculated 

global kk wac q15 slp v kst rho s taucrit fminout
kk=kk+1;

% 1.: h is a given input: calculate the resulting hydraulic radius 
Rh=(h*wac)/(2*h+wac);

% 2.: Estimate what grain size would result from flow depth and
% local morphology, using a fixed value for taucrit. 
d90=(rho*h*slp)/((s-rho)*taucrit); 

%3.: Calculate an analytic value for the friction factor based on the grain size
kst_analytic=21.6/(d90^(1/6)); 




% 4.: Calculate a flow velocity
v=kst_analytic*Rh.^(2/3)*slp^0.5;

% 5.: Estimate resulting discharge
Q_calc=v*wac*h;

v=kst*(Rh)^(2/3)*slp^0.5;

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