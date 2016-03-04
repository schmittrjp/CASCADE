function [dQ] = hydraulic_solver(h)
% Hydraulic solver in the case, that the d50 is already known. 
global kk wac q15 slp v kst rho s taucrit fminout d50

kk=kk+1;
% 1.: Strickler EQ with initial Strickler value
Rh=(h*wac)/(2*h+wac);

v=kst*(Rh)^(2/3)*slp^0.5;
% 2.:  Calculate initial Q=v*A
Q_calc=v*wac*h;
% 3.: initial d50 
% d50=(rho*h*slp)/((s-rho)*taucrit);

% 4.: estimate new kst
kst_analytic=29/(d50^(1/6)); % use the mean d50 of the reach

%5.: Calculate new flow velocity and discharge  

v1=kst_analytic*(Rh).^(2/3)*slp.^0.5;
Q_calc1=v1*wac*h;

% Objective function
dQ=abs(Q_calc1-q15);

fminout(kk,1)=h;
fminout(kk,2)=v1; 
fminout(kk,3)=d50;
fminout(kk,4)=dQ;

end

% 
% syms rho h s tcrit w h S Q kst
% 
% solve(Q==21/(rho*h*S/(1.6*tcrit))^(1/6)*(w*h/(w+2*h))^(2/3)*sqrt(S)*w*h,h)
% 
% 
% solve(Q/w/h==kst*(h*w/(2*h+w))^(2/3)*S^0.5,h)