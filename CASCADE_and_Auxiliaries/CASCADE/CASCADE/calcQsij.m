function [Qsij] = calcQsij(dii,hjj,vjj,Sjj,Wjj,taucrit,t_prc)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%% The transport capacity for a grain size in a  downstream reach. 
            
            % hydraulic parameters of the reach, 
            
                if Sjj<10^-4
                    Sjj=10^-4;
                end
                if vjj<10^-5
                    vjj=10^-5;
                end

               
                    if dii<=2*10^(-3) % for sand                     
                       Qsij=Engelund_Hansen(hjj,dii,vjj,Wjj,Sjj,9.81)*t_prc; % sand transported  during the flow conditions of t_prc
                    else % for gravel
                        Qsij=Wong_Parker(hjj,dii,vjj,Wjj,Sjj,9.81,taucrit)*t_prc; % gravel transported  during the flow conditions of t_prc
                    end   
                                        
                    if isfinite(Qsij)==0 || isreal(Qsij)==0; Qsij=0; end
           
          
end

