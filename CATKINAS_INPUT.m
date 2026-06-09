%(x): for note Energy input: [Ga(BEP)          |deltaG(scalling)  ]   %#ok
 (1): NH3  + # <-> NH3#      [ 1.00000  0.00000  0.18520  0.54830 ]   %#ok
 (2): NH3# + # <-> NH2# + H# [ 0.51940  1.16230  0.28470 -0.26480 ]   %#ok
 (3): NH2# + # <->  NH# + H# [ 0.75310  1.21310  0.48710  0.25830 ]   %#ok
 (4): NH#  + # <->   N# + H# [ 0.57900  1.14580  0.54470  0.20630 ]   %#ok
 (5): 2N#      <->  N2# +  # [ 0.24400  1.46070 -1.83220 -1.21170 ]   %#ok
 (6): 2H#      <->  H2# +  # [ 0.32410  0.46700 -0.21440  0.31810 ]   %#ok
 (7): N2#      <->   N2 +  # [ 0.00000  0.00000 -0.22280 -0.80860 ]   %#ok
 (8): H2#      <->   H2 +  # [ 0.00000  0.00000 -0.12500 -0.77790 ]   %#ok%Variation range(s) of descriptor(s):
E1 = linspace(-2.0,2.0,401);                                               % sampling energies for descriptor E1 (E_N(111)) 
%Reaction conditions:
T = 773;                                                                   % reaction temperature 
P_H2_FROZ  = 0.375;                                                        % H2 pressure (bar)
P_N2_FROZ  = 0.125;                                                        % N2 pressure (bar)
P_NH3_FROZ = 0.500;                                                        % NH3 pressure(bar)
Q_v_INIT = 1;                                                              % initial coverage of vacancy
%Control parameters:
CalcDRC = 1;                                                               % calculate the degree of rate control
ThermoMode = 0;                                                            % no correction            
BarrierMode = 1;                                                           % deal the adsorption with collision theory            
BEPMode = 1;                                                               % the energy input of barrier is by scalling 
Mr_NH3 = -2;                                                               % the barrier of the NH3 desorption is used by given                                                                                                                      
npar = 6;                                                                  % the process number for parallel computing                                              
PlotMode = [1:50:401];                                                     % plot every 50 points in [1,401] range 
