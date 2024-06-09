% Currents that modulate spike shape dynamics on multiple timescales induce
% ramping bursts in model respiratory neurons
% Code written by : Dr. Victor Matveev, NJIT

function Cost = BurstModel_dspk(Params, plotFlag )

Params   = abs(Params);  % Physical parameters can't be negative
gNaP     = Params(1); 
gsyn     = Params(2);
gK       = Params(3);
gNa      = Params(4);
Vh       = Params(5);
sh       = Params(6);
Vht      = Params(7);
sht      = Params(8);
EK       = Params(9);
ENa      = Params(10);
VhNaP    = Params(11);
khNaP    = Params(12);
VhNaPt   = Params(13);
khtNaP   = Params(14);

fprintf('Pars = [%g %g %g %g %g %g %g %g %g] ', Params);

% Conductances (nS) 
gL  =4;
thNaPmax=5250;

%INa and INaP Time Constants (ms) 
tmNamax = 0.25; 
thNamax  = 8.46;
tmNaPmax = 1;   
thNa2=1010;

%Constants
C = 36; 

%Sodium Current Constants (mV)
VmNa  = -43.8; VhNa = -68; VmNat  = -43.8; VhNat = -67.5;
kmNa  = 6;     khNa = -11.9;
kmtNa = 14;   khtNa = -12.8;

%Persistent Sodium Current
VmNaP  = -47.1; 
VmNaPt  = -47.1; 
kmNaP  = 3.1;   
kmtNaP = 6.2;   

%Constants for n
nA = 0.011; nAV = 44; nAk = 5;
nB = 0.17;  nBV = 49; nBk = 40;

%Reversal Potentials (mV)
EL  = -62.5;  Esyn = -10;

%Currents
INa  = @(m,h,h2,V)    gNa  * m.^3 .* h .* h2 .* (V - ENa);
INaP = @(m,h,V)    gNaP * m    .* h .* (V - ENa);
IK   = @(n,V) gK   * n.^4      .* (V + EK);
IL   = @(V)        gL   * (V - EL);
ISyn = @(V)        gsyn * (V - Esyn);

%Intermediate Functions
k1    = @(V)  nA *(nAV + V) ./ (1-exp((-nAV-V)/nAk));
k2    = @(V)  nB * exp((-V-nBV)/nBk);
tau_n = @(V)     1 ./ (k1(V) + k2(V));
n_inf = @(V) k1(V) ./ (k1(V) + k2(V));

DF = @(t, x) [ (-INa(x(4), x(2), x(7), x(1)) - INaP(x(5), x(3), x(1)) - IK(x(6), x(1)) - IL(x(1)) - ISyn(x(1))) / C; ...
               (1 ./ (1 + exp(-(x(1) - VhNa )/khNa )) - x(2)) .* (cosh((x(1)-VhNat) /khtNa )) / thNamax ; ...
               (1 ./ (1 + exp(-(x(1) + VhNaP)/(-1*khNaP))) - x(3)) .* (cosh((x(1)+VhNaPt)/khtNaP)) / thNaPmax; ...
               (1 ./ (1 + exp(-(x(1) - VmNa )/kmNa )) - x(4)) .* (cosh((x(1)-VmNat) /kmtNa )) / tmNamax ; ...
               (1 ./ (1 + exp(-(x(1) - VmNaP)/kmNaP)) - x(5)) .* (cosh((x(1)-VmNaPt)/kmtNaP)) / tmNaPmax; ...
               (n_inf(x(1)) - x(6)) / tau_n(x(1)); ...
               ((1./(1 + exp(-(x(1)+Vh)/(-1*sh)))) - x(7)) .* ((cosh((x(1)+Vht)/sht))/thNa2)];
IC = [-57; 0.3; 0.3; 0.1; 0; 0; 1];

T1 = 25000;                                 % Inital equilibration time
T2 = 10000;                                 % Duration of the main run
options = odeset('RelTol', 1e-5);

[T, Y] = ode15s( DF, [0 T1], IC, options);  % Pre-equilibrate duration = T1
Vlist  = [Y(end,1)];                        % Initialize voltage array / list
Tlist  = [0];                               % Initialize time array / list
Cost   = 0;                                 % Initialize Cost variable
cnt    = 1;                                 % Initialize integration counter

                                            % Cost=0 returned if not enough
while Cost == 0 && cnt < 5                  % bursts detected; integrate at most 4 times
    cnt = cnt + 1;                          % Increment integraiton counter
    [T, Y] = ode15s( DF, [0 T2], Y(end, :), options);  % Integrate
    Vlist = [Vlist; Y(:,1)];                           % Append V(t) list
    Tlist = [Tlist; T + Tlist(end)];                   % Append time list
    Cost = ComputeCost_dspk(Tlist, Vlist, plotFlag);    % Compute Cost
end

if cnt > 4
    Cost = 100;
end


end
 
