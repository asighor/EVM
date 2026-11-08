function varargout = EVM_models(model, quantity, h, p, varargin)
%EVM_MODELS  Reference implementation of the soil hydraulic functions compared in the
% manuscript: the van Genuchten (VG), Peters-Durner-Iden (PDI) and extended van
% Genuchten-Mualem (EVM) models. The same code is used in Scripts S1-S6.
%
%   theta = EVM_models(model, 'theta', h, p)          water content (cm3 cm-3)
%   C     = EVM_models(model, 'capacity', h, p)       C = -dtheta/dh (cm-1)
%   [Sc, Snc] = EVM_models(model, 'saturation', h, p) capillary and non-capillary saturation
%   [K, Kc, Knc, Kv] = EVM_models(model, 'K', h, p, 'Ksc', Ksc)   conductivity (cm d-1)
%
%   model : 'VG', 'PDI' or 'EVM'
%   h     : suction magnitude (cm), h >= 0 (pF = log10 h)
%   p     : [alpha (cm-1), n (-), theta_r, theta_s]
%   options (name-value): 'Ksc' saturated capillary conductivity (cm d-1, required for 'K'),
%     'h0' oven-dry suction (1e7 cm), 'c' EVM closure coefficient (0.1; beta = 1 + c/n),
%     'KncForm' EVM non-capillary conductivity, 'powerlaw' = Eq. (20) (default) or
%     'integral' = Eq. (B2) (integral with Eq. 19), 'T' temperature for Kv (20 C).
%
% Equations (ha = 1/alpha, m = 1 - 1/n, G = [1 + (alpha h)^n]^-m, G0 = G(h0)):
%   VG  : theta = thr + (ths - thr) G                                          (van Genuchten, 1980)
%         K = Ksc G^0.5 [1 - (1 - G^(1/m))^m]^2                                 (Mualem, 1976a)
%   PDI : theta = (ths - thr) Sc + thr Snc,  Sc = (G - G0)/(1 - G0)              (Iden and Durner, 2014)
%         Snc = 1 - [max(u,0) + b ln(1 + exp(-|u|/b))]/log10(alpha h0),  u = log10(alpha h)
%         b = 0.1 + 0.2/n^2 [1 - exp(-(thr/(ths - thr))^2)]
%         K = Kc + Knc + Kv,  Knc = Ks,nc (alpha h0)^(-1.5 (1 - Snc))           (Eq. 16, Peters, 2013)
%   EVM : theta = (ths - thr) Sc + thr Snc,  Snc = 1 - ln(1 + (alpha h)^beta)/ln(1 + (alpha h0)^beta)
%         beta = 1 + c/n                                                         (Eqs. 6-7)
%         K = Kc + Knc + Kv,  Knc = Ks,nc (alpha h0)^(-1.5 (1 - Snc))           (Eq. 20, l = 0, eta = 2, g = 0.75)
%         ('integral': Knc = Ks,nc [((alpha h0)^(g Snc) - 1)/((alpha h0)^g - 1)]^2, Eq. B2)
%   PDI and EVM: Kc = Ksc Sc^0.5 {1 - [(1 - G^(1/m))/(1 - G0^(1/m))]^m}^2    (Eq. 14, Ks,c once)
%         Ks,nc = 86.2 theta(1e5 cm) alpha^1.5 (cm d-1)                        (Eq. 17, Peters et al., 2021)
%         Kv = (rho_sv/rho_w) D (Mw g/(R T)) exp(-h Mw g/(R T)), D = Da (ths - theta)^(10/3)/ths^2  (Appendix A)
% The capillary Mualem term is evaluated in logarithms, so that it does not underflow as n -> 1.
%
% Example:
%   h = logspace(-2, 7, 200)'; p = [0.02 1.5 0.1 0.45];
%   th = EVM_models('EVM', 'theta', h, p);  C = EVM_models('PDI', 'capacity', h, p);
%   [K, Kc, Knc, Kv] = EVM_models('EVM', 'K', h, p, 'Ksc', 50);

o = struct('Ksc', NaN, 'h0', 1e7, 'c', 0.1, 'KncForm', 'powerlaw', 'T', 20);
for k = 1:2:numel(varargin), o.(varargin{k}) = varargin{k + 1}; end
model = upper(char(model)); h = double(h);
switch lower(char(quantity))
    case 'theta'
        varargout{1} = shp_theta(h, p, model, o);
    case 'capacity'
        varargout{1} = shp_capacity(h, p, model, o);
    case 'saturation'
        [varargout{1}, varargout{2}] = shp_sat(h, p, model, o);
    case 'k'
        if ~isfinite(o.Ksc), error('EVM_models: option ''Ksc'' is required for K.'); end
        [lK, Kc, Knc, Kv] = shp_logK(h, p, model, o.Ksc, o);
        varargout = {10.^lK, Kc, Knc, Kv};
    otherwise
        error('EVM_models: unknown quantity ''%s''.', quantity);
end
end

function [Sc, Snc] = shp_sat(h, p, model, o)
al = p(1); n = p(2); m = 1 - 1/n;
G = (1 + (al*h).^n).^(-m);
if strcmp(model, 'VG'), Sc = G; Snc = zeros(size(h)); return, end
G0 = (1 + (al*o.h0)^n)^(-m);
Sc = (G - G0)/(1 - G0);
Snc = noncap_sat(h, p, model, o);
end

function Snc = noncap_sat(h, p, model, o)
al = p(1); n = p(2); thr = p(3); ths = p(4);
if strcmp(model, 'PDI')
    if ths == thr, r2 = Inf; else, r2 = (thr/(ths - thr))^2; end
    b = 0.1 + 0.2/n^2*(1 - exp(-r2));
    u = log10(al*h);
    Snc = 1 - (max(u, 0) + b*log1p(exp(-abs(u)/b)))/log10(al*o.h0);
else
    beta = 1 + o.c/n;
    Snc = 1 - log1p((al*h).^beta)/log1p((al*o.h0)^beta);
end
end

function th = shp_theta(h, p, model, o)
thr = p(3); ths = p(4);
[Sc, Snc] = shp_sat(h, p, model, o);
if strcmp(model, 'VG'), th = thr + (ths - thr)*Sc; return, end
th = (ths - thr)*Sc + thr*Snc;
end

function C = shp_capacity(h, p, model, o)
al = p(1); n = p(2); thr = p(3); ths = p(4); m = 1 - 1/n; x = al*h;
dG = m*n*al*x.^(n - 1)./(1 + x.^n).^(m + 1);                 % -dGamma/dh
if strcmp(model, 'VG'), C = (ths - thr)*dG; return, end
G0 = (1 + (al*o.h0)^n)^(-m);
if strcmp(model, 'PDI')
    if ths == thr, r2 = Inf; else, r2 = (thr/(ths - thr))^2; end
    b = 0.1 + 0.2/n^2*(1 - exp(-r2));
    u = log10(x);
    dSnc = 1./(1 + exp(-u/b))./(h*log(10)*log10(al*o.h0));   % -dSnc/dh
else
    beta = 1 + o.c/n;
    dSnc = beta*al*x.^(beta - 1)./((1 + x.^beta)*log1p((al*o.h0)^beta));
end
C = (ths - thr)*dG/(1 - G0) + thr*dSnc;
end

function [lK, Kc, Knc, Kv] = shp_logK(h, p, model, Ksc, o)
al = p(1); n = p(2); thr = p(3); ths = p(4); m = 1 - 1/n;
lG = -m*log1p((al*h).^n);
if strcmp(model, 'VG')
    lG0 = -Inf; lSc = lG;
else
    lG0 = -m*log1p((al*o.h0)^n);
    lSc = lG + log1p(-exp(lG0 - lG)) - log1p(-exp(lG0));
end
le = lG/m; ld = lG0/m;
u  = (exp(le) - exp(ld))/(1 - exp(ld));
lB = log(-expm1(m*log1p(-u)));
t  = u < 1e-10;
lB(t) = log(m) + le(t) + log1p(-exp(ld - le(t))) - log1p(-exp(ld));
lKc = log(Ksc) + 0.5*lSc + 2*lB;
Kc = exp(lKc);
if strcmp(model, 'VG')
    lK = lKc/log(10); Knc = zeros(size(h)); Kv = zeros(size(h));
    return
end
[Sc, Snc] = shp_sat(h, p, model, o); Sc = max(Sc, 0); Snc = max(Snc, 0);
Ksnc = 86.2*shp_theta(1e5, p, model, o)*al^1.5;
Rv = al*o.h0;
if strcmp(model, 'EVM') && strcmpi(o.KncForm, 'integral')
    Knc = Ksnc*(expm1(0.75*Snc*log(Rv))/expm1(0.75*log(Rv))).^2;     % Eq. (B2), integral with Eq. (19)
else
    Knc = Ksnc*Rv.^(-1.5*(1 - Snc));                                  % Eq. (16) / Eq. (20)
end
th = (ths - thr)*Sc + thr*Snc;
Kv = vapor_K(h, th, ths, o.T);
lK = log10(Kc + Knc + Kv);
end

function Kv = vapor_K(h, th, ths, Tc)
T = 273.15 + Tc; Mw = 0.018015; g = 9.81; Rg = 8.314;
rho_sv = 1e-3*exp(31.3716 - 6014.79/T - 7.92495e-3*T)/T;     % kg m-3
rho_w  = 998.2;                                              % kg m-3
Da  = 2.12e-5*(T/273.15)^2*1e4*86400;                        % cm2 d-1
cH  = Mw*g/(Rg*T)/100;                                       % cm-1
tau = max(ths - th, 0).^(10/3)/ths^2;                        % Eq. (A4)
Kv  = rho_sv/rho_w*tau*Da*cH.*exp(-cH*h);                    % Eqs. (A1)-(A3)
end
