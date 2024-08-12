function [Z, RE, mRE] = ROSSmain_v1_1(N,shswt,shinp,tvec,v0,delta,alpha,beta,lambda,U)
% ROSS Main Code
%   Z = ROSSmain_v1_1(N,shswt,shinp,tvec,v0,delta,alpha,beta,lambda,tau,U,s)
%   produces the cell array Z of size 1x(steps+1) of the free boundary evolution. 
%
%   [Z, RE, mRE] = ROSSmain_v1_1(N,shswt,shinp,tvec,v0,delta,alpha,beta,lambda,tau,U,s)
%   produces the cell arrays Z, the vector RE for the relative error of the
%   rate of change of area law at each time step and the scalar mRE: the
%   maxiumum relative error for the entire motion.
%
% INPUTS
%   N           = series truncation of the conformal map power series c_k zeta^k.
%   shswt       = shape switch = 0 (Laurent shape) , 1 (non-Laurent shape).
%   shinp       = shape input -- LA_shape and SC_shape functions.
%   tvec        = amplitude of noise.
%                   Number of unknowns = n = 2N+3.
%   v0          = basic rate of spread (ROS).
%   delta       = curvature parameter -- set to 0 in [1].
%   alpha       = ratio of radiative to convective ROS effects = [0,1].
%   beta        = pyrogenic wind parameter.
%   lambda      = ambient wind parameter.
%   U           = ambient wind vector.
%
% OUTPUTS
%   Z           = cell array of size 1x(steps+1) of free boundary data.
%                   Each entry (t-1) gives a list of complex z=x+iy values 
%                   for points with coordinates (x,y) on the free boundary 
%                   at time step t. (t=0: initial shape.)
%
% OPTIONAL OUTPUTS
%   RE       = vector of the relative error the relative error of the
%               rate of change of area law at each time step. 
%
%   mRE      = maximum relative error for the entire free boundary motion.
%
% REFERENCES
%
% [1]   Wildfire + Spotfire paper
%
% END OF DOCUMENTATION
%
%Code 
% Setup
n=2*N+3; t0 = tvec(1); % number of unknowns; initial time t0.
cdot0_est = 0.0*ones(1,n); % assume initial derivatives of c0 are zero.

if shswt==0 % Laurent shape
    c0_est = shape_LA_v1(shinp,N,n); % initial coefficients of cmap series.
    g = @(zeta) zeta.^(-1); % initial cmap (zeta^{-1} for Laurent shape).
else % non-Lauent (Schwarz-Christoffel) shape
    c0_est = 0.0*ones(1,n); c0_est(1)=1; % initial coefficients of cmap series.
    g = shape_SC_v1(shinp); % initial cmap.
end
 
% Eqns, decic and ode15i
f = @(t,c,cdot)ROSSeqns(t,c,cdot,g,N,n,v0,delta,alpha,beta,lambda,U); % defining ftn f based on eqns.
rtoly=1e-7; atoly=1e-7; % define relative and absolute tolerances.
options = odeset('RelTol', rtoly, 'AbsTol', atoly); % tolerances for decic.
options2 = odeset('RelTol', rtoly, 'AbsTol', atoly); % tolerances for ode15i.

[c0,cdot0] = decic(f,t0,c0_est,1,cdot0_est,[],options); % initial conditions for ode15i.
[t,coef] = ode15i(f,tvec,c0,cdot0,options2); % run ode15i.
coeff=coef.'; t=t.'; tsep = size(t,2); % flip coef and t; tsep = true final step (if error then tsep<size(tvec)).

% Constructing the free boundary
Pts=1000; theta=linspace(-pi,pi,Pts); zeta = exp(1i*theta); % unit zeta-disk.
gd = gradient(g(zeta))./gradient(zeta); % derivative of g: used in RCA law.
Area = 0.0*ones(1,tsep); RHS = 0.0*ones(1,tsep); % area vector, right hand side of RCA law vector.

for k=1:tsep
    z=coeff(1,k).*g(zeta); % this is the free boundary z.
    zfz = coeff(1,k).*zeta.*gd; %this is zeta.*df/dzeta.
    for kk=2:N+2
        z=z+(coeff(kk,k)+1i*coeff(N+1+kk,k))*zeta.^((kk-2));
        zfz=zfz+((kk-2)).*(coeff(kk,k)+1i*coeff(N+1+kk,k)).*zeta.^((kk-2));
    end

    area = polyarea(real(z),imag(z)); % area enclosed by curve z.
    rhs = RCArhs(z,theta,zfz,delta,alpha,beta,lambda,U); % RHS of RCA law.

    Z{k}=z; Area(k) = area; RHS(k) = rhs; % update cell array/vector.
end

dArea = gradient(Area)./gradient(t); % vector of rate of change of area in time.
RE = abs((dArea-RHS)./dArea); mRE = max(RE); % relative error, max rel error
end

%% Appendix A: ROSSeqns
function res = ROSSeqns(t,c,cdot,g,N,n,v0,delta,alpha,beta,lambda,U)
theta = linspace(-pi,pi,n+1); zeta = exp(1i*theta); % unit zeta-disk.

gd = gradient(g(zeta))./gradient(zeta); % first deriavtive of initial map g.
gdd = gradient(gd)./gradient(zeta); % second derivative of initial map g.

sum1=cdot(1).*g(zeta); %This is dot(z)
sum2=c(1).*zeta.*gd; %This is zeta*dz/dzeta
sum3=c(1).*zeta.*(zeta.*gdd+gd); %This is zeta*d/dzeta(zeta*dz/dzeta)

for m=2:N+2
sum1=sum1+(cdot(m)+1i*cdot(N+1+m)).*zeta.^((m-2));
sum2=sum2+((m-2)).*(c(m)+1i*c(N+1+m)).*zeta.^((m-2));
sum3=sum3+((m-2)).^2.*(c(m)+1i*c(N+1+m)).*zeta.^((m-2));
end

for k=1:n % take n points on zeta-disk, n equations for n unknowns.
radterm = max(0,alpha.*v0.*abs(sum2(k))+delta*(real(sum3(k).*conj(sum2(k)))./(abs(sum2(k)))^2));
convterm = max(0,(1-alpha).*v0.*abs(sum2(k))-beta-lambda.*real(sum2(k).*conj(U)));
res(k)=-real(sum1(k).*conj(sum2(k)))-radterm-convterm; res = res.';
end

if rem(n,2)~=0 % flip matrix if needed.
res=res.';
end
end

%% Appendix B: RHS of RCA Law
function rhs = RCArhs(z,theta,zfz,delta,alpha,beta,lambda,U)
L = perimeter(polyshape(real(z),imag(z)));
awterm = lambda.*real(zfz.*conj(U)); %ambient wind term
intgr = max(0,((1-alpha).*abs(zfz)-beta-awterm));
cterm = trapz(theta,intgr);

rhs = alpha.*L-2*pi*delta+cterm;
end
