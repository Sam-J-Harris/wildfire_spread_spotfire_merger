function [Z, RE, mRE] = ROSSOmain_v7_1(N,problem,input,tvec,v0,delta,alpha,beta,lambda,tau,U,s)
% Setup
n=2*N+3; c0_est = 0.0*ones(1,n); c0_est(1)=1; cdot0_est = 0.0*ones(1,n); t0 = tvec(1);
g=@(zeta) zeta.^(-1); rtoly=1e-7; atoly=1e-7; %tolerance set up

if problem==0 %Laurent shape
    c0_est = LA_shape_map_v1(input,N,n);
elseif problem==1 %Schwarz-Christoffel shape
    g = SC_shape_map_v7(input);
else    %Derrida-Hakim shape
    g = DH_shape_map_v1(input);
end
 
% Eqns, decic and ode15i
f = @(t,c,cdot)ROSSOeqns(t,c,cdot,g,N,n,v0,delta,alpha,beta,lambda,tau,U,s); %Defining ftn f based on eqns
options = odeset('RelTol', rtoly, 'AbsTol', atoly); %tolerances for decic
options2 = odeset('RelTol', rtoly, 'AbsTol', atoly); %tolerances for decic

[c0,cdot0] = decic(f,t0,c0_est,1,cdot0_est,[],options); %initial conditions for ode15i
[t,coef] = ode15i(f,tvec,c0,cdot0,options2); coeff=coef.'; t=t.'; tsep = size(t,2);

% Constructing the free boundary
Pts=1000; theta=linspace(-pi,pi,Pts); zeta = exp(1i*theta);
gd = gradient(g(zeta))./gradient(zeta); %g' : used in RCA law
RHS = 0.0*ones(1,tsep); Area = 0.0*ones(1,tsep);

for k=1:tsep
    z=coeff(1,k).*g(zeta);
    sum2 = coeff(1,k).*zeta.*gd; %this is zeta.*df/dzeta
    for kk=2:N+2
        z=z+(coeff(kk,k)+1i*coeff(N+1+kk,k))*zeta.^((kk-2));
        sum2=sum2+((kk-2)).*(coeff(kk,k)+1i*coeff(N+1+kk,k)).*zeta.^((kk-2));
    end
    rhs = RCArhs(z,theta,sum2,delta,alpha,beta,lambda,U); %THIS IS WRONG, FIX
    area = polyarea(real(z),imag(z));
    Z{k}=z; RHS(k) = rhs; Area(k) = area;
end

dArea = gradient(Area)./gradient(t);
RE = abs((dArea-RHS)./dArea); mRE = max(RE); % relative error, max rel error
end

%% Appendix A: ROSSOeqns
function res = ROSSOeqns(t,c,cdot,g,N,n,v0,delta,alpha,beta,lambda,tau,U,s)
theta = linspace(-pi,pi,n+1); zeta = exp(1i*theta);

gd = gradient(g(zeta))./gradient(zeta); %gd2 = evaldiff(g,zeta); %alt gd
gdd = gradient(gd)./gradient(zeta);

sum1=cdot(1).*g(zeta); %This is dot(z)
sum2=c(1).*zeta.*gd; %This is zeta*dz/dzeta
sum3=c(1).*zeta.*(zeta.*gdd+gd); %This is zeta*d/dzeta(zeta*dz/dzeta)

for m=2:N+2
sum1=sum1+(cdot(m)+1i*cdot(N+1+m)).*zeta.^((m-2));
sum2=sum2+((m-2)).*(c(m)+1i*c(N+1+m)).*zeta.^((m-2));
sum3=sum3+((m-2)).^2.*(c(m)+1i*c(N+1+m)).*zeta.^((m-2));
end

for k=1:n
radterm = max(0,alpha.*v0.*abs(sum2(k))+delta*(real(sum3(k).*conj(sum2(k)))./(abs(sum2(k)))^2));
convterm = max(0,(1-alpha).*v0.*abs(sum2(k))-beta-lambda.*real(sum2(k).*conj(U)));
res(k)=-real(sum1(k).*conj(sum2(k)))-radterm-convterm;  res = res.';
end

if rem(n,2)~=0
res=res.';
end
end

%% Appendix B: RHS of RCA Law
function rhs = RCArhs(z,theta,sum2,delta,alpha,beta,lambda,U)
L = perimeter(polyshape(real(z),imag(z)));
awterm = lambda.*real(sum2.*conj(U)); %ambient wind term
intgr = max(0,(alpha-beta-awterm));
cterm = trapz(theta,intgr);

rhs = L-2*pi*delta+cterm;
end
