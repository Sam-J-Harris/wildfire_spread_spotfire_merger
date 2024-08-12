%% Appendix A: Fire Timestepping function
function [bigz,tmax] = ROSAtstep_v1_4(bigz,bigc,tmax,merdata,mcnt,J,v0,delta,alpha,beta,lambda,U,tstep,rkswt,pcswt,resl,imswt)
% = calculates normal velocity step delta z, then timesteps to find z_{t+1}.
% Code:
k1 = firestep(bigz,bigc,J,v0,delta,alpha,beta,lambda,U,pcswt,imswt); % finding delta Z
[bigz, tmax] = fireRK(k1,bigz,bigc,tmax,merdata,mcnt,J,v0,delta,alpha,beta,lambda,U,tstep,rkswt,pcswt,resl,2,imswt); % computes RK0, RK2 or RK4 timestepping
end

%% Appendix A1: Single Timestep
function bigDel = firestep(bigz,bigc,J,v0,delta,alpha,beta,lambda,U,pcswt,imswt)
% = finds normal velocity step delta z.
% Code:
for j=1:J
    dz = bigz{j} - circshift(bigz{j},-1); v = 1i.*(dz)./abs(dz); 
    nVec{j} = -(v+circshift(v,1))./abs(v+circshift(v,1)); % normal vector for each fire
    curv{j} = LineCurv(bigz{j}); % curvature
end

Pwind = nVec; % if no pyrowind effect
if beta ~=0
    Pwind = fireaaa(bigz,bigc,J,nVec,pcswt,imswt); % find pyrowind dphi/dn (AAA-LS)
end

for j=1:J % find delta Z for each fire line
    bigDel{j}=(alpha.*v0 - delta.*curv{j} + max(0, (1-alpha).*v0 + beta.*Pwind{j}+lambda.*cdot(U,nVec{j}))).*nVec{j};
end  
end

%% Appendix A2: AAA-least squares function
function Pwind=fireaaa(bigz,bigc,J,nVec,pcswt,imswt)
% = calculates pyrogenic wind effect dphi/dn by solving for phi in the fire line exterior.
% Uses multiply connected AAA-LS method - see Costa & Trefethen 2023.
% Code:
inpolygonc = @(z,w) inpolygon(real(z),imag(z),real(w),imag(w)); % ftn to determine if point is inside polygon
rtoly = pcswt*1e-13; ztoly = pcswt*1e-2; % residue/pole-zero tolerances (set pcswt=0 for no pole control)

Z=[]; for j =1:J, Z=[Z; bigz{j}]; end % setting up Dirichlet boundary condition
H=log(abs(h(Z,bigc,J))); %see ftn "h" below. Set phi=H on bdry then subtract H at the end.

Allpols=[]; % list of all poles from AAA
A1=[]; % real part of log terms
A2=[]; % real part of Runge (series) and Neuman (singularity) terms
A3=[]; % 0 vector (= imag part of log terms)
A4=[]; % imag part of Runge (series) and Neuman (singularity) terms

[~,polk,resk] = aaa(H,Z,'cleanup',1,'tol',rtoly); % global AAA poles
for j=1:J %for each spotfire j
    pol_in = polk(inpolygonc(polk,bigz{j})); res_in = resk(inpolygonc(polk,bigz{j})); % determines which poles are inside polygon j
    
    Pol = pol_in(abs(res_in)>rtoly).'; % poles of residue < rtoly eliminated (pole control)
    D = min(abs((bigz{j})-Pol),[],1); D(D<ztoly*max(abs(bigz{j}))) = 0; D(D~=0) =1; % dist. b/ween poles and bdry
    pol=Pol.*D; pol(abs(pol)==0)=[]; % poles close to the boundary < ztoly eliminated (pole control)
    Pols{j}=pol; Allpols=[Allpols Pols{j}]; %remaining poles put into pole list

    dvals{j}=min(abs(bigz{j}-Pols{j}),[],1); % closest pole to each point on polygon j
    [Hes,P] = VAorthog(1./(Z-bigc{j}),20); % Hes = Hessenberg matrices used in VAeval later; P = bases of Runge (series) terms
    Q=dvals{j}./(Z-Pols{j}); % Q = bases of Neuman (singularity) terms
    bigHes{j}=Hes; bigP{j}=P; bigQ{j}=Q; % big lists of Hes, P and Q
    
    if j~=J % updating log basis (care taken if j==J)
        A1 = [A1 log(abs((Z-bigc{j})./(Z-bigc{j+1})))];
    else
        A1 = [A1 log(abs((Z-bigc{j})./(Z-bigc{1})))];
    end
    A3 = [A3 0*Z];
    A2 = [A2 real(bigP{j}) real(bigQ{j})]; A4 = [A4 -imag(bigP{j}) -imag(bigQ{j})]; % update Runge and Newman bases
end

A = [A1 A2 A3 A4]; % matrix A of basis vectors created
c = reshape(A\H,[],2)*[1; 1i]; % c = vector of coefficients. On the boundary, phi = real(f(z)) = A*c = H => c = A\H.

phi = @(z) -log(abs(h(z,bigc,J)))+real(f(z,c,bigc,J,bigHes,dvals,Pols)); % phi with H subtracted off - see  ftn "f" below (for psi, see v1_6).
ROSApimage_v1_1(Allpols,bigz,J,phi,imswt); % contour plot of phi at timestep (if nec)

eps = 1e-8; % compute dphi/dn = phi(z+ndz) - phi(z) / eps
for j=1:J
    z = bigz{j}; znew = z+eps*nVec{j};
    Pwind{j} = (phi(znew)-phi(z))./abs(eps); %pyrogenic wind term
end

%%% Additional Functions
function bigprod=h(Z,bigc,J) % phi->-H=-log|h| as z->infty
    bigprod = 1;
    for k=1:J
        bigprod=bigprod.*(Z-bigc{k});
    end
end

function F = f(z,c,bigc,J,bigHes,dvals,Pols) % phi = real(F(z))
    B1=[]; B2=[];
    for n=1:J % i.e. for fires 1,2,...,J
        if n~=J %correct log term used
            B1 = [B1 log((z(:)-bigc{n})./(z(:)-bigc{n+1}))];
        else
            B1 = [B1 log((z(:)-bigc{n})./(z(:)-bigc{1}))];
        end
        B2 = [B2 VAeval(1./(z(:)-bigc{n}),bigHes{n}) dvals{n}./(z(:)-Pols{n})]; % evaluated Runge (series) and Neuman (singularity) bases exterior to the spotfires
    end
    B=[B1 B2]; F = reshape(B*c,size(z)); % B=matrix of log, Runge and Neuman bases, then F = B*c
end
end

%% Appendix A3: Runge Kutta Function - NEED TO FIX
function [bigz, tmax] = fireRK(k1,bigz,bigc,tmax,merdata,mcnt,J,v0,delta,alpha,beta,lambda,U,tstep,rkswt,pcswt,resl,inswt,imswt)
% = timesteps using Runge Kutta approach, either Euler's method (RK1), second order RK (RK2) or fourth order RK (RK4).
% Code:
tstepa = tstep; tstepm = 0.00075; % actual tstep value (changes if emergency RK1 used); modified tstep for RK1
if rkswt==0 % RK1
    for j = 1:J, bigz{j} = bigz{j}+tstep*k1{j}; end % computes RK1 timestep
else % either RK2 or RK4
    olap=1; sint = 1; % overlapping and self intersect counters
    for j = 1:J, bigz1t{j} = bigz{j} + (tstep/2).*k1{j}; end % update fire line in prep for k2
    for j = 1:J, sint = sintchk(bigz1t{j},sint); end % see if any self intersects have been deleted
    for j = 1:J-1, olap = olapchk(bigz1t{j}, bigz1t{j+1},olap); end % see if any fire lines overlap others
    if olap==0 % if fire lines overlap, do RK1 with modified tstep
        tstepa = tstepm; for j = 1:J, bigz{j} = bigz{j}+tstepa*k1{j}; end % computes RK1 timestep 
    elseif sint==0 % if fire line self intersected, do RK1 with modified tstep
        tstepa = tstepm; for j = 1:J, bigz{j} = bigz{j}+tstepa*k1{j}; end % computes RK1 timestep 
    else
        bigz1 = ROSAsmooth_v1_4(bigz1t,mcnt,J,resl,inswt,imswt); % fire line smoothing
        k2 = firestep(bigz1,bigc,J,v0,delta,alpha,beta,lambda,U,pcswt,imswt); % compute k2
        for j = 1:J, zrk2{j} = bigz{j}+tstep*k2{j}; end % computes RK2 timestep
        if rkswt ==2 %RK2
            bigz = zrk2;
        else    % RK4
            for j = 1:J, bigz2t{j} = bigz{j} + (tstep/2).*k2{j}; end % update fire line in prep for k3
            
            for j = 1:J, sint = sintchk(bigz2t{j},sint); end % see if any self intersects have been deleted
            for j = 1:J-1, olap = olapchk(bigz2t{j}, bigz2t{j+1},olap); end % see if any fire lines overlap others
            
            if olap==0 || sint==0, bigz = zrk2; % if a merge has happened, just do RK2
            else
                bigz2 = ROSAsmooth_v1_4(bigz2t,mcnt,J,resl,inswt,imswt); % fire line smoothing
                k3 = firestep(bigz2,bigc,J,v0,delta,alpha,beta,lambda,U,pcswt,imswt); % compute k3
                for j = 1:J, bigz3t{j} = bigz{j} + tstep.*k3{j}; end % update fire line in prep for k4
                
                for j = 1:J, sint = sintchk(bigz3t{j},sint); end % see if any self intersects have been deleted
                for j = 1:J-1, olap = olapchk(bigz3t{j}, bigz3t{j+1},olap); end % see if any fire lines overlap others
                
                if olap==0 || sint==0, bigz = zrk2; % if a merge has happened, just do RK2
                else
                    bigz3 = ROSAsmooth_v1_4(bigz3t,mcnt,J,resl,inswt,imswt); % fire line smoothing
                    k4 = firestep(bigz3,bigc,J,v0,delta,alpha,beta,lambda,U,pcswt,imswt); % compute k4
                    for j = 1:J
                        bigz{j} = bigz{j}+(tstep/6)*(k1{j}+2.*k2{j}+2.*k3{j}+k4{j}); end % computes RK4 timestep
                end
            end
        end
    end
end
tmax = tmax + tstepa; % updates time
end

%% Appendix A4: Complex Dot Product
function dotprod = cdot(z1,z2) % dot product of two complex numbers
x1 = real(z1); y1 = imag(z1); x2 = real(z2); y2 = imag(z2);
dotprod = x1*x2 + y1*y2;
end

%% Appendix A5: Check if fires self-intersect
function sint = sintchk(z, sinti)
% = determine if a fire line intersects itself. Uses selfintersect function from [REF online].
zt = z(1:end-1); xt = real(zt); yt=imag(zt);
[x0,y0,segments]=selfintersect(xt,yt); znew = z; % uses selfintersect function - see REF.

if size(segments,2)~=0 && size(segments,1)~=0 % ie there are overlapping segments
    sint = 0;
else
    sint = sinti;
end
end

%% Appendix A6: Check if fires overlap
function olap = olapchk(z1, z2, olapi)
x1 = real(z1); y1 = imag(z1); x1 = x1(~isnan(x1)); y1 = y1(~isnan(y1));
x2 = real(z2); y2 = imag(z2); x2 = x2(~isnan(x2)); y2 = y2(~isnan(y2)); % remove any NaN points from data
poly1 = polyshape(x1,y1); poly2 = polyshape(x2,y2); % convert to polyshape objects

if overlaps(poly1,poly2)
    olap=0; 
else
    olap = olapi;
end % determine if the fire lines overlap
end