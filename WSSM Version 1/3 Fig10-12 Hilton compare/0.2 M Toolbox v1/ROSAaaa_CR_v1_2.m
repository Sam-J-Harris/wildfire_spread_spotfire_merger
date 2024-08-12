%% Appendix A: Fire Timestepping function
function bigz = ROSAaaa_CR_v1_2(bigz,bigc,merdata,J,v0,delta,alpha,beta,lambda,tau,U,s,tstep,rkswit,pcswit,resl,intswit,imswit)
    k1 = firestep(bigz,bigc,J,v0,delta,alpha,beta,lambda,tau,U,s,pcswit,imswit); %finding delta Z
    for j = 1:J
        zrk0{j} = bigz{j}+tstep*k1{j};
    end
    if rkswit==0 % RK0
        bigz = zrk0;
    else
        for j = 1:J
             z1temp = bigz{j} + (tstep/2).*k1{j};
             z1sm = smoothdata(z1temp,'gaussian',5); z1sm(end) = z1sm(1);
             bigz1{j} = z1sm;
        end
        [~, ~, merdataz1, ~] = ROSAmerger_CR_v1_2(bigz1,bigc,merdata,J,resl,intswit,imswit);
        if max(size(merdataz1))~=0
            bigz = zrk0;
        else
            k2 = firestep(bigz1,bigc,J,v0,delta,alpha,beta,lambda,tau,U,s,pcswit,imswit);
            for j = 1:J
                zrk2{j} = bigz{j}+tstep*k2{j};
            end
            if rkswit ==2 %RK2
                bigz = zrk2;
            else    % RK4
                for j = 1:J
                    z2temp = bigz{j} + (tstep/2).*k2{j};
                    z2sm = smoothdata(z2temp,'gaussian',5); z2sm(end) = z2sm(1);
                    bigz2{j} = z2sm;
                end
                [~, ~, merdataz2, ~] = ROSAmerger_CR_v1_2(bigz2,bigc,merdata,J,resl,intswit,imswit);
                if max(size(merdataz2))~=0
                    bigz = zrk0;
                else
                    k3 = firestep(bigz2,bigc,J,v0,delta,alpha,beta,lambda,tau,U,s,pcswit,imswit);
                    for j = 1:J
                        z3temp = bigz{j} + tstep.*k3{j};
                        z3sm = smoothdata(z3temp,'gaussian',5); z3sm(end) = z3sm(1);
                        bigz3{j} = z3sm;
                    end
                    [~, ~, merdataz3, ~] = ROSAmerger_CR_v1_2(bigz3,bigc,merdata,J,resl,intswit,imswit);
                    if max(size(merdataz3))~=0
                        bigz = zrk0;
                    else
                        k4 = firestep(bigz3,bigc,J,v0,delta,alpha,beta,lambda,tau,U,s,pcswit,imswit);
                        for j = 1:J
                            bigz{j} = bigz{j}+(tstep/6)*(k1{j}+2.*k2{j}+2.*k3{j}+k4{j});
                        end
                    end
                end
            end
        end
    end
end

%% Appendix A1: Single Timestep
function bigDel = firestep(bigz,bigc,J,v0,delta,alpha,beta,lambda,tau,U,s,pcswit,imswit)
    for j=1:J
        dz = bigz{j} - circshift(bigz{j},-1); v = 1i.*(dz)./abs(dz); 
        nVec{j} = -(v+circshift(v,1))./abs(v+circshift(v,1)); %normal vector for each fire
        curv{j} = LineCurv(bigz{j}); %curvature
    end
    
    %Pyrogenic Wind
    if beta ~=0
        Pwind = timestep_AAA(bigz,bigc,J,nVec,pcswit,imswit); %find pyrowind dphi/dn (AAA-LS)
    else
        Pwind = nVec;
    end

    %Find deltaz for each fire
    for j=1:J
        bigDel{j}=(alpha.*v0 - delta.*curv{j} + max(0, (1-alpha).*v0 + beta.*Pwind{j}+lambda.*cdot(U,nVec{j}))).*nVec{j};
    end  
end

%% Appendix A2: AAA-least squares function
function Pwind=timestep_AAA(bigz,bigc,J,nVec,pcswit,imswit)
    % Pole tolerance setup (for pole control)
    if pcswit ==0
        rtoly=0; ztoly=0;
    else
        rtoly = 1e-13; ztoly = 1e-2; % residue/pole-zero tolerances (set =0 for no pole control)
    end
    
    % Dirichlet boundary condition
    Z=[];
    for j =1:J
        Z=[Z; bigz{j}];
    end
    H=log(abs(h(Z,bigc))); %see ftn "h" below. Set phi=H on bdry then subtract H at the end.

    % AAA poles 
    [~,polk,resk] = aaa(H,Z,'cleanup',1,'tol',rtoly); %global AAA
    inpolygonc = @(z,w) inpolygon(real(z),imag(z),real(w),imag(w)); % ftn to determine if pole is inside polygon
    Allpols=[]; 
    % vectors of bases for each spotfire
    A1=[]; % real part of log terms
    A2=[]; % real part of Runge (series) and Neuman (singularity) terms
    A3=[]; % 0 vector (= imag part of log terms)
    A4=[]; % imag part of Runge (series) and Neuman (singularity) terms
    
    % Pole control and least squares
    for j=1:J %for each spotfire j
        % Remove Froissart doublets ie pole control %%%%%%%%%%%
        pol_in = polk(inpolygonc(polk,bigz{j})); res_in = resk(inpolygonc(polk,bigz{j})); %determines which poles are inside polygon j
        
        Pol = pol_in(abs(res_in)>rtoly).'; % poles of small residue < rtoly eliminated
        D = min(abs((bigz{j})-Pol),[],1); D(D<ztoly*max(abs(bigz{j}))) = 0; D(D~=0) =1; % dist. b/ween poles and bdry
        pol=Pol.*D; pol(abs(pol)==0)=[]; % poles close to the boundary < ztoly eliminated
        Pols{j}=pol; Allpols=[Allpols Pols{j}]; %remaining poles put into pole list
    
        % Least Squares %%%%%%%%%%
        dvals{j}=min(abs(bigz{j}-Pols{j}),[],1); % closest pole to each point on polygon j
        [Hes,P] = VAorthog(1./(Z-bigc{j}),20); % Hes = Hessenberg matrices used in VAeval later; P = bases of Runge (series) terms
        Q=dvals{j}./(Z-Pols{j}); % Q = bases of Neuman (singularity) terms
        bigHes{j}=Hes; bigP{j}=P; bigQ{j}=Q; % big lists of Hes, P and Q
        
        % Updating A1,A2,A3 and A4 with correct bases
        if j~=J
            A1 = [A1 log(abs((Z-bigc{j})./(Z-bigc{j+1})))];
        else
            A1 = [A1 log(abs((Z-bigc{j})./(Z-bigc{1})))];
        end
        A3 = [A3 0*Z];
        A2 = [A2 real(bigP{j}) real(bigQ{j})]; A4 = [A4 -imag(bigP{j}) -imag(bigQ{j})];    
    end

    % Evaluate Solution %%%%%%%%%
    A = [A1 A2 A3 A4]; % matrix A created
    c = reshape(A\H,[],2)*[1; 1i]; % c = vector of coefficients. On the boundary, phi = real(f(z)) = A*c = H => c = A\H.
    phi = @(z) -log(abs(h(z,bigc)))+real(f(z,c,bigc,bigHes,dvals,Pols)); %pot phi with H subtracted off - see  ftn "f" below (for psi, see v1_6).

    eps = 1e-8; 
    for j=1:J
        z = bigz{j}; znew = z+eps*nVec{j};
        Pwind{j} = (phi(znew)-phi(z))./abs(eps); %pyrogenic wind term
    end

    
    %%%%%%% Plots %%%%%%%%%%%
    if imswit==1
        axl = 7; LW = 'linewidth'; MS = 'markersize'; ms = 6; PO = 'position'; FS = 'fontsize';
        figure(1)
        axes(PO,[.5 .3 .5 .5]);
        plot(Allpols,'.r',MS,ms), hold on
        x = linspace(-axl,axl,500); [xx,yy] = meshgrid(x,x); zz = xx+1i*yy; 
        for j=1:J
            if j==1
                cond = inpolygonc(zz,bigz{j});
            else
                cond = cond|inpolygonc(zz,bigz{j});
            end
        end
        uu = phi(zz); uu(cond) = NaN;
        for j=1:J % plotting boundaries
            zb = bigz{j}; xb = real(zb); yb = imag(zb); plot(xb,yb,'k',LW,.9), 
        end
        contour(x,x,uu,30,LW,1), hold off, 
        axis([-axl axl -axl axl]), axis square, title('Laplace solution');
    end
    
    %%%%%%% Additional Functions %%%%%%%%%%%
    function bigprod=h(Z,bigc) % phi->-H=-log|h| as z->infty
    bigprod = 1;
    for k=1:size(bigc,2)
        bigprod=bigprod.*(Z-bigc{k});
    end
    end
    
    function F = f(z,c,bigc,bigHes,dvals,Pols) % phi = real(F(z))
    B1=[]; B2=[];
    for n=1:size(bigc,2) % i.e. for fires 1,2,...,J
        if n~=size(bigc,2) %correct log term used
            B1 = [B1 log((z(:)-bigc{n})./(z(:)-bigc{n+1}))];
        else
            B1 = [B1 log((z(:)-bigc{n})./(z(:)-bigc{1}))];
        end
        B2 = [B2 VAeval(1./(z(:)-bigc{n}),bigHes{n}) dvals{n}./(z(:)-Pols{n})]; % evaluated Runge (series) and Neuman (singularity) bases exterior to the spotfires
    end
    B=[B1 B2]; F = reshape(B*c,size(z)); % B=matrix of log, Runge and Neuman bases, then F = B*c
    end
end

%% Appendix A3: Complex Dot Product
function dotprod = cdot(z1,z2)
    x1 = real(z1); y1 = imag(z1); x2 = real(z2); y2 = imag(z2);
    dotprod = x1*x2 + y1*y2;
end
