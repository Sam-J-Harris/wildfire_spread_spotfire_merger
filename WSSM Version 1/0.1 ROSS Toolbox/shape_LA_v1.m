function c0_est = shape_LA_v1(input,N,n)
% Outputs the c coefficients of the Laurent conformal map 
%   z = g(zeta) = c_{-1}*zeta.^{-1} + c_k*zeta^k.
%
% SHAPE DIRECTORY
%   11: Lemniscate
%   12: thin ellipse
%
%   The following correspond to figure numbers from [1], e.g. 21 = figure 2a.
%   21: circle
%   22: slanted ellipse
%   23: triangle 1
%   24: triangle 2
%   25: corrugated circle
%   26: irregular object
%   31: diamond
%   41: regular ellipse
%   51: irregular pentagon
%   71: dumbell
%   91: three-pronged object (fidget spinner)
%
% REFERENCES
%
%   [1]     Rycroft, C. H., & Bazant, M. Z. (2016). "Asymmetric collapse by 
%           dissolution or melting in a uniform flow". Proc. R. Soc. A, 
%           472(2185), 20150531.
%
% END OF DOCUMENTATION
%
%Code
c0_est=0.0*ones(1,n); A0 = 1; c0_est(1)=real(A0); %circle (if input==21)

if input==11 %lemniscate
    ee = 1.01; M=3;
    k=1; cff = 1/M; c0_est(1)=ee;
    while (k*M-1)<N
        q = (cff*ee^(-k*M+1))/factorial(k); pwr = k*M-1;
        c0_est(pwr+2) = real(q); c0_est(pwr+2+N+1)=imag(q);
        cff = cff*(1/M-k); k=k+1;
    end

elseif input==12 %thin ellipse
    C1 = 0.9;
    c0_est(3)=real(C1); c0_est(N+1+3)=imag(C1);


% --- Shapes from Rycroft and Bazant 2016 ---
elseif input==22 %slanted ellipse
    C1 = 0.3+0.2*1i;
    c0_est(3)=real(C1); c0_est(N+1+3)=imag(C1);

elseif input==23 %triangle 1
    C2 = -0.35;
    c0_est(4)=real(C2); c0_est(N+1+4)=imag(C2);

elseif input==24 %triangle 2
    C2 = 0.35*1i;
    c0_est(4)=real(C2); c0_est(N+1+4)=imag(C2);

elseif input==25 %corrugated circle
    C15 = 0.05*1i;
    c0_est(17)=real(C15); c0_est(N+1+17)=imag(C15);

elseif input ==26 %irregular object
    C1 = -0.28 + 0.2*1i; 
    c0_est(3)=real(C1); c0_est(N+1+3)=imag(C1); 
    C6 = 0.1;
    c0_est(8)=real(C6); c0_est(N+1+8)=imag(C6);

elseif input == 31 %diamond
    C3 = 0.25;
    c0_est(5)=real(C3); c0_est(N+1+5)=imag(C3);

elseif input == 41 %regular ellipse
    C1 = 0.3;
    c0_est(3)=real(C1); c0_est(N+1+3)=imag(C1);

elseif input == 51 %irregular pentagon
    C1 = (1/10)+(3/20)*1i; 
    c0_est(3)=real(C1); c0_est(N+1+3)=imag(C1); 
    C4 = (1/10)+(1/20)*1i;
    c0_est(6)=real(C4); c0_est(N+1+6)=imag(C4);

elseif input == 71 %dumbell
    C1 = (-7/10);
    c0_est(3)=real(C1); c0_est(N+1+3)=imag(C1);
    C3 = -0.25;
    c0_est(5)=real(C3); c0_est(N+1+5)=imag(C3);

elseif input == 91 %three pronged object
    C2 = -49/100;
    c0_est(4)=real(C2); c0_est(N+1+4)=imag(C2);
    C5 = -17/100;
    c0_est(7)=real(C5); c0_est(N+1+7)=imag(C5);
    C8 = -3/40;
    c0_est(10)=real(C8); c0_est(N+1+10)=imag(C8);
    C11 = -27/1000;
    c0_est(13)=real(C11); c0_est(N+1+13)=imag(C11);
    C14 = -3/500;
    c0_est(16)=real(C14); c0_est(N+1+16)=imag(C14);
end

% --- Scaling the shape ---
if input~=11 && 00
Pts=1000; Theta=linspace(0,2*pi,Pts); Zeta = exp(1i*Theta); 

z1=c0_est(1)*Zeta; x1 = real(z1); y1=imag(z1); %circle of area pi
z2=c0_est(1)*Zeta; %unscaled shape
    for kk=2:N+1
    z2=z2+(c0_est(kk)+1i*c0_est(N+1+kk))*Zeta.^(-(kk-2)); 
    end
x2 = real(z2); y2=imag(z2); %x and y of unscaled shape
A1 = polyarea(x1,y1); A2 = polyarea(x2,y2); %areas of the two shapes
c0_est = sqrt(A1/A2).*c0_est; %scaled shape
end
end