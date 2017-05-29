%example code to demonstrate the propagation of all kinds of weird beams
%that could be made inside a TEM (and originate in original work in optics)
%(c) 2017 Jo Verbeeck
% EMAT, University of Antwerp, Belgium
% jo.verbeeck@uantwerp.be

%use your imagination to expand and surprise me with new ideas/applications!

close all
clear
kpoints=512; %the sampling

%choose the type of wave (or better, expand with new types)
wavetype='circular'
%wavetype='circularvortex'
%wavetype='bessel'
%wavetype='airy'
%wavetype='axicon'
%wavetype='helicon'
%wavetype='snake'

alpha=21.4e-3; %the half angle of the final probe in the sample plane [rad] (can be typically varied in the condensor settings for STEM mode, bigger gives smaller probes as long as aberrations allow)

%some parameters of the microscope
lambda=2e-12; %wavelength of the fast electrons [m]
k0=2*pi/lambda; %wavevector of the fast electrons [1/m]

%set up a grid in reciprocal space (the space of the grating, typically the condensor plane)
oversample=4; %the higher this is, the better the sampling in real space, but the worse in k space
kmax=k0*oversample*alpha; %highest transverse k vector that we will take into acount
k=linspace(-kmax,kmax,kpoints);
[Kx,Ky]=meshgrid(k,k); %create 2D grid of k points
K=sqrt(Kx.^2+Ky.^2);
Theta=K/k0; %and convert to angles
PHI=atan2(Kx,Ky); %angle in the 2D plane

%set up a grid in real space (the space of the final probe, typically the sample plane)
xmax=kpoints/kmax;
x=linspace(-xmax,xmax,kpoints);
[X,Y]=meshgrid(x,x); %and create a 2D grid with it
R=sqrt(X.^2+Y.^2);
rdisc=2/(k0*alpha);

outermask=Theta<alpha; %clip any reciprocal components above Kmax, needed as the grating will enforce this maximum anyway

%create a wave in the condensor plane
switch lower(wavetype)
    case 'circular'
        %conventional idealised (no aberrations) STEM beam from a circular aperture
        %easy to add aberrations as well
        wave=outermask;
    case 'circularvortex'
        %conventional idealised (no aberrations) STEM beam from a circular aperture
        %easy to add aberrations as well
        m=1;
        wave=outermask.*exp(1i*m*PHI);
    case 'airy'
        %approximation to an Airy wave (no amplitude modulation)
        a=2*(1/(k0*alpha)^3); %strength factor determining the focal range and size of the probe
        wave=outermask.*exp(1i*2*pi*a*((Kx+Ky).^3+(Kx-Ky).^3));
    case 'bessel'
        %an approximate Bessel beam (a ring of finite width in k space, mathematically should be infinitely thin ring)        
        kouter=k0*alpha;
        kinner=0.9*kouter;
        m=0; %topological charge
        wave=outermask.*((K>kinner).*(K<kouter)).*exp(i*m*PHI); %a ring like illumination
    case 'snake'
        %try to make a snake bessel wave
        kouter1=k0*alpha;
        kinner1=0.95*kouter1;
        kouter2=0.8*k0*alpha;
        kinner2=0.95*kouter2;
        PHIp=PHI+pi/2;
        m1=2; %topological charge
        m2=1; %topological charge
        bessel1=(K>kinner1).*(K<kouter1).*(exp(i*m1*PHIp)+exp(-i*m1*PHIp));
        bessel2=(K>kinner2).*(K<kouter2).*(exp(i*m2*PHIp)+exp(-i*m2*PHIp));
        wave=outermask.*(bessel1+bessel2); %a ring like illumination        
    case 'axicon'
        %quasi Bessel produced with an Axicon (doesn't seem to do what I
        %want
        a=2/(k0*alpha); %strength of the axicon
        wave=outermask.*exp(1i*2*pi*a*K);
    case 'helicon'        
        %A helicon wave can be generated as a superposition of 2 bessel
        %waves, rotation speed of the helix is dependent on (k1^2-k2^2)/(m1-m2)
        %see Patterson paper        
        kouter1=k0*alpha;
        kinner1=0.95*kouter1;
        kouter2=0.8*k0*alpha;
        kinner2=0.95*kouter2;        
        m1=1; %topological charge
        m2=-1; %topological charge
        bessel1=(K>kinner1).*(K<kouter1).*(exp(i*m1*PHI));
        bessel2=(K>kinner2).*(K<kouter2).*(exp(i*m2*PHI));
        wave=outermask.*(bessel1+bessel2); 
    otherwise
        disp('unknown wave type, exit');
        break;
end

%show the phase of the incoming probe
figure
phplot(wave,abs(wave).^2);
axis image;
title('wave k space')
xlabel('kx [1/m]')
ylabel('ky [1/m]')


%now calculate the distribution near focus in real space
%shift up and down the focal plane using a defocus term in the condensor
%plane (Fresnel propagator method). This only works for a 'reasonable'
%defocus range depending heavily on sampling
focusrange=100e-9;
focussteps=100;
def=linspace(-focusrange/2,focusrange/2,focussteps);

I=zeros(kpoints,kpoints,focussteps); %preallocate memory for I (speed)
figure
for k=1:focussteps   
    defocus=exp(1i*k0*(1/2)*Theta.^2.*def(k));          
    resultwave=fftshift(fft2(fftshift(wave.*defocus))); %Fourier transforming brings us to the observation plane
    I(:,:,k)=abs(resultwave).^2; %store it in a 3D dataset
    imagesc(x,x,I(:,:,k));
    title('intensity in plane')
    pause(0.1); %to force redrawing of figure while in the loop
end

%make cross sections through this 3D cube
figure
slice(x,x,def,I,[],[0],[]); %create a set of slices defined by cuts to the points defined by the last 3 parameters e.g. [],[0],[] means a cut through y=0 plane.
title('slice through 3D intensity')
xlabel('x');
ylabel('y');
zlabel('z');
shading flat
view(0,0)
