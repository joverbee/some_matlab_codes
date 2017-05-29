%example code for the holographic reconstruction method in the TEM
%(c) 2017 Jo Verbeeck
% EMAT, University of Antwerp, Belgium
% jo.verbeeck@uantwerp.be

%from a desired wave, the program calculates a binary grating that will result
%in the far field of the plane of the grating in a holographic
%reconstruction of the desired wave (and its complex conjugate and an
%unwanted central beam)
%the binary grating can be produced with a FIB and put in one of the
%condensor planes to recreate the desired wave (and its unwanted
%sideeffects) in the sample plane in STEM mode.

%different choices of reference waves and desired waves are implemented
%use your imagination to expand and surprise me with new ideas/applications!

close all
clear
kpoints=2*4096; %the sampling

%set up the desired wave with its parameters (expand this as needed)

%desiredwave='vortexsuperposition'
%param1=3;
%param2=-3;

desiredwave='abbcorr'
param1=-1e-3; %negative Cs to correct for condensor lens aberrations, not that it is unrealistic to do this with alpha high, as other aberrations than Cs will enter and the bars will get too close to each other

%desiredwave='vortex'
%param1=1; %simple vortex wave



%choose a reference wave
%refwave='quadratic'; %params see in definition of the ref wave
refwave='tiltedplane';

%choose an illumination wave
%illumwave='plane';
illumwave='aberrated'; % a more realistic wave in the condensor plane with aberrations

flip=0; %can flip black and white in aperture, sometimes one is easier to make mechanical as the other, wave is the same due to Babinet?

alpha=21.4e-3; %the half angle of the final probe in the sample plane [rad] (can be typically varied in the condensor settings for STEM mode, bigger gives smaller probes as long as aberrations allow)

%some parameters of the microscope
lambda=2e-12; %wavelength of the fast electrons [m]
k0=2*pi/lambda; %wavevector of the fast electrons [1/m]

%set up a grid in reciprocal space (the space of the grating, typically the condensor plane)
kmax=k0*10*alpha; %highest transverse k vector that we will take into acount (choose somewhat larger as alfa to avoid aliasing effects and get better sampling in real space)
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

switch lower(refwave)
    case 'tiltedplane'
        %define a tilted plane reference wave to interfere with the desired wave
        gridbars=8; %number of gridbars, higher is typically better but harder to make in the FIB
        dx=gridbars/(k0*alpha); %this determines the distance between the diffracted beams
        dy=0;
        refwave=exp(i*2*pi*(Kx*dx+Ky*dy));
    case 'quadratic'
        %quadratic phase wave is nice as it puts the sidebeams on the optical axis at
        %different defocus. In these planes the unwanted sidebeams and the central beam will be out of focus        
        def=100e-9;%defocus needed to bring the desired wave in focus (the higher this is, the more fine the features in the mask will need to become- hard to make and sampling issues)
        refwave=exp(i*k0*(1/2)*Theta.^2.*def);
    otherwise
        disp('unknown ref wave type, exit');
        break;
end


%define the desired wave in reciprocal space (the Fourier transform of the
%real space desired wave)
outermask=Theta<alpha; %clip any reciprocal components above Kmax, needed as the grating will enforce this maximum anyway


switch lower(desiredwave)
    case 'vortex'
        %a simple vortex with m=param1
        flip=1;
        m=param1; %take - otherwise the pattern is inversed and mechanically not connected
        desire=outermask.*exp(+i*m*PHI);              
    case 'vortexsuperposition'    
        %a superposition of 2 vortices with m=param1 and m=param2       
        m1=param1;
        m2=param2;
        desire=outermask.*(exp(+i*m1*PHI)+exp(+i*m2*PHI));   
    case 'abbcorr'
        %a way to compensate for Cs in a non-corrected microscope! a poor
        %mans aberration corrector, Cs is given as param1
        %allows you to make a finer probe than what you thought possible in
        %a given microscope, but at the expense of the unwanted beams which
        %are broader.
        cs=param1;
        chi=k0*(1/4)*Theta.^4.*cs; %simple aberration function with only Cs, can easily be expanded to the full aberration function as below
        %note that correcting the higher order makes not much sense as the
        %aberrations drift over time and the aperture is fixed.
        %chi=(2*pi/lambda)*(...
        %         A0.*Theta.*cos(PHI-Phi11)+...
        %        (1/2)*Theta.^2.*(A1*cos(2*(PHI-Phi22))+C1)+...
        %        (1/3)*Theta.^3.*(A2*cos(3*(PHI-Phi33))+B2*cos(PHI-Phi31))+...
        %        (1/4)*Theta.^4.*(A3*cos(4*(PHI-Phi44))+S3*cos(2*(PHI-Phi42))+C3)+...
        %        (1/5)*Theta.^5.*(A4*cos(5*(PHI-Phi55))+B4*cos(PHI-Phi51)+D4.*cos(3*(PHI-Phi53)))+...
        %        (1/6)*Theta.^6.*(A5*cos(6*(PHI-Phi66))+R5*cos(4*(PHI-Phi64))+S5*cos(2*(PHI-Phi62))+C5)...
        %        );
        desire=outermask.*exp(i*chi);
    case 'airywave'
        %a special wave that has a bent trajectory 
        
    otherwise        
        disp('unknown desired wave type, exit');
        break;
end

%interfere both reference and desired wave
wave=refwave+desire;
I=abs(wave).^2; %intensity pattern of the interference in the grating plane
thresh=0.5*max(I(:)); %clip at half intensity level
%binarise this pattern
if flip==1
    ap=I<thresh;
else
    ap=I>thresh;    
end
%make circular outside rim (needed?)
ap=ap.*(Theta<alpha); 

%now assume a plane wave illumination of this aperture and go to real space
%with a simple Fourier transform
switch lower(illumwave)
    case 'plane'
        illumination=ones(size(ap));
    case 'aberrated'
        A0=0;       
        Phi11 =0;
        A1 = 0;
        Phi22 = 0;
        C1 = 0;
        A2 = 0;
        Phi33 = 0;
        B2 = 0;
        Phi31 = 0;
        A3 = 0;
        Phi44 = 0;
        S3 = 0;
        Phi42 = 0;
        C3 = 1e-3;
        A4 = 0;
        Phi55 = 0;
        B4 = 0;
        Phi51 = 0;
        D4 = 0;
        Phi53 = 0;
        A5 = 0;
        Phi66 = 0;
        R5 = 0;
        Phi64 = 0;
        S5 = 0;
        Phi62 = 0;
        C5 = 0;
        chi=k0*(...
                 A0.*Theta.*cos(PHI-Phi11)+...
                (1/2)*Theta.^2.*(A1*cos(2*(PHI-Phi22))+C1)+...
                (1/3)*Theta.^3.*(A2*cos(3*(PHI-Phi33))+B2*cos(PHI-Phi31))+...
                (1/4)*Theta.^4.*(A3*cos(4*(PHI-Phi44))+S3*cos(2*(PHI-Phi42))+C3)+...
                (1/5)*Theta.^5.*(A4*cos(5*(PHI-Phi55))+B4*cos(PHI-Phi51)+D4.*cos(3*(PHI-Phi53)))+...
                (1/6)*Theta.^6.*(A5*cos(6*(PHI-Phi66))+R5*cos(4*(PHI-Phi64))+S5*cos(2*(PHI-Phi62))+C5)...
                );
        illumination=exp(i*chi);
    otherwise
        disp('unknown desired wave type, exit');
        break;
end

resultwave=fftshift(fft2(fftshift(ap.*illumination)));
Iresultwave=abs(resultwave).^2;

%OPTIONAL
%convolve intensity with gaussian source profile to simulate finite source
%size effect
%sig=0.3e-10;% source size broadening [m]
%fwhm=2*sqrt(2*log(2))*sig
%gausspeak=exp(-((K*lambda/alpha)/(sqrt(2)*sig)).^2);
%gausspeak=gausspeak/sum(gausspeak(:));%normalise
%Iresultwave=ifftshift(ifft2((fft2(fftshift(Iresultwave)).*fft2(fftshift(gausspeak)))));


%show the continuous interference pattern
figure
imagesc(k,k,I);
title('interference pattern')
xlabel('kx [1/m]')
ylabel('ky [1/m]')

%show binarised mask
figure
imagesc(ap);
title('binary mask')
xlabel('kx [1/m]')
ylabel('ky [1/m]')

%show resulting wave
figure;
subplot(1,2,1);
phase=angle(resultwave);
imagesc(x,x,phase);
title('phase')
xlabel('x [m]')
ylabel('y [m]')

subplot(1,2,2);
imagesc(x,x,abs(resultwave).^2);
title('intensity');
xlabel('x [m]')
ylabel('y [m]')

%create nicer image by clipping out center to prepare BMP to send of to FIB
%note the below uses exportfig (get from matlabcentral) to export eps
%figures that can nicely be used in latex.
width=kpoints/8;
height=kpoints/8;
xstart=round(kpoints/2-width/2);
xstop=round(kpoints/2+width/2);
ystart=round(kpoints/2-height/2);
ystop=round(kpoints/2+height/2);
apclip=ap(xstart:xstop,xstart:xstop); %a clipped aperture
Iclip=Iresultwave(ystart:ystop,xstart:xstop); %a clipped resulting intensity
resultwaveclip=resultwave(ystart:ystop,xstart:xstop); %a clipped resulting wave
rclipx=x(xstart:xstop); %clipped x coordinates
rclipy=x(ystart:ystop); %clipped y coordinates

%dpclippsf=absp(ystart:ystop,xstart:xstop);
%dpclippsfsource=abspsource(ystart:ystop,xstart:xstop);

h1=figure
imagesc(apclip)
colormap('gray');
axis image;
axis off
%title('binarised aperture')
imwrite(apclip,'aperture.bmp','BMP'); %can be imported in FIB

%show the intensity of the probe with a certain grayscaling
h2=figure
imagesc(Iclip,[0, 0.1*max(Iclip(:))]);
axis image;
axis off
colormap('gray')

%and another one, oversaturated
h3=figure
imagesc(Iclip,[0, 0.005*max(Iclip(:))]);
axis image;
axis off
colormap('gray')

%and one with color coding for the phase using the nice phplot routine
h4=figure
phplot(resultwaveclip,Iclip);
axis image;
axis off

%and a line trace through the center
plot(rclipx,Iclip(round(height/2),:),'b')
title('linescan through intensity')

%export to EPS to be used in Latex (need to download the exportfig package
%from matlabcentral and place the files in the same folder)
%use in a publication and thank me for the code ;)
exportfig(h1,'aperture.eps'); 
exportfig(h2,'probeintensity.eps')
exportfig(h3,'probeintensityoversat.eps')
exportfig(h3,'probeintensitywithphase.eps')

%loop through focus to see evolution of the wave near the focal point

