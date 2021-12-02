
%This file reports some functions that can be used for the analysis of the
%linear response. 
%This software has been produced for the H2020 project TiPES and is distributed under the GNU licence agreement
%by Valerio Lucarini
%email: v.lucarini@reading.ac.uk

%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Green,response]=green(T,obs,f)

%This script gives as output the Green function computed for an
%observable and the prediction of the response. We assume that the system is perturbed according by a focing having 
% as time modulation a Heaviside distribution, as discussed in the attached document. 
% We consider an ensemble of N simulations. 
% Each observable is measured for a time frame T (given as input). We
% assume that T is equispaced (otherwise please use first MATLAB functions
% like interp.). We give as input the modulating function we want to use
% for prediction purposes
% the input T and the input f are column vectors of lenght N. obs is a matrix tensor with
% dimensions (T,N). The output Green has dimension (T-1,M).


Ta=T-T(1); %first element of time is set to zero
Text=[-Ta(end:-1:2);Ta] % time is extended symmetrically to the negative domain

dT=T(2)-T(1);

ensobs=squeeze(mean(obs,2));

Green=diff(ensobs)/dT; 

Greenext=[zeros(size(T(1:end-1)));Green;0]; 
% The domain of the Green function is extended symmetrically to the negative domain of times. 
%The Green function vanishes  in the negative domain by definition.

fext=[zeros(size(T(1:end-1)));f]; % The domain of the forcing function is extended symmetrically to the negative domain of times. 
%The forcing function vanishes in the negative domain.

response=conv(Greenext,fext,'same')*dT;

%plot(Text,response) to see the result

omegamax=pi/dT %Nyquist frequency (frequenct is meant as angular frequency) 
N=1000; %number of subdivision of the frequency vector; needs to be adjusted case by case

domega=omegamax/N

omega=[-omegamax:domega:omegamax]; %frequency vector

susceptibility=zeros(size(omega)); %susceptibility vector

for j=1:max(size(susceptibility));
    
    susceptibility(j)=trapz(Text,Greenext*exp(i*omega(j)*Text))*dT;

end
%here we have computed the frequency dependent susceptibility of the
%system.

realchi=real(susceptibility); 
imagchi=imag(susceptibility);

%real and imaginary part of the susceptibility

%The functions below allow for testing the Kramers-Kronig relations 

%%%%%%%%%%%%%%%%%%%%%%%%%%%


function imchi=kkim(omega,rechi)
%The program inputs are the vector of the frequency
%components and the vector of the real part of the susceptibility
%under examination.
%The two vectors must have the same length 
%and the frequency vector omega must be equispaced. 
%If not, apply MATLAB functions such as interp.
%See the book 
%"Kramers-Kronig Relations in Optical Materials Research"
%by Lucarini, V., Saarinen, J.J., Peiponen, K.-E., Vartiainen, E.M. 
%Springer, Heidelberg, 2005
%for further details.
%The output is the estimate of the imaginary part as obtained
%with K-K relations.


if size(omega,1)>size(omega,2);
omega=omega';
end; if size(rechi,1)>size(rechi,2);
rechi=rechi';
end;
%Here the program rearranges the two vectors so that,
%whichever their initial shape, they become row vectors.

g=size(omega,2);
%Size of the vectors.%

imchi=zeros(size(rechi));
%The output is initialized.

a=zeros(size(rechi));
b=zeros(size(rechi));
%Two vectors for intermediate calculations are initialized

deltaomega=omega(2)-omega(1);
%Here we compute the frequency (or energy) interval

j=1;
beta1=0;
for k=2:g;
b(1)=beta1+rechi(k)*omega(k)/(omega(k)^2-omega(1)^2);
beta1=b(1);
end;
imchi(1)=-2/pi*deltaomega*b(1)*omega(1)^(1);
%First element of the output: the principal part integration
%is computed by excluding the first element of the input

j=g;
alpha1=0;
for k=1:g-1;
a(g)=alpha1+rechi(k)*omega(k)/(omega(k)^2-omega(g)^2);
alpha1=a(g);
end;
imchi(g)=-2/pi*deltaomega*a(g)*omega(g)^(1);
%Last element of the output: the principal part integration
%is computed by excluding the last element of the input.

for j=2:g-1; ;
%Loop on the inner components of the output vector.
alpha1=0;
beta1=0;
for k=1:j-1;
a(j)=alpha1+rechi(k)*omega(k)/(omega(k)^2-omega(j)^2);
alpha1=a(j);
end;
for k=j+1:g;
b(j)=beta1+rechi(k)*omega(k)/(omega(k)^2-omega(j)^2);
beta1=b(j);
end;
imchi(j)=-2/pi*deltaomega*(a(j)+b(j))*omega(j)^(1);
%Last element of the output: the principal part integration
%is computed by excluding the last element of the input
end;

%%%%%%%%%%%%%%%%%

function rechi=kkre(omega,imchi)
%The program inputs are the vector of the frequency
%components, the vector of the imaginary
%part of the susceptibility under examination. 
%The two vectors must have the same length 
%and the frequency vector omega must be equispaced. 
%If not, apply MATLAB functions such as interp.
%See the book 
%"Kramers-Kronig Relations in Optical Materials Research"
%by Lucarini, V., Saarinen, J.J., Peiponen, K.-E., Vartiainen, E.M. 
%Springer, Heidelberg, 2005
%for further details.
%The output is the estimate of the real part as obtained
%with K-K relations.

if size(omega,1)>size(omega,2);
omega=omega';
end; if size(imchi,1)>size(imchi,2);
imchi=imchi';
end;
%Here the program rearranges the two vectors so that,
%whichever their initial shape, they become row vectors.
g=size(omega,2);
%Size of the vectors.%
rechi=zeros(size(imchi));
%The output is initialized.
a=zeros(size(imchi));
b=zeros(size(imchi));
%Two vectors for intermediate calculations are initialized
deltaomega=omega(2)-omega(1);
%Here we compute the frequency (or energy) interval
j=1;
beta1=0;
for k=2:g;
b(1)=beta1+imchi(k)*omega(k)^(1)/(omega(k)^2-omega(1)^2);
beta1=b(1);
end;
rechi(1)=2/pi*deltaomega*b(1)*omega(1)^(-2*0);
%First element of the output: the principal part integration
%is computed by excluding the first element of the input
j=g; 
alpha1=0; 
for k=1:g-1;
a(g)=alpha1+imchi(k)*omega(k)^(1)/(omega(k)^2-omega(g)^2);
alpha1=a(g);
end;
rechi(g)=2/pi*deltaomega*a(g)*omega(g)^(-2*0);
%Last element of the output: the principal part integration
%is computed by excluding the last element of the input
for j=2:g-1; ;
%Loop on the inner components of the output vector.
alpha1=0;
beta1=0;
for k=1:j-1;
a(j)=alpha1+imchi(k)*omega(k)^(1)/(omega(k)^2-omega(j)^2);
alpha1=a(j);
end;
for k=j+1:g;
b(j)=beta1+imchi(k)*omega(k)^(1)/(omega(k)^2-omega(j)^2);
beta1=b(j);
end;
rechi(j)=2/pi*deltaomega*(a(j)+b(j))*omega(j)^(-2*0);
end;
%Last element of the output: the principal part integration
%is computed by excluding the last element of


