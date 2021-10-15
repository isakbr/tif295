function [lambda,peak] = get_peaks(data)
%Input data is the intensity of a set of measurements
%for all the particles

%Number of particles
nr=length(data(1,1,:));

%Number of measurements
meas=length(data(1,:,1));

%Number of wavelengths
int_nr=length(data(:,1,1));

%Percentile
%  Peaks=prctile(data,[16 50 84],1);
%
 [peak,lambda] = max(data);
 
end
