% process_MP

% track 
kernel = 0;
pxtom = 1;
fs = 1;
searchrad = 300;
[tracks0,tracklength0] = ptvProcess2(centers,kernel,pxtom,fs,searchrad); %,angles);
save('tracks.mat','tracks0','tracklength0','kernel')