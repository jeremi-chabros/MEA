S = sparseDownSample(mSpikes, 720 * 10, 'sum'); 
Stot = sum(S, 2); 
plot(Stot) 