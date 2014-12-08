function ConvergenceCRREven(S,K,r,T,sigma,q,N,IsCall)
	% Calculate the continuous time price using the BSM Model
	black_scholes_price = BSPrice(S,K,T,r,sigma,q,IsCall);
	LatticeC = zeros(1,N/2);
	% Calculate the CRR discrete approximation for various timesteps
	for i = (1:N/2)
		LatticeC(i) = EuroCRR(S,K,r,T,sigma,q,i*2-1,IsCall);
	end
	vector = (1:N/2)*2;
	plot(vector, ones(1,N/2)*black_scholes_price);
	hold on;
	plot(vector, LatticeC);