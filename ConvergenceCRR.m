function ConvergenceCRR(S,K,r,T,sigma,q,N,IsCall)
	% Calculate the continuous time price using the BSM Model
	black_scholes_price = BSPrice(S,K,T,r,sigma,q,IsCall);
	LatticeC = zeros(1,N);
	% Calculate the CRR discrete approximation for various timesteps
	for i = (1:N)
		LatticeC(i) = EuroCRR(S,K,r,T,sigma,q,i,IsCall);
	end
	plot(1:N, ones(1,N)*black_scholes_price);
	hold on;
	plot(1:N, LatticeC);