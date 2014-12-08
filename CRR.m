function [price, lattice] = CRR(S,K,r,T,sigma,q,N,IsCall,IsAmer)	
	deltaT = T/N; % Discrete time step.
	u = exp(sigma * sqrt(deltaT)); % Up movement factor.
	d = 1/u; % Down movement factor.
	p = (exp((r-q)*deltaT) - d)/(u-d); % Probability of up movement.
	lattice = zeros(N+1, N+1); % Empty lattice for our values. Must add 1 for one-based indexing.

	if IsCall
		Intrinsic = @(S) max(0,S-K);
	else
		Intrinsic = @(S) max(0,K-S);
	end

	% Our lattice is lattice(value,time) or lattice(y,x) if thinking visually.
	for i = 0:N
		% Compute intrinsic value at maturity
		S_t = S*(u^i)*(d^(N-i));
		lattice(i+1, N+1) = Intrinsic(S_t);
	end

	for j = N-1:-1:0 % Moving backward in time from the second to last time step because we already calculated the final time step.
		for i = 0:j % Moving up and down at a single time.
			OptionVal = exp(-r*deltaT) * (p * lattice(i+2,j+2) + (1-p) * lattice(i+1,j+2));
			if IsAmer
				IntrinsicVal = Intrinsic(S*(u^(i))*(d^(j-i))); % Is this j-i?
				lattice(i+1,j+1) = max(IntrinsicVal,OptionVal);
			else
				lattice(i+1, j+1) = OptionVal;
			end
		end
	end
	price = lattice(1,1);