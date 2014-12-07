function [price, lattice] = EuroCRR(S,K,r,T,sigma,q,N,IsCall,IsAmer)	
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
		% Example:
			% i=0; intrisic = S * u^0*d^(N) = lowest value final node.
			% i=1; intrisic = S * u^1*d^(N-1) = second from lowest final node.
			% i=2; intrisic = S * u^2*d^(N-2) = third from lowest final node.
			% ...
			% u^(N)*d^0 = highest value final node.
		% Why is the lattice upside down?  We are assigning the lowest intrisic...
		% value to the top node in the lattice.  It is arbitrary, but makes it a little less intuitive.
		% So lattice(1,N+1) = lowest value final node.
		% And lattice(N+1, N+1) = highest value final node.
		S_t = S*(u^i)*(d^(N-i));
		lattice(i+1, N+1) = Intrinsic(S_t);
	end

	for j = N-1:-1:0 % Moving backward in time from the second to last time step because we already calculated the final time step.
		for i = 0:j % Moving up and down at a single time.
			% First iteration is lattice(0+N-1, N-1+1) = lattice(N-1, N).
			% This is the second to last timestep, and the highest value at that time.
			% You might think we would use the N slot which is 1 above the N+1 slot,
			% But our tree is not symetrical visually, so we are actually going to move up two slots.
			% In other words, by the final iteration, we want to get back to lattice(1,1)
			% So we move up by two slots on each iteration.
			% Visually, the tree sort of looks like this:
				%        - - - - - 
				%          - - - -
				%            - - -
				%              - -
				%                -
			% The lowest slot to the right is lattice(N+1,N+1).  Moving backwards in time,
			% and up one slot in the tree, we get to lattice(N-1, N).
			% And the final iteration is lattice(1,1)
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