function [price, lattice] = Binomial(S,K,r,T,sigma,q,N,IsCall,IsAmer,Method)	
	deltaT = T/N; % Discrete time step.
	if strcmp(Method,'EQP') % compare two strings
		% EQP Specification
		u=exp((r-q-(sigma^2)/2)*deltaT + sigma*sqrt(deltaT));
		d=exp((r-q-(sigma^2)/2)*deltaT - sigma*sqrt(deltaT));
		p=0.5;
	elseif strcmp(Method,'TIAN') % compare two strings
		% Strike exactly on node - TIAN
		u = exp(sigma * sqrt(deltaT)); % Up movement factor.
		d = 1/u; % Down movement factor.
		j=ceil((log(K/S)- N*log(d))/(log(u/d)));
		tilt = (K/(S*(u^j)*(d^(N-j))))^(1/N);
		u=u*tilt; % tilt the tree
		d=d*tilt;
		p=(exp((r-q)*deltaT) - d)/(u-d);
	elseif strcmp(Method,'LR')
		% Leisen-Reimer Specification
		d1 = (log(S/K) + (r-q+(sigma^2)/2)*T) / (sigma*sqrt(T));
		d2 = (log(S/K) + (r-q-(sigma^2)/2)*T) / (sigma*sqrt(T));

		% Preizer-Pratt inversion method 1
		h = @(x) 0.5+sign(x)*(0.25-0.25*(exp( -1*((x/(N+1/3))^2)*(N+1/6) ) ) )^0.5;

		% Preizer-Pratt inversion method 2
		% h = @(x) 0.5+sign(x)*(0.25-0.25*(exp( -1*((x/(N+1/3+0.1/(N+1)))^2)*(N+1/6) ) ) )^0.5;

		p = h(d2);
		u = exp((r-q)*deltaT)*(h(d1)/h(d2));
		d = (exp((r-q)*deltaT) - p*u) / (1-p);
	else
		% CRR Specification
		u = exp(sigma * sqrt(deltaT)); % Up movement factor.
		d = 1/u; % Down movement factor.
		p = (exp((r-q)*deltaT) - d)/(u-d); % Probability of up movement.
	end

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
				IntrinsicVal = Intrinsic(S*(u^(i))*(d^(j-i)));
				lattice(i+1,j+1) = max(IntrinsicVal,OptionVal);
			else
				lattice(i+1, j+1) = OptionVal;
			end
		end
	end
	price = lattice(1,1);