FIN 5601 - Financial Technology - Fall 2014
===========================================


1.b) CRR Method
---------------

	>> EuroCRR(50,50,.05,1,.3,.01,100,1)

	ans =

		6.7936


2.b) CRR Oscillatory Convergence
--------------------------------

	>> ConvergenceCRR(50,50,.05,1,.3,.01,50,1)
	% Oscillates above and below.
![ConvergenceCRR](https://raw.githubusercontent.com/matthewfieger/fin_5601/master/ConvergenceCRR.png)

 
3.a) CRR Convergence - Even
---------------------------

	>> ConvergenceCRREven(50,50,.05,1,.3,.01,50,1)
	%Converges smoothly from below.
![ConvergenceCRREven](https://raw.githubusercontent.com/matthewfieger/fin_5601/master/ConvergenceCRREven.png)
 
3.b) CRR Convergence - Odd
--------------------------

	>> ConvergenceCRROdd(50,50,.05,1,.3,.01,50,1)
	% Converges smoothly from above.
![ConvergenceCRROdd](https://raw.githubusercontent.com/matthewfieger/fin_5601/master/ConvergenceCRROdd.png)

 
5) Generalized Binomial Methods
-------------------------------

	>> Binomial(50,47,.05,1,.3,.01,100,1,0,'EQP')

	ans =

		8.3236

	>> Binomial(50,47,.05,1,.3,.01,100,1,0,'TIAN')

	ans =

		8.3186

	>> Binomial(50,47,.05,1,.3,.01,100,1,0,'CRR')

	ans =

		8.3217

	>> Binomial(50,47,.05,1,.3,.01,100,1,0,'LR')

	ans =

		8.3067

References
----------
* Brandimarte, Paolo. Numerical Methods in Finance and Economics: A MATLAB-based Introduction. Hoboken, NJ: Wiley Interscience, 2006. Print.

* Haug, Espen, Gaarder. The Complete Guide to Option Pricing Formulas. New York, NY: McGraw-Hill, 2006. Print.


