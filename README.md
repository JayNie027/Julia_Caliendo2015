# Julia_Caliendo2015

  This repo include Julia code to solve general equilibrium model with multi-countries, multi-sectors, and sectoral linkage. The details of model can be obatin from Caliendo and Parro(2015).

  In the program, I set N = 2 (2 countries) and J = 2 (2 sectors). All parameters are symmetric across countries and sectors, e.g. production coeffcients, value added shares, producticities, and etc. There are couple of differences to Caliendo (2015). 
  
  * I assume trade are balanced. One can change that by adding trade deficit in eqlibrium conditions.
  * In their paper, they use hat algebra so there is no need to calibrate full parameters such as bilateral trade cost and sectoral productivities. In this Julia program, you need to calibrate those parameters since it not served for hat algebra.


If you see anything wired. Please feel free to let me know. 
