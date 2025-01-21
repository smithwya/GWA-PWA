Comments to the input file:
#the type of kmatrix for each amplitude can be: "kmat-nominal", "kmat-CDD"
#the type of rhoN for each amplitude can be: "rhoN-nominal", "rhoN-Qmodel"
#"Fit", "Polesearch" and "Plot" are mutually exclusive actions
#"polesearchzero(-n)" means that the Polesearch algorithm will consider numbers < 10^-7 as zeros 
#grid(double resx, double redx, int repts, double imsx, double imdx, int impts)
#in fitparams for each call we specify the parameters SetMaxFunctionCalls, SetMaxIterations, SetTolerance, SetPrintLevel;

Remember:
1- CDD Kmat has positive coeffs whose positions in the list of fit parameters is hardcoded counted; the same is for the Kmat poles.
