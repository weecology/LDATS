now doing with a constant geometric series (ultimate, not penultimate)
this is all done with 1000 iterations, so mondo salt grains for the time being

6 temps

2^2
0.896 0.915 0.894 0.898 0.870

2^6
0.744 0.680 0.784 0.908 0.950

2^10
0.557 0.433 0.827 0.963 0.991

2^15
0.325 0.658 0.953 0.993 0.999

2^20
0.322 0.825 0.984 1.000 1.000


4 temps

2^6
0.594 0.456 0.862



3 temps

2^4
0.567 0.639

2^6
0.457 0.676

2^8
0.319 0.736

2^10
0.144 0.903



2 temps
2^2
0.616

2^4
0.311

2^4.5
0.41

2^5
0.162

2^6
0.121


2^10
0.098


it really looks like we don't want/need > 2 temps./ that we can't get
constant swap rates with more



trying now with 3 temps 2^4 1e5 iterations

0.54994 0.60110


also like oh yeah totes look at this with other nchangepoints
it could be that with fewer change points we need fewer temps


ok with 3 change points back to 1e3 nit for now


3 temps 
2^4
0.268 0.631

2^6
0.177 0.810

2 temps
2^6
0.201


now with 4 change points and 1e4 nit
2 temps 
2^6
0.2304

3 temps 
2^6
0.2749 0.8288


5 change points and 1 e4

2 temps 
2^6
0.2557

3 temps 
2^6
0.3111 0.8253

