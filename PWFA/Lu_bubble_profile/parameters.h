const double clight=2.9979e+8;//m/s
const double pi=3.1415926535897931;
const double l0=1.e-6;// m ==1mum
const double t0=l0/clight;// s ~ 3.33 fs
const double me=0.510998928;//MeV
const double qe=1.602176565e-19;// Coulomb
const double re=2.8179403267e-15;//electron classical radius m
const double pC=1.e-12;//pico

double n0  ;//cm^-3
double wp  ;//plasma frequency rad/s
double kp  ;//omegap/c
double lp  ;// plasma length
double Sz_w;// sigma Z gaussian beam witness
double Sz_d;// sigma Z gaussian beam driver
double Sr_d;// sigma R gaussian beam driver
double tm_w;// mean position gaussiam beam witness
double Qd  ;// charge driver
double Qw  ;// charge witness
