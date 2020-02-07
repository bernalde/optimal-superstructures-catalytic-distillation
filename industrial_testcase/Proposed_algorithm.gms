$ontext
================================================================================
Caso de estudio industrial

          David Esteban Bernal Neira-David Alejandro Liñan Romero
                                 2019
================================================================================
$offtext
*-------------------------------------------------------------------------------
*                                Sección 1
*           Conjuntos para definir operacion en estado estable
*-------------------------------------------------------------------------------
*Se usa solo el elemento inicial con el punto de colocacion inicial: estado estable
*Esto es equivalente a quitar j y N de la formualcion
sets j "1 punto de colocacion (0)" /1/
     N "1 elemento finito" /1/;
*-------------------------------------------------------------------------------
*                                Sección 2
*       Conjuntos, variables, parámetros y ecuaciones principales del sistema
*-------------------------------------------------------------------------------
*Conjuntos
set comp "lista de componentes que intervienen en el sistema" /iButene, Ethanol, nButene, ETBE/;
set Net "Todas las etapas de la columna reales y sobrantes incluyendo con y rev" /1*45/;
*Alias para Net, que se usa en las restricciones binarias:
alias(Net,Net1);
set F "1 representa etanol y 2 representa butenos" /1,2/;

*Variables principales
positive variables
L(N,j,Net)  "Flujo de líquido [mol/min]"
V(N,j,Net)  "Flujo de vapor [mol/min]"
x(N,j,comp,Net)     "Porcentaje molar en el líqudio [%]"
y(N,j,comp,Net)     "Porcentaje molar en el vapor [%]"
Temp(N,j,Net)       "Temperatura de operación [K]"
P(N,j,Net)  "Presión por etapa [bar]"
Z(N,j,Net)  "Coeficiente de compresibilidad [-]"
RR(N,j)      "Relación molar de reflujo [-]"
Qc(N,j)      "Carga térmica del condensador [kJ/min]"
Qr(N,j)      "Carga térmica del rehervidor [kJ/min]"
BR(N,j)       "Boil up [-]"
;


*Parámetros hidráulicos
parameter
da      "Diámetro de los agujeros [m]"  /2E-3/
ep      "Espesor del plato [m]" /0.002/
pitch   "Distancia entre agujeros [m]"  /0.009/
Sfactor "Factor de seguridad altura de la columna [-]" /0.15/
poro  "Porosidad del plato [-]"
K0    "Coeficiente de orificio [-]"
;
poro=0.907*sqr(da/pitch);
K0=(880.6-(67.7*da/ep)+(7.32*((da/ep)**2))-(0.338*((da/ep)**3)))*1E-3;

*Variables hidraulicas
positive variable D "Diámetro de la columna [m]";
positive variable hw      "Weir height [m]";
positive variable  HS      "Altura de cada plato [m]";
positive variable Htotal "Altura total de la columna [m]";
positive variable At "Area activa [m2]";
positive variable Ad  "Area de derramadero [m2]";
positive variable Lw "Weir length [m]";
positive variable A0 "Area agujerada [m2]";

*Ecuaciones hidraulicas
Equation EqHwmin;
EqHwmin.. hw=g=0.05*HS;
Equation EqHwmax;
EqHwmax.. hw=l=HS/3;
equation EqAt;
EqAt.. At=e=sqr(D/2)*(pi-(1.854590-0.96));
equation EqAd;
EqAd.. Ad=e=sqr(D/2)*(0.5*(1.854590-0.96));
equation EqLw;
EqLw.. Lw=e=0.8*D;
equation EqA0;
EqA0.. A0=e=At*poro;







*Alimentacion  butenos
parameter
FB "Flujo de alimentación de butenos [mol/min]" /29458.33333/
zb(N,j,comp) "Porcentaje molar en la alimentación de butenos";
zb(N,j,'iButene')=40;
zb(N,j,'nButene')=100-zb(N,j,'iButene');
zb(N,j,'Ethanol')=0;
zb(N,j,'ETBE')=0;

*Alimentacion etanol


parameter
FE  "Flujo de alimentación de etanol [mol-h]"   /11783.33333/
ze(comp)        "Porcentaje molar en la alimentación de etanol"
/
iButene 0
Ethanol 100
nButene 0
ETBE 0
/
;
*Parametros de operacion
parameter
Pop     "Presión de operación condensador [bar]"        /7.6/
TaliB   "Temperatura de alimentación de butenos [K]"    /323/
TaliE   "Temperatura de alimentación etanol [K]"        /342.38/
xBetbe  "Composición molar de ETBE en fondos deseada"   /99.1/
cR      "Constante de los gases [m3*bar/K*mol]" /0.00008314/
;
*-------------------------------------------------------------------------------
*                                Sección 3
*                    Parametro de conversion de unidades
*-------------------------------------------------------------------------------
parameter
hora    "Si estamos en análisis por minutos u hora [s]" /60/;

*-------------------------------------------------------------------------------
*                                Sección 4
*                          Restricciones de pureza
*-------------------------------------------------------------------------------

equations pureza0(Net);
pureza0(Net)$(ord(Net) eq card(Net)).. x('1','1','ETBE',Net)=g=xBetbe;

*-------------------------------------------------------------------------------
*                                Sección 5
*    Cálculo de presiones de saturación por medio de la ecuación de Antoine
*-------------------------------------------------------------------------------
*Constantes de la ecuación de Antoine expandida
parameters
C1a(comp)
/
iButene 66.4970745
Ethanol 61.7910745
nButene 40.3230745
ETBE    52.67507454
/
C2a(comp)
/
iButene -4634.1
Ethanol -7122.3
nButene -4019.2
ETBE    -5820.2
/
C3a(comp)
/
iButene 0
Ethanol 0
nButene 0
ETBE    0
/
C4a(comp)
/
iButene 0
Ethanol 0
nButene 0
ETBE    0
/
C5a(comp)
/
iButene -8.9575
Ethanol -7.1424
nButene -4.5229
ETBE    -6.1343
/
C6a(comp)
/
iButene 1.3413E-5
Ethanol 2.8853E-6
nButene 4.8833E-17
ETBE    2.1405E-17
/
C7a(comp)
/
iButene 2
Ethanol 2
nButene 6
ETBE    6
/
;
positive variables Psat(N,j,comp,Net) presión de saturación (bar);
equations EqPsat(N,j,comp,Net);
EqPsat(N,j,comp,Net).. Psat(N,j,comp,Net)=e=exp( C1a(comp) + (C2a(comp)/(Temp(N,j,Net)+C3a(comp))) + (C4a(comp)*Temp(N,j,Net)) + (C5a(comp)*log(Temp(N,j,Net)) + (C6a(comp)*(Temp(N,j,Net)**C7a(comp)))) );

*-------------------------------------------------------------------------------
*                                Sección 6
*     Cálculo de densidades de líquido por medio de la ecuación IK-CAPI
*     Cálculo de densidades de líquido por medio de la ecuación DIPPR crítica
*     Cálculo de densidades de gas por medio de ecuación de gas ideal corregida
*-------------------------------------------------------------------------------
*Constantes de la ecuación DIPPR
parameters
MW(comp) "Peso molecular [kg/kmol]"
/
iButene 56.10752
Ethanol 46.06904
nButene 56.10752
ETBE    102.17656
/
Tcrit(comp) "Temperatura crítica [K]"
/
iButene 417.9
Ethanol 516.2
nButene 419.6
ETBE    509.4
/
Pcrit(comp) "Presión crítica [bar]"
/
iButene 38.98675
Ethanol 60.35675
nButene 39.18675
ETBE    28.32675
/
C1rh(comp)
/
iButene        8.9711123119
Ethanol        -2.932961888E-2
nButene        5.956235579
ETBE           -1.323678817E-1
/
C2rh(comp)
/
iButene        0
Ethanol        6.9361857406E-4
nButene        0
ETBE           2.1486345729E-3
/
C3rh(comp)
/
iButene 0
Ethanol -1.962897037E-6
nButene 0
ETBE    -6.092181735E-6
/
C4rh(comp)
/
iButene 0
Ethanol 2.089632106E-9
nButene 0
ETBE    6.4627035532E-9
/
C5rh(comp)
/
iButene 0
Ethanol 0
nButene 0
ETBE    0
/
C6rh(comp)
/
iButene -1.4666609E-10
Ethanol 0
nButene -9.3717935E-11
ETBE    0
/
C7rh(comp)
/
iButene 1.286186216E-12
Ethanol 0
nButene 8.150339357E-13
ETBE    0
/
C8rh(comp)
/
iButene -4.33826109E-15
Ethanol 0
nButene -2.72421122E-15
ETBE    0
/
C9rh(comp)
/
iButene 6.619652613E-18
Ethanol 0
nButene 4.115761136E-18
ETBE    0
/
C10rh(comp)
/
iButene -3.8362103001E-21
Ethanol 0
nButene -2.3593237507E-21
ETBE    0
/
C1r(comp)
/
iButene        1.1446
Ethanol        1.6288
nButene        1.0877
ETBE        0.66333
/
C2r(comp)
/
iButene        0.2724
Ethanol        0.27469
nButene        2.6454E-01
ETBE        2.6135E-01
/
C3r(comp)
/
iButene        0.28172
Ethanol        0.23178
nButene        0.2843
ETBE        0.28571
/
C4r(comp)
/
iButene 0
Ethanol 0
nButene 0
ETBE    0
/
;
positive variable Tcritm(N,j,Net);

equation EqTcritm(N,j,Net);
EqTcritm(N,j,Net).. Tcritm(N,j,Net) =e= (sqr(sum(comp,(x(N,j,comp,Net)/100)*Tcrit(comp)/(Pcrit(comp)**0.5))))/(sum(comp,(x(N,j,comp,Net)/100)*Tcrit(comp)/Pcrit(comp)));
positive variables rho(N,j,comp,Net) "Densidad molar por componente de líquido [mol/m^3]";
equation Eqrho(N,j,comp,Net);
Eqrho(N,j,comp,Net).. rho(N,j,comp,Net)=e=( C1r(comp)/(C2r(comp)**(1+((1-(Temp(N,j,Net)/Tcritm(N,j,Net)))**C4r(comp)))) )*1000;

positive variable rhoV(N,j,Net) "Densidad molar de vapor [mol/m^3]";
equation EqurhoV(N,j,Net);
EqurhoV(N,j,Net).. rhoV(N,j,Net)=e=P(N,j,Net)/(0.00008314*Temp(N,j,Net)*(Z(N,j,Net)));

*-------------------------------------------------------------------------------
*                                Sección 7
*     Cálculo de tensión superficial por medio de la ecuación DIPPR crítica
*-------------------------------------------------------------------------------
*Constantes de la ecuación DIPPR
parameters
C1sig(comp)
/
iButene        0.05544
Ethanol        0.03764
nButene        0.055945
ETBE        0.071885
/
C2sig(comp)
/
iButene        1.2453
Ethanol        -2.157E-5
nButene        1.2402
ETBE        2.1204
/
C3sig(comp)
/
iButene        0.0
Ethanol        1.025E-7
nButene        0
ETBE        -1.5583
/
C4sig(comp)
/
iButene 0
Ethanol 0
nButene 0
ETBE    0.76657
/
;
positive variables sigma(N,j,Net) "Tensión superficial líquido vapor [N/m]";
equation Eqsigma(N,j,Net);
Eqsigma(N,j,Net).. sigma(N,j,Net)=e=sum(comp,(x(N,j,comp,Net)/100)*C1sig(comp)*(1-(Temp(N,j,Net)/Tcritm(N,j,Net)))**(C2sig(comp)+C3sig(comp)*(Temp(N,j,Net)/Tcritm(N,j,Net))+C4sig(comp)*((Temp(N,j,Net)/Tcritm(N,j,Net)))**2));
*-------------------------------------------------------------------------------
*                                Sección 8
*          Cálculo de coeficientes de actividad por medio del modelo NRTL
*-------------------------------------------------------------------------------
table a_nrtl(comp,comp) Parámetro a de NRTL
                       iButene            Ethanol            nButene            ETBE
iButene                0.0                0.0                0.0                0.0
Ethanol                0.0                0.0                0.0                0.0
nButene                0.0                0.0                0.0                0.0
ETBE                   0.0                0.0                0.0                0.0
;

table b_nrtl(comp,comp) Parámetro b de NRTL
                iButene                Ethanol            nButene            ETBE
iButene         0.0                    623.5810010        107.526499         219.73407
Ethanol         141.9632130            0.0                164.57256          187.104064
nButene         -93.24546420           595.5299820        0.0                226.373398
ETBE            -172.59152             344.481315         -177.88565         0.0
;

table c_nrtl(comp,comp) Parámetro c de NRTL
                       iButene            Ethanol            nButene            ETBE
iButene                0.0                0.3                0.3                0.3
Ethanol                0.3                0.0                0.3                0.3
nButene                0.3                0.3                0.0                0.3
ETBE                   0.3                0.3                0.3                0.0
;
alias (comp,comp1);
parameter alfa_nrtl(comp,comp);
alfa_nrtl(comp,comp1)$(ord(comp) ne ord(comp1))=c_nrtl(comp,comp1);

*Parámetros G y Tao
variables tao_nrtl(N,j,comp,comp1,Net);
equations Eq_tao_nrtl(N,j,comp,comp1,Net);
Eq_tao_nrtl(N,j,comp,comp1,Net).. tao_nrtl(N,j,comp,comp1,Net)=e=a_nrtl(comp,comp1) + (b_nrtl(comp,comp1)/Temp(N,j,Net));

variables g_nrtl(N,j,comp,comp1,Net);
equations Eq_g_nrtl(N,j,comp,comp1,Net);
Eq_g_nrtl(N,j,comp,comp1,Net).. g_nrtl(N,j,comp,comp1,Net)=e=exp( -alfa_nrtl(comp,comp1)*tao_nrtl(N,j,comp,comp1,Net));

*Coeficiente de actividad (gamma)
alias (comp,comp2,comp3);
variables gamma(N,j,comp,Net);
equations Eqgamma(N,j,comp,Net);
Eqgamma(N,j,comp,Net).. gamma(N,j,comp,Net)=e=
        exp(sum(comp1,x(N,j,comp1,Net)*tao_nrtl(N,j,comp1,comp,Net)*
        g_nrtl(N,j,comp1,comp,Net))/sum(comp1,x(N,j,comp1,Net)*
        g_nrtl(N,j,comp1,comp,Net))+sum(comp1,x(N,j,comp1,Net)*
        g_nrtl(N,j,comp,comp1,Net)/sum(comp2,x(N,j,comp2,Net)*
        g_nrtl(N,j,comp2,comp1,Net))*(tao_nrtl(N,j,comp,comp1,Net)-
        sum(comp2,x(N,j,comp2,Net)*tao_nrtl(N,j,comp2,comp1,Net)*
        g_nrtl(N,j,comp2,comp1,Net))/sum(comp3,x(N,j,comp3,Net)*
        g_nrtl(N,j,comp3,comp1,Net)))));

*-------------------------------------------------------------------------------
*                                Sección 9
*                      Cálculo de reacción química
*-------------------------------------------------------------------------------
Parameter
Nu(comp) "Coeficientes estequiométricos en la reacción"
/
iButene -1
Ethanol -1
nButene 0
ETBE    1
/
;
parameter mcat  /1000/;

variable Ketbe(N,j,Net) "Constante de equilibrio [-]";
$ontext
$offtext
equation EqKetbe(N,j,Net);
EqKetbe(N,j,Net)$((ord(Net) ne card(Net)) and (ord(Net) ne 1)).. Ketbe(N,j,Net) =e=
                exp(10.387+4060.59/(Temp(N,j,Net))
                -2.89055*log(Temp(N,j,Net))
                -0.01915144*Temp(N,j,Net)
                +0.0000528586*power(Temp(N,j,Net),2)
                -0.0000000532977*power(Temp(N,j,Net),3));

positive variable Krate(N,j,Net) "Tasa de avance de reacción [mol/(kg_cat.min)]";
equation EqKrate(N,j,Net);
EqKrate(N,j,Net)$((ord(Net) ne card(Net)) and (ord(Net) ne 1))..  Krate(N,j,Net) =e= 7.41816E15*exp(-60400.0/(8.314*Temp(N,j,Net)))*hora/3600;
positive variable Ka(N,j,Net) "Tasa de adsorción";
equation EqKa(N,j,Net);
EqKa(N,j,Net)$((ord(Net) ne card(Net)) and (ord(Net) ne 1))..  Ka(N,j,Net) =e= exp(-1.0707+1323.1/Temp(N,j,Net));
variable Rx(N,j,Net) "Tasa de reacción [mol/(kg_cat.min)]";
equation EqRx(N,j,Net);
EqRx(N,j,Net)$((ord(Net) ne card(Net)) and (ord(Net) ne 1))..  Rx(N,j,Net)*(power(1+Ka(N,j,Net)*gamma(N,j,'Ethanol',Net)*x(N,j,'Ethanol',Net)/100,3))*Ketbe(N,j,Net) =e=
                        (Krate(N,j,Net)*(gamma(N,j,'Ethanol',Net)*x(N,j,'Ethanol',Net)/100))
                        *((Ketbe(N,j,Net)*gamma(N,j,'iButene',Net)*x(N,j,'iButene',Net)/100*gamma(N,j,'Ethanol',Net)*x(N,j,'Ethanol',Net)/100)
                        -(gamma(N,j,'ETBE',Net)*x(N,j,'ETBE',Net)/100));
*-------------------------------------------------------------------------------
*                                Sección 10
*                           Ecuación de estado (calculo de phi)
*-------------------------------------------------------------------------------
parameter
Omega(comp) "Factor acéntrico [-]"
/
iButene 0.19484
Ethanol 0.643558
nButene 0.184495
ETBE    0.316231
/
TcritSRK(comp) "Temperatura crítica de Soave-Redlich-Kwong [K]"
/
iButene 417.9
Ethanol 514
nButene 419.5
ETBE    509.4
/
mEOS(comp) "Parameter m in EOS"
biEOS(comp) "Parameter bi in EOS"
;
mEOS(comp)=0.48508+1.55171*Omega(comp)-0.15613*sqr(Omega(comp));
biEOS(comp)=0.08664*0.00008314*TcritSRK(comp)/Pcrit(comp);
positive variable alphaEOS(N,j,comp,Net);
equation EqAlphaEOS(N,j,comp,Net);
EqAlphaEOS(N,j,comp,Net).. alphaEOS(N,j,comp,Net) =e= sqr(1+mEOS(comp)*(1-(Temp(N,j,Net)/Tcritm(N,j,Net))**(1/2)));
positive variable aiEOS(N,j,comp,Net);
equation EqaiEOS(N,j,comp,Net);
EqaiEOS(N,j,comp,Net).. aiEOS(N,j,comp,Net) =e= alphaEOS(N,j,comp,Net)*0.42747*(sqr(0.00008314*TcritSRK(comp)))/Pcrit(comp);

positive variable bEOS(N,j,Net);
equation EqbEOS(N,j,Net);
EqbEOS(N,j,Net).. bEOS(N,j,Net) =e= sum(comp,(y(N,j,comp,Net)/100)*biEOS(comp));

positive variable aEOS(N,j,Net);
equation EqaEOS(N,j,Net);
EqaEOS(N,j,Net).. aEOS(N,j,Net) =e= sum(comp,sum(comp1, (y(N,j,comp,Net)/100)*(y(N,j,comp1,Net)/100)*(aiEOS(N,j,comp,Net)*aiEOS(N,j,comp1,Net))**0.5));

equation VaporZ(N,j,Net);
VaporZ(N,j,Net).. (Z(N,j,Net))**3-(Z(N,j,Net))**2+(Z(N,j,Net))
                *((aEOS(N,j,Net)*P(N,j,Net)/((0.00008314*Temp(N,j,Net))**2))
                -(bEOS(N,j,Net)*P(N,j,Net)/(0.00008314*Temp(N,j,Net)))
                -(bEOS(N,j,Net)*P(N,j,Net)/(0.00008314*Temp(N,j,Net)))**2)
                -((aEOS(N,j,Net)*P(N,j,Net)/((0.00008314*Temp(N,j,Net))**2)))
                *(bEOS(N,j,Net)*P(N,j,Net)/(0.00008314*Temp(N,j,Net))) =e= 0;

positive variable phi(N,j,comp,Net);
equation EqPhi(N,j,comp,Net);
EqPhi(N,j,comp,Net).. phi(N,j,comp,Net) =e= exp(((Z(N,j,Net))-1)*biEOS(comp)/bEOS(N,j,Net)
                                        -log((Z(N,j,Net))-bEOS(N,j,Net))
                                        -(aEOS(N,j,Net)/bEOS(N,j,Net))
                                        *(2*((aiEOS(N,j,comp,Net)/aEOS(N,j,Net))**(1/2))
                                        -biEOS(comp)/bEOS(N,j,Net))*log(((Z(N,j,Net))
                                        -bEOS(N,j,Net))/(Z(N,j,Net))));

*-------------------------------------------------------------------------------
*                                Sección 11
*                           Cálculo de entalpías
*-------------------------------------------------------------------------------
*Constantes de Cp (kJ/mol.K) gas ideal
parameters
C1c(comp)
/
iButene        0.016052191
Ethanol        0.00901418
nButene        -0.00299356
ETBE        -0.014651654
/
C2c(comp)
/
iButene        0.000280432
Ethanol        0.000214071
nButene        0.000353198
ETBE        0.000698631
/
C3c(comp)
/
iButene        -0.00000010914988
Ethanol        -0.000000083903472
nButene        -0.00000019904047
ETBE           -0.00000044791741
/
C4c(comp)
/
iButene        0.0000000000090979164
Ethanol        0.0000000000013732704
nButene        0.000000000044631288
ETBE           0.00000000011636811
/
C5c(comp)
/
iButene        0
Ethanol        0
nButene        0
ETBE           0
/
C6c(comp)
/
iButene        0
Ethanol        0
nButene        0
ETBE           0
/
;
parameter
Tref "Temperatura de referencia [K]" /298.15/
Hform(comp) "Entalpía de formación (kJ/mol)"
/
iButene -16.9147
Ethanol -234.963
nButene -0.125604
ETBE        -313.9
/
Tb      "Temperatura de ebullición de los componentes a P=9.5bar [K]"
/
iButene 341.7
Ethanol 421.9
nButene 342.6
ETBE    438.8
/
;
*Entalpía de la fase vapor (kJ/mol) -- Int(CpdT)
variable HVi(N,j,comp,Net),HV(N,j,Net);
equations EqHVi(N,j,comp,Net),EqHV(N,j,Net);
EqHVi(N,j,comp,Net).. HVi(N,j,comp,Net)=e=( (C1c(comp)*(Temp(N,j,Net)-Tref)) + ((C2c(comp)/2)*((Temp(N,j,Net)**2)-(Tref**2)))
                                   + ((C3c(comp)/3)*((Temp(N,j,Net)**3)-(Tref**3))) + ((C4c(comp)/4)*((Temp(N,j,Net)**4)-(Tref**4)))
                                   + ((C5c(comp)/5)*((Temp(N,j,Net)**5)-(Tref**5))) + ((C6c(comp)/6)*((Temp(N,j,Net)**6)-(Tref**6))) + Hform(comp)
                                   + (8.314/1000)*Temp(N,j,Net)*(Z(N,j,Net)-1)+(1+mEOS(comp))*((aEOS(N,j,Net)**0.5)/bEOS(N,j,Net))*log(Z(N,j,Net)/(Z(N,j,Net)+(bEOS(N,j,Net)*P(N,j,Net)/(0.00008314*Temp(N,j,Net))))));

EqHV(N,j,Net).. HV(N,j,Net)=e=sum(comp,HVi(N,j,comp,Net)*y(N,j,comp,Net)/100);

*Constantes de entalpía de vaporización (kJ/mol)
parameter
C1v(comp)
/iButene        32.614
Ethanol        55.789
nButene        33.774
ETBE        45.29
/
C2v(comp)
/iButene        0.38073
Ethanol        0.31245
nButene        0.5107
ETBE        0.27343
/
C3v(comp)
/iButene        0
Ethanol        0
nButene        -0.17304
ETBE        0.21645
/
C4v(comp)
/iButene        0
Ethanol        0
nButene        0.05181
ETBE        -0.11756
/
C5v(comp)
/iButene        0
Ethanol        0
nButene        0
ETBE        0
/
;

*Temperaturas reducidas
parameter Tred(comp);
Tred(comp)=Tb(comp)/Tcrit(comp);


parameter alphaEOSb(comp), aiEOSb(comp);
alphaEOSb(comp)=(1+mEOS(comp)*(1-(Tb(comp)/Tcrit(comp))**(1/2)))**2;
aiEOSb(comp)=alphaEOSb(comp)*0.42747*((0.00008314*TcritSRK(comp))**2)/Pcrit(comp);
positive variable Zboil(N,j,comp,Net);
equation VaporZb(N,j,comp,Net);
VaporZb(N,j,comp,Net).. (Zboil(N,j,comp,Net))**3-(Zboil(N,j,comp,Net))**2+(Zboil(N,j,comp,Net))
                        *((aiEOSb(comp)*P(N,j,Net)/((0.00008314*Tb(comp))**2))
                        -(biEOS(comp)*P(N,j,Net)/(0.00008314*Tb(comp)))
                        -(biEOS(comp)*P(N,j,Net)/(0.00008314*Tb(comp)))**2)
                        -((aiEOSb(comp)*P(N,j,Net)/((0.00008314*Tb(comp))**2)))
                        *(biEOS(comp)*P(N,j,Net)/(0.00008314*Tb(comp))) =e= 0;

*Entalpía de vaporización (kJ/mol)
parameter DHvap(comp), Hvib(comp);
DHVap(comp)=( C1v(comp)*( (1-Tred(comp))**( C2v(comp) + (C3v(comp)*Tred(comp)) + (C4v(comp)*(Tred(comp)**2)) + (C5v(comp)*(Tred(comp)**3)) ) ) );
HVib(comp)=( (C1c(comp)*(Tb(comp)-Tref)) + ((C2c(comp)/2)*((Tb(comp)**2)-(Tref**2))) + ((C3c(comp)/3)*((Tb(comp)**3)-(Tref**3))) + ((C4c(comp)/4)*((Tb(comp)**4)-(Tref**4))) + ((C5c(comp)/5)*((Tb(comp)**5)-(Tref**5))) + ((C6c(comp)/6)*((Tb(comp)**6)-(Tref**6))) + Hform(comp));
variable depHvib(N,j,comp,Net);
equation EqdepHvib(N,j,comp,Net);
EqdepHvib(N,j,comp,Net).. depHvib(N,j,comp,Net) =e= (8.314/1000)*Tb(comp)*(Zboil(N,j,comp,Net)-1)
                                                +(1+mEOS(comp))*((aiEOSb(comp)**0.5)/biEOS(comp))
                                                *log(Zboil(N,j,comp,Net)/(Zboil(N,j,comp,Net)+(biEOS(comp)*P(N,j,Net)/(0.00008314*Tb(comp)))));

*Constantes de Cp (kJ/mol.K) de líquido
parameter
C1l(comp)
/
iButene         0.08768
Ethanol         0.10264
nButene         0.18205
ETBE            0.11096
/
C2l(comp)
/iButene        0.0002171
Ethanol         -0.00013963
nButene         -0.001611
ETBE            0.00031422
/
C3l(comp)
/iButene        -9.15300E-07
Ethanol         -3.03410E-08
nButene         1.19630E-05
ETBE            1.74800E-07
/
C4l(comp)
/iButene        2.2660E-09
Ethanol         2.0386E-09
nButene         -3.7454E-08
ETBE            0
/
C5l(comp)
/iButene        0
Ethanol         0
nButene         4.5027E-11
ETBE            0
/
;

*Entalpía de la fase liquida (kJ/mol)
variable HLi(N,j,comp,Net),HL(N,j,Net);
equation EqHLi(N,j,comp,Net),EqHL(N,j,Net);
EqHLi(N,j,comp,Net).. HLi(N,j,comp,Net)=e=HVib(comp)-DHVap(comp)
                        +((C1l(comp)*(Temp(N,j,Net)-Tb(comp))) + ((C2l(comp)/2)*((Temp(N,j,Net)**2)-(Tb(comp)**2)))
                        +((C3l(comp)/3)*((Temp(N,j,Net)**3)-(Tb(comp)**3))) + ((C4l(comp)/4)*((Temp(N,j,Net)**4)-(Tb(comp)**4)))
                        +((C5l(comp)/5)*((Temp(N,j,Net)**5)-(Tb(comp)**5))))+depHvib(N,j,comp,Net);

EqHL(N,j,Net).. HL(N,j,Net)=e=sum(comp,HLi(N,j,comp,Net)*x(N,j,comp,Net)/100);
*-------------------------------------------------------------------------------
*                                Sección 12
*                  Cálculo de entalpía de alimentación
*-------------------------------------------------------------------------------
*Entalpía de la alimentación de butenos
parameter HV_b(comp)    "Entalpía de vapor de la alimentación [kJ/mol]"
          Tred_b(comp)  "Temperatura reducida alimentación [-]"
          DHVap_b(comp) "Entalpía de vaporización alimentación [kJ/mol]"
          HL_b(comp)    "Entalpía de líquido de la alimentación [kJ/mol]";
HV_b(comp)=( (C1c(comp)*(TaliB-Tref)) + ((C2c(comp)/2)*((TaliB**2)-(Tref**2))) + ((C3c(comp)/3)*((TaliB**3)-(Tref**3))) + ((C4c(comp)/4)*((TaliB**4)-(Tref**4))) + ((C5c(comp)/5)*((TaliB**5)-(Tref**5))) + ((C6c(comp)/6)*((TaliB**6)-(Tref**6))) + Hform(comp));
Tred_b(comp)=TaliB/Tcrit(comp);
DHVap_b(comp)=( C1v(comp)*( (1-Tred_b(comp))**( C2v(comp) + (C3v(comp)*Tred_b(comp)) + (C4v(comp)*(Tred_b(comp)**2)) + (C5v(comp)*(Tred_b(comp)**3)) ) ) );
HL_b(comp)=HV_b(comp)-DHVap_b(comp);
parameter alphaEOSbut(comp), aiEOSbut(comp), aEOSbut(N,j), bEOSbut(N,j);
alphaEOSbut(comp)=(1+mEOS(comp)*(1-(TaliB/Tcrit(comp))**(1/2)))**2;
aiEOSbut(comp)=alphaEOSbut(comp)*0.42747*((0.00008314*TcritSRK(comp))**2)/Pcrit(comp);
bEOSbut(N,j)=sum(comp,(zb(N,j,comp)/100)*biEOS(comp));
aEOSbut(N,j)=sum(comp,sum(comp1, (zb(N,j,comp)/100)*(zb(N,j,comp1)/100)*(aiEOSbut(comp)*aiEOSbut(comp1))**0.5));


*Zbut se calcula para todas las etapas internas
positive variable Zbut(N,j,Net);
equation VaporZbut(N,j,Net);
VaporZbut(N,j,Net)$((ord(Net) ne card(Net)) and (ord(Net) ne 1)).. (Zbut(N,j,Net))**3-(Zbut(N,j,Net))**2+(Zbut(N,j,Net))
                        *((aEOSbut(N,j)*P(N,j,Net)/((0.00008314*TaliB)**2))
                        -(bEOSbut(N,j)*P(N,j,Net)/(0.00008314*TaliB))
                        -(bEOSbut(N,j)*P(N,j,Net)/(0.00008314*TaliB))**2)
                        -((aEOSbut(N,j)*P(N,j,Net)/((0.00008314*TaliB)**2)))
                        *(bEOSbut(N,j)*P(N,j,Net)/(0.00008314*TaliB)) =e= 0;

*HFB se calcula para todas las etapas internas
variable HFB(N,j,Net) "Entalpía de la alimentación de butenos";
equation EqHFB(N,j,Net);
EqHFB(N,j,Net)$((ord(Net) ne card(Net)) and (ord(Net) ne 1)).. HFB(N,j,Net) =e= sum(comp,(zb(N,j,comp)/100)*(HL_b(comp)+(8.314/1000)*TaliB*(Zbut(N,j,Net)-1)
                        +(1+mEOS(comp))*((aEOSbut(N,j)**0.5)/bEOSbut(N,j))
                        *log(Zbut(N,j,Net)/(Zbut(N,j,Net)+(bEOSbut(N,j)*P(N,j,Net)/(0.00008314*TaliB))))));

*Entalpía de la alimentación de etanol
parameter HV_e(comp)    "Entalpía de vapor de la alimentación [kJ/mol]"
          Tred_e(comp)  "Temperatura reducida alimentación [K]"
          DHVap_e(comp) "Entalpía de vaporización alimentación [kJ/mol]"
          HL_e(comp)    "Entalpía de líquido de la alimentación [kJ/mol]";
HV_e(comp)=( (C1c(comp)*(TaliE-Tref)) + ((C2c(comp)/2)*((TaliE**2)-(Tref**2))) + ((C3c(comp)/3)*((TaliE**3)-(Tref**3))) + ((C4c(comp)/4)*((TaliE**4)-(Tref**4))) + ((C5c(comp)/5)*((TaliE**5)-(Tref**5))) + ((C6c(comp)/6)*((TaliE**6)-(Tref**6))) + Hform(comp));
Tred_e(comp)=TaliE/Tcrit(comp);
DHVap_e(comp)=( C1v(comp)*( (1-Tred_e(comp))**( C2v(comp) + (C3v(comp)*Tred_e(comp)) + (C4v(comp)*(Tred_e(comp)**2)) + (C5v(comp)*(Tred_e(comp)**3)) ) ) );
HL_e(comp)=HV_e(comp)-DHVap_e(comp);
parameter alphaEOSeth(comp), aiEOSeth(comp), aEOSeth, bEOSeth;
alphaEOSeth(comp)=(1+mEOS(comp)*(1-(TaliE/Tcrit(comp))**(1/2)))**2;
aiEOSeth(comp)=alphaEOSeth(comp)*0.42747*((0.00008314*TcritSRK(comp))**2)/Pcrit(comp);
bEOSeth=sum(comp,(ze(comp)/100)*biEOS(comp));
aEOSeth=sum(comp,sum(comp1, (ze(comp)/100)*(ze(comp1)/100)*(aiEOSeth(comp)*aiEOSeth(comp1))**0.5));

*Zeth se calcula para todas las etapas internas
positive variable Zeth(N,j,Net);
equation VaporZeth(N,j,Net);
VaporZeth(N,j,Net)$((ord(Net) ne card(Net)) and (ord(Net) ne 1)).. (Zeth(N,j,Net))**3-(Zeth(N,j,Net))**2+(Zeth(N,j,Net))
                        *((aEOSeth*P(N,j,Net)/((0.00008314*TaliE)**2))
                        -(bEOSeth*P(N,j,Net)/(0.00008314*TaliE))
                        -(bEOSeth*P(N,j,Net)/(0.00008314*TaliE))**2)
                        -((aEOSeth*P(N,j,Net)/((0.00008314*TaliE)**2)))
                        *(bEOSeth*P(N,j,Net)/(0.00008314*TaliE)) =e= 0;

*HFE se alcula para todas las etapas internas
variable  HFE(N,j,Net)   "Entalpía de la alimentación de etanol [kJ/mol]";
equation EqHFE(N,j,Net);
EqHFE(N,j,Net)$((ord(Net) ne card(Net)) and (ord(Net) ne 1)).. HFE(N,j,Net) =e= sum(comp,(ze(comp)/100)*(HL_e(comp)+(8.314/1000)*TaliE*(Zeth(N,j,Net)-1)
                        +(1+mEOS(comp))*((aEOSeth**0.5)/bEOSeth)
                        *log(Zeth(N,j,Net)/(Zeth(N,j,Net)+(bEOSeth*P(N,j,Net)/(0.00008314*TaliE))))));
*-------------------------------------------------------------------------------
*                                Sección 13
*                Definicion de parametros binarios
*-------------------------------------------------------------------------------

*Parametro que determina si las etapas de rxn deben considerarse como etapas de equilibrio
scalar CASE "0 indica que en las etapas de rxn si hay equilibrio"/0/;

*Existencia de catalizador
parameter yc(Net) "1 indica que en la etapa si hay catalizador";

*Existencia de reflujo
parameter yr(Net) "1 indica que en la etapa si hay reflujo";
yr(Net)$((ord(Net) ne card(Net)) and (ord(Net) ne 1))=0;
yr('2')=1;
*Existencia de boil up
parameter yb(Net) "1 indica que en la etapa si hay boil up";

*Permite saber si la etapa es real o sobrante
parameter par(Net) "1 indica que la etapa es real fisicamente";

*por definicion el rehervidor y el condensador existen en la columna
par('1')=1;
par(Net)$(ord(Net) eq card(Net))=1;

*Existencia de relaciones de equilibrio
parameter ye(Net) "1 indica que la etapa es de equilibrio";

*Existencia de alimentacion
parameter yf(Net,F) "1 indica que en la etapa hay alimentacion de F";

*-------------------------------------------------------------------------------
*                                Sección 14
*                           Restricciones logicas
*-------------------------------------------------------------------------------
*Esta seccion no se requiere en este caso
$ontext
scalar cmej /1/;
scalar NCmax "numero maximo de etapas reactivas" /15/;

equation logic1(Net) "The boil up stage is below the reflux stage";
logic1(Net)$(ord(Net)>1 and ord(Net)<card(Net))..cmej*(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yr(Net1)))=g=cmej*(yb(Net));

equation logic2 "There is one reflux stage";
logic2..cmej*(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le (card(Net1)-1))),yr(Net1)))=e=cmej*1;

equation logic3 "There is one boil up stage";
logic3..cmej*(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le (card(Net1)-1))),yb(Net1)))=e=cmej*1;

equation logic4(F)"There is one feed stage of EtOH and there is one feed stage of butenes";
logic4(F)..cmej*(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le (card(Net1)-1))),yf(Net1,F)))=e=cmej*1;

equation logic6 "There is a maximum number of  catalytic stages";
logic6.. cmej*(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le (card(Net1)-1))),yc(Net1))) =e=cmej*NCmax;

equation logic7(Net,F) "Both feed stages are below the reflux";
logic7(Net,F)$(ord(Net)>1 and ord(Net)<card(Net))..cmej*(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yr(Net1)))=g=cmej*yf(Net,F);

equation logic8(Net,F) "The boil up stage is below the feed stages";
logic8(Net,F)$(ord(Net)>1 and ord(Net)<card(Net))..cmej*(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yf(Net1,F)))=g=cmej*yb(Net);

equation logic9(Net)  "The EtOH feed is above the butenes feed";
logic9(Net)$(ord(Net)>1 and ord(Net)<card(Net))..cmej*(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yf(Net1,'1')))=g=cmej*yf(Net,'2');

equation logic10(Net) "The catalytic stages are below the EtOH feed stage";
logic10(Net)$(ord(Net)>1 and ord(Net)<card(Net))..cmej*(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yf(Net1,'1')))=g=cmej*yc(Net);

equation logic11(Net) "The catalytic stages are above the butenes  feed stage";
logic11(Net)$(ord(Net)>1 and ord(Net)<card(Net))..cmej*((sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yf(Net1,'2')))-(yf(Net,'2')))=l=cmej*(1-yc(Net));

equation logic12(Net) "The catalytic stages are below the reflux stage";
logic12(Net)$(ord(Net)>1 and ord(Net)<card(Net))..cmej*(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yr(Net1)))=g=cmej*yc(Net);

equation logic13(Net) "The catalytic stages are above the boil up stage";
logic13(Net)$(ord(Net)>1 and ord(Net)<card(Net))..cmej*((sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yb(Net1)))-(yb(Net)))=l=cmej*(1-yc(Net));
$offtext

*-------------------------------------------------------------------------------
*                                Sección 15
*                         Ecuaciones del condensador
*-------------------------------------------------------------------------------

*Condiciones iniciales (operación en estado estable)
equation BalMasaC0,BalMasaParcialC0(comp),SumaC0,EquilibrioC0(comp),BalEnergiaC0;
BalMasaC0.. 0=e=V('1','1','2')-V('1','1','1')*(1+RR('1','1'));
BalMasaParcialC0(comp).. 0=e=V('1','1','2')*y('1','1',comp,'2')-V('1','1','1')*x('1','1',comp,'1')*(1+RR('1','1'));
SumaC0.. sum(comp,y('1','1',comp,'1')-x('1','1',comp,'1'))=e=0;
EquilibrioC0(comp).. y('1','1',comp,'1')*P('1','1','1')*phi('1','1',comp,'1')=e=Psat('1','1',comp,'1')*gamma('1','1',comp,'1')*x('1','1',comp,'1');
BalEnergiaC0.. 0=e=V('1','1','2')*HV('1','1','2')-V('1','1','1')*(1+RR('1','1'))*HL('1','1','1')-QC('1','1');

*Flujo de liquido fijo
equation fixedL(N,j);
fixedL(N,j)..L(N,j,'1')=e=0;
*-------------------------------------------------------------------------------
*                                Sección 16
*                    Ecuaciones de la columna (Punto Inicial)
*-------------------------------------------------------------------------------

*Condiciones iniciales (operación en estado estable)
equations BalMasa0(Net,Net1),BalMasaParcial0(comp,Net,Net1),Suma0(Net),BalEnergia0(Net,Net1);
BalMasa0(Net,Net1)$((ord(Net)>1 and ord(Net)<card(Net)) and (ord(Net1) eq card(Net1)))..0=e=yf(Net,'1')*FE+yf(Net,'2')*FB+RR('1','1')*V('1','1','1')*yr(Net)+BR('1','1')*L('1','1',Net1)*yb(Net)+L('1','1',Net-1)+V('1','1',Net+1)-L('1','1',Net)-V('1','1',Net)+yc(Net)*(sum(comp,Nu(comp))*mcat*Rx('1','1',Net)) ;
BalMasaParcial0(comp,Net,Net1)$((ord(Net)>1 and ord(Net)<card(Net)) and (ord(Net1) eq card(Net1)))..0=e=yf(Net,'1')*FE*ze(comp)+yf(Net,'2')*FB*zb('1','1',comp)+RR('1','1')*V('1','1','1')*yr(Net)*x('1','1',comp,'1')+BR('1','1')*L('1','1',Net1)*yb(Net)*y('1','1',comp,Net1)+L('1','1',Net-1)*x('1','1',comp,Net-1)+V('1','1',Net+1)*y('1','1',comp,Net+1)-L('1','1',Net)*x('1','1',comp,Net)-V('1','1',Net)*y('1','1',comp,Net)+100*yc(Net)*(Nu(comp)*mcat*Rx('1','1',Net));
Suma0(Net)$(ord(Net)>1 and ord(Net)<card(Net)).. sum(comp,x('1','1',comp,Net)-y('1','1',comp,Net))=e=0;
BalEnergia0(Net,Net1)$(ord(Net)>1 and ord(Net)<card(Net) and ord(Net1) eq card(Net1))..0=e=yf(Net,'1')*FE*HFE('1','1',Net)+yf(Net,'2')*FB*HFB('1','1',Net)+RR('1','1')*V('1','1','1')*yr(Net)*HL('1','1','1')+BR('1','1')*L('1','1',Net1)*yb(Net)*HV('1','1',Net1)+L('1','1',Net-1)*HL('1','1',Net-1)+V('1','1',Net+1)*HV('1','1',Net+1)-L('1','1',Net)*HL('1','1',Net)-V('1','1',Net)*HV('1','1',Net);

*Relaciones de equilibrio
equations Equilibrio10(comp,Net);
Equilibrio10(comp,Net)$(ord(Net)>1 and ord(Net)<card(Net))..0=e=ye(Net)*((y('1','1',comp,Net)*P('1','1',Net)*phi('1','1',comp,Net))-(Psat('1','1',comp,Net)*gamma('1','1',comp,Net)*x('1','1',comp,Net)));
equations Equilibrio20(comp,Net);
Equilibrio20(comp,Net)$(ord(Net)>1 and ord(Net)<card(Net) and ord(comp) ne 1)..0=e=(sum(Net1$(ord(Net1) ge 2 and ord(Net1) le ord(Net)),yr(Net1)))*(1-ye(Net))*(y('1','1',comp,Net)-y('1','1',comp,Net+1));
equation Equilibrio30(comp,Net);
Equilibrio30(comp,Net)$(ord(Net)>1 and ord(Net)<card(Net) and ord(comp) ne 1)..0=e=(1-sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yr(Net1)))*(1-ye(Net))*(x('1','1',comp,Net)-x('1','1',comp,Net-1));
Equation Equilibrio40(Net,Net1);
Equilibrio40(Net,Net1)$(ord(Net)>1 and ord(Net)<card(Net) and ord(Net1) eq card(Net1))..0=e=(1-ye(Net))*(V('1','1',Net)-V('1','1',Net+1)-BR('1','1')*L('1','1',Net1)*yb(Net));
*-------------------------------------------------------------------------------
*                                Sección 17
*                        Ecuaciones del rehervidor
*-------------------------------------------------------------------------------

*Condiciones iniciales (operación en estado estable)
equation BalMasaR0(Net),BalMasaParcialR0(comp,Net),SumaR0(Net),EquilibrioR0(comp,Net),BalEnergiaR0(Net);
BalMasaR0(Net)$(ord(Net) eq card(Net))..0=e=L('1','1',Net-1)-L('1','1',Net)*(1+BR('1','1'));
BalMasaParcialR0(comp,Net)$(ord(Net) eq card(Net))..0=e=L('1','1',Net-1)*x('1','1',comp,Net-1)-L('1','1',Net)*(x('1','1',comp,Net)+BR('1','1')*y('1','1',comp,Net));
SumaR0(Net)$(ord(Net) eq card(Net)).. sum(comp,y('1','1',comp,Net)-x('1','1',comp,Net))=e=0;
EquilibrioR0(comp,Net)$(ord(Net) eq card(Net))..y('1','1',comp,Net)*P('1','1',Net)*phi('1','1',comp,Net)=e=Psat('1','1',comp,Net)*gamma('1','1',comp,Net)*x('1','1',comp,Net);
BalEnergiaR0(Net)$(ord(Net) eq card(Net))..0=e=QR('1','1')+L('1','1',Net-1)*HL('1','1',Net-1)-L('1','1',Net)*HL('1','1',Net)-BR('1','1')*L('1','1',Net)*HV('1','1',Net);

*Variable fija del flujo de vapor en la ultima etapa
equation fixedV(N,j,Net);
fixedV(N,j,Net)$(ord(Net) eq card(Net))..V(N,j,Net)=e=0;
*-------------------------------------------------------------------------------
*                                Sección 18
*               Relaciones hidráulicas para todas las etapas internas
*-------------------------------------------------------------------------------
*Caracteristicas del catalizador
scalar fracvol /0.3/;
scalar fracEnvelop /0.5/;

*Definición de velocidad de vapor
positive variables far(N,j,Net) "Factor de areación [-]";
equations Eqfa(N,j,Net);
Eqfa(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net)).. par(Net)*(far(N,j,Net))=e=par(Net)*(0.981*exp(-0.411*((V(N,j,Net)/(rhoV(N,j,Net))/hora)*(rhoV(N,j,Net)*sum(comp,MW(comp)*y(N,j,comp,Net)/100)/1000)**(0.5))/At));

positive variable hD(N,j,Net)   "Altura del líquido por encima del divisor [m]";
equations EqhD(N,j,Net);
EqhD(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net)).. (hD(N,j,Net))=e=(0.6*(((((L(N,j,Net)/sum(comp,rho(N,j,comp,Net)*x(N,j,comp,Net)/100))/hora)/Lw))**(2/3)));

positive variable uhv(N,j,Net) "Velocidad del vapor por los agujeros [m/s]";
equations Equhv(N,j,Net);
Equhv(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net)).. par(Net)*(uhv(N,j,Net))=e=par(Net)*((V(N,j,Net)/(rhoV(N,j,Net))/hora)/A0);

positive variable unv(N,j,Net) "Velocidad del vapor por el plato [m/s]";
equations Equnv(N,j,Net);
Equnv(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net)).. par(Net)*unv(N,j,Net)=e=par(Net)*((V(N,j,Net)/(rhoV(N,j,Net))/hora)/At);

*Definicion de velocidad del liquido
positive variable ul(N,j,Net) "Velocidad del líquido en el derramadero [m/s]";
equations Equl(N,j,Net);
Equl(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net)).. par(Net)*ul(N,j,Net)=e=par(Net)*((L(N,j,Net)/(sum(comp,rho(N,j,comp,Net)*x(N,j,comp,Net)/100))/hora)/Ad);

*Carga de liquido
positive variable hcl(N,j,Net)  "Altura del líquido libre en régimen de spray [m]"
equation Eqhcl(N,j,Net);
scalar consmach /1e-20/;
Eqhcl(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net)).. par(Net)*hcl(N,j,Net)=e=par(Net)*((0.157*(poro**(-0.791))/(1+1.04E-4*(((((L(N,j,Net)+consmach)/sum(comp,rho(N,j,comp,Net)*x(N,j,comp,Net)/100))/hora)/Lw)**(-0.59))
                                                        *(poro**(-1.791))))*(da**0.833)
                                                        *(996/(sum(comp,rho(N,j,comp,Net)*x(N,j,comp,Net)/100)*sum(comp,MW(comp)*x(N,j,comp,Net)/100)/1000))**(0.5*(1-0.91*da/poro)));
positive variable Csbf(N,j,Net);
equation EqCsbf(N,j,Net);
EqCsbf(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net)).. par(Net)*(Csbf(N,j,Net))=e=par(Net)*(0.37*(((sqr(da)*sigma(N,j,Net)/(sum(comp,rho(N,j,comp,Net)*x(N,j,comp,Net)/100)*sum(comp,MW(comp)*x(N,j,comp,Net)/100)/1000)))**0.125)
                                                        *((((rhoV(N,j,Net))*sum(comp,MW(comp)*y(N,j,comp,Net)/100)/1000)/(sum(comp,rho(N,j,comp,Net)*x(N,j,comp,Net)/100)*sum(comp,MW(comp)*x(N,j,comp,Net)/100)/1000))**0.1)
                                                        *((HS/hcl(N,j,Net))**0.5));
*Carga de liquido en etapas cataliticas
positive variable Lload(N,j,Net) "carga de liquido en etapas cataliticas [m_s]";
equation eqLload(N,j,Net);
eqLload(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net))..(1-fracvol)*((3.1415926/4)*(D**2))*Lload(N,j,Net)=e=ul(N,j,net)*Ad ;

*Factor de flujo de vapor en etapas cataliticas
positive variable Ffactor(N,j,Net) "Factor de flujo de vapor para etapas cataliticas [Pa**0.5]"
equation eqFfactor(N,j,Net);
eqFfactor(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net)).. (1-fracvol)*(3.1415926/4)*(D**2)*(((rhov(N,j,net))*(sum(comp,(y(N,j,comp,Net)/100)*MW(comp)))*(1/1000))**(1/2))*(Ffactor(N,j,Net))=e=(V(N,j,net)*(1/60))*(sum(comp,(y(N,j,comp,Net)/100)*MW(comp)*(1/1000)));

*Caida de presion
positive variables DPL(N,j,Net) "Caída de presión por la presencia de líquido [bar]";
equations EqDPL(N,j,Net);
EqDPL(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net)).. (DPL(N,j,Net))=e=((far(N,j,Net)*9.81*(sum(comp,rho(N,j,comp,Net)*x(N,j,comp,Net)/100)*sum(comp,MW(comp)*x(N,j,comp,Net)/100)/1000)*(hD(N,j,Net)+hw))/100000);

positive variables DPS(N,j,Net) "Caída de presión debido a la presencia de los agujeros - seco [bar]";
equations EqDPS(N,j,Net);
EqDPS(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net)).. (DPS(N,j,Net))=e=((1/(2*sqr(K0)))*( (((sqr(V(N,j,Net)/(rhoV(N,j,Net))/hora)/A0)) )*((rhoV(N,j,Net))*sum(comp,MW(comp)*y(N,j,comp,Net)/100)/1000)*(1-sqr(poro)))/100000);

positive variable DPq(N,j,Net)      "Caída de presión en el derramadero [bar]";
equations EqDPq(N,j,Net);
EqDPq(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net)).. DPq(N,j,Net)=e=(1/(100000))*1.62*((sum(comp,rho(N,j,comp,Net)*x(N,j,comp,Net)/100)*sum(comp,MW(comp)*x(N,j,comp,Net)/100)/1000))/(sqr(Lw*hw))*(sqr((L(N,j,Net)/sum(comp,rho(N,j,comp,Net)*x(N,j,comp,Net)/100))/hora)+sqr((V(N,j,Net)/(rhoV(N,j,Net))/hora)));

positive variables DP(N,j,Net)  "Caída de presión total [bar]";
positive variable dPcat(N,j,net)    "caida de presion por catalizador en etapas cataliticas [bar]";
equations EqDP(N,j,Net),EqDPR(N,j,Net),EqdPcat(N,j,net),EqP(N,j,Net),EqPC(N,j,Net),EqPR(N,j,Net) "Definición de presión por etapa [bar]";

EqDPR(N,j,Net)$(ord(Net) eq card(Net)).. DP(N,j,Net)=e=DP(N,j,Net-1);
EqdPcat(N,j,net)$(ord(Net)>1 and ord(Net)<card(Net))..dPcat(N,j,net)=e=hs*fracEnvelop*(0.001)*(   (5.69228924748553E-06)*((Lload(N,j,Net)*60*60)**3.05308055949085)*((Ffactor(N,j,Net))**7.851695947) + 1.367015225*((Ffactor(N,j,Net))**1.764157687)    );
EqDP(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net)).. (DP(N,j,Net))=e=(DPS(N,j,Net)+DPL(N,j,Net))+yc(Net)*dPcat(N,j,net);
EqP(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net)).. P(N,j,Net)=e=P(N,j,Net-1)+par(Net)*DP(N,j,Net);
EqPC(N,j,Net)$(ord(Net) eq 1).. P(N,j,Net)=e=Pop;
EqPR(N,j,Net)$(ord(Net) eq card(Net)).. P(N,j,Net)=e=P(N,j,Net-1);

*Efectos indeseados en la columna
*Downflow flooding (inundación en los derramaderos)
equation DownFlood(N,j,Net);
DownFlood(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net))..0=g=((HD(N,j,Net)+((DP(N,j,Net)+DPq(N,j,Net))*100000)
                                                        /(9.81*(((sum(comp,rho(N,j,comp,Net)*x(N,j,comp,Net)/100)*sum(comp,MW(comp)*x(N,j,comp,Net)/100)/1000))
                                                        -(rhoV(N,j,Net)*sum(comp,MW(comp)*y(N,j,comp,Net)/100)/1000))))-(HS))*par(Net);
*Entrainment flooding (inundación por arrastre de líquido)
equation EntrainFloodV(N,j,Net), EntrainFloodL(N,j,Net);
EntrainFloodV(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net))..par(Net)*((unv(N,j,Net))-(Csbf(N,j,Net)*(((((sum(comp,rho(N,j,comp,Net)*x(N,j,comp,Net)/100)*sum(comp,MW(comp)*x(N,j,comp,Net)/100)/1000))
                                                        -(rhoV(N,j,Net)*sum(comp,MW(comp)*y(N,j,comp,Net)/100)/1000)))
                                                        /(rhoV(N,j,Net)*sum(comp,MW(comp)*y(N,j,comp,Net)/100)/1000))**0.5))=l=0;
EntrainFloodL(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net))..par(Net)*((ul(N,j,Net))-((sigma(N,j,Net)*9.81*(((sum(comp,rho(N,j,comp,Net)*x(N,j,comp,Net)/100)*sum(comp,MW(comp)*x(N,j,comp,Net)/100)/1000))
                                                        -(rhoV(N,j,Net)*sum(comp,MW(comp)*y(N,j,comp,Net)/100)/1000))
                                                        /((sum(comp,rho(N,j,comp,Net)*x(N,j,comp,Net)/100)*sum(comp,MW(comp)*x(N,j,comp,Net)/100)/1000)**2))**(1/4)))=l=0
*Weeping (lloriqueo)
equation Weep(N,j,Net);
Weep(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net)).. 0=g=(((0.68-0.12)/(((rhoV(N,j,Net)*sum(comp,MW(comp)*y(N,j,comp,Net)/100)/1000)
                                                                /((sum(comp,rho(N,j,comp,Net)*x(N,j,comp,Net)/100)*sum(comp,MW(comp)*x(N,j,comp,Net)/100)/1000)
                                                                *9.81*far(N,j,Net)*(hw+hd(N,j,Net))))**0.5))-(uhv(N,j,Net)))*par(Net);
*Catalyst flooding (inundación del empaque del catalizador)
equation catflood(N,j,net);
catflood(N,j,net)..yc(net)*(dPcat(N,j,net)-(12E-3)*hs*fracEnvelop)=l=0;

*Construcción de la columna
equation Size "Tamaño del equipo";
Size.. 1*Htotal =e= 1*((1+Sfactor)*sum(Net$(ord(Net)>1 and ord(Net)<card(Net)),HS*par(Net)));

equation  ammountcat "Espacio disponible para el catalizador";
ammountcat..mcat=l=(fracvol)*((3.1415926/4)*(D**2))*(hs*fracEnvelop)*770;

equation DtoLratio "Relacion entre el diametro y la altura";
DtoLratio..htotal/D=l=20;
*-------------------------------------------------------------------------------
*                                Sección 19
*                            Función objetivo
*-------------------------------------------------------------------------------
parameters alfa1 "Peso para pureza de ETBE" /1e4/
           alfa2 "Peso para carga térmica del rehervidor" /500/
           alfa3 "Peso para relación de reflujo" /100/
           CostQr "Costo de la carga térmica del rehervidor [$/yr]"
           CostQc "Costo de la carga térmica del condensador [$/yr]"
           CostB  "Ganacia del ETBE en fondos [$/yr]"
           CostEth"Costo de alimentación de etanol [$/yr]"
           CostBut"Costo de alimentación de butanos [$/yr]"
           year   "Operational hours per year [hr/yr]" /8000/
           CostCat"Costo del catalizador [$/kg]" /7.7/
           AF     "Factor de anualización (5 años, 5% de tasa de interés) [1/yr]"
           MS     "Marshall & Swift coefficient" /1050/
           FM     "Material factor (carbon steel)" /1/
           FPres  "Pressure factor (up to 200psi=13.78bar)" /1.15/
           C0     "Costo de inversión inicial AF(Cr1+Cc1) [$]" /10000/
           CT     "Costo de inversión para etapas [$]"
           Csh    "Costo de inversión para coraza [$]"
           Fcol   "Factor de costo de la columna [-]";
Fcol=FM*FPres;
AF=(0.05/(1-1/(1+0.05)**5));
CT=AF*(MS/280)*4.7*Fcol;
Csh=AF*(MS/280)*101.9*(2.18+Fcol);
CostQr=146.8/(hora);
CostQc=24.5/(hora);
CostB=25.3*3600*year/(1000*hora);
CostEth=15*3600*year/(1000*hora);
CostBut=8.25*3600*year/(1000*hora);

variables zobj;

equation Fobj(Net);


Fobj(Net)$(ord(Net) eq card(Net)).. zobj=e=1*((((CostEth*FE+CostBut*FB+(CostQr*Qr('1','1'))+(CostQc*Qc('1','1')))/year)*year)
+(C0+AF*(mcat*CostCat*sum(Net1,yc(Net1))+CT*((D/0.3048)**1.55)*(sum(Net1,HS*par(Net1))/0.3048)+Csh*((D/0.3048)**1.066)*((Htotal/0.3048)**0.802)))
-((((CostB*L('1','1',Net)))/year)*year));
*-------------------------------------------------------------------------------
*                                Sección 20
*                            Cotas en las variables
*-------------------------------------------------------------------------------
*bounds
D.up=10;
D.lo=0.5;

hw.up=1;
hw.lo=0.0001;

hs.up=100;
hs.lo=0.1;

htotal.up=53;

at.lo=0.01;
at.up=100;

ad.lo=0.0001;
ad.up=20;

lw.up=50;

A0.lo=1e-12;
A0.up=100;

RR.lo(N,j)=1;
RR.up(N,j)=100;

Qr.up(N,j)=4000000000;
Qr.lo(N,j)=1000;

Qc.lo(N,j)=1000;
Qc.up(N,j)=900000000;

L.up(N,j,Net)=1000000;
V.up(N,j,Net)=1000000;

BR.up(N,j)=100;
P.lo(N,j,Net)=Pop;

zboil.up(N,j,comp,Net)=1.5;
zboil.lo(N,j,comp,Net)=0.7;

Zbut.up(N,j,Net)= 1.5;
Zbut.lo(N,j,Net)= 0.5;

Zeth.up(N,j,Net)= 1.5;
Zeth.lo(N,j,Net)= 0.5;

Psat.up(N,j,comp,Net)=100;

tao_nrtl.up(N,j,comp,comp1,Net)    =5;
tao_nrtl.lo(N,j,comp,comp1,Net)    =-5;

g_nrtl.up(N,j,comp,comp1,Net)=2;
g_nrtl.lo(N,j,comp,comp1,Net)=0;

gamma.up(N,j,comp,Net)=50;
gamma.lo(N,j,comp,Net)=0;

Ketbe.lo(N,j,Net)=0;
Ketbe.up(N,j,Net)=100;

Krate.up(N,j,Net)=1000000000000;
Krate.lo(N,j,Net)=-10000000;
Krate.scale(N,j,Net)=10000000;

Ka.up(N,j,Net)=100;

Rx.up(N,j,Net)=10000000;
Rx.lo(N,j,Net)=-1000000000;
Rx.scale(N,j,Net)=100000;

alphaEOS.up(N,j,comp,Net)=5;
aiEOS.up(N,j,comp,Net)=1e-3;
aiEOS.lo(N,j,comp,Net)=1e-6;
aiEOS.scale(N,j,comp,Net)=1e-4;

bEOS.up(N,j,Net)=0.01;

aEOS.up(N,j,Net)=1e-3;

phi.up(N,j,comp,Net)=2;

HVi.up(N,j,comp,Net)=1000;
HVi.lo(N,j,comp,Net)=-1000;

HV.up(N,j,Net)=1000;
HV.lo(N,j,Net)=-1000;

depHvib.lo(N,j,comp,Net)=-100;
depHvib.up(N,j,comp,Net)=100;

HLi.lo(N,j,comp,Net)=-1000;
HLi.up(N,j,comp,Net)=1000;

HL.lo(N,j,Net)=-1000;
HL.up(N,j,Net)=1000;

HFB.lo(N,j,Net)=-50;
HFB.up(N,j,Net)=1;

HFE.lo(N,j,Net)=-500;
HFE.up(N,j,Net)=0;

far.up(N,j,Net)=2;
far.lo(N,j,Net)=0.1;

hD.up(N,j,Net)=0.1;
hD.lo(N,j,Net)=0.0001;

DPL.up(N,j,Net)=0.1;

DPS.up(N,j,Net)=0.01;

DP.up(N,j,Net)=2;

DPq.up(N,j,Net)=0.5;

uhv.up(N,j,Net)=1000;
uhv.lo(N,j,Net)=0.4;

hcl.up(N,j,Net)=1;
hcl.lo(N,j,Net)=1e-6;


x.lo(N,j,comp,Net)= 0;
x.up(N,j,comp,Net)= 100;

y.lo(N,j,comp,Net)= 0;
y.up(N,j,comp,Net)= 100;

z.lo(N,j,Net)= 0.5 ;
z.up(N,j,Net)=  1.3  ;

Temp.lo(N,j,Net)= 200;
Temp.up(N,j,Net)= 441.89;

Tcritm.up(N,j,Net)=600;
Tcritm.lo(N,j,Net)=417.9;

bEOS.lo(N,j,Net)=min(biEOS('iButene'),biEOS('Ethanol'),biEOS('nButene'),biEOS('ETBE'));
aEOS.lo(N,j,Net)=1e-6;

unv.lo(N,j,Net)=0.01;
unv.up(N,j,Net)=89;

ul.up(N,j,Net)=3000;
ul.lo(N,j,Net)=0.001;

rho.up(N,j,comp,Net)=25000;
rho.lo(N,j,comp,Net)=8000;

rhoV.up(N,j,Net)=500;
rhoV.lo(N,j,Net)=Pop/(0.00008314*Temp.up(N,j,Net)*(Z.up(N,j,Net)));

Csbf.up(N,j,Net)=0.2;
Csbf.lo(N,j,Net)=0.1;

sigma.up(N,j,Net)=0.03;
sigma.lo(N,j,Net)=0.005;

*-------------------------------------------------------------------------------
*                                Sección 21
*                            Solución del modelo
*-------------------------------------------------------------------------------
model Rx_OCP_index1 /all/;
option nlp=conopt4;

Rx_OCP_index1.threads=0;
Rx_OCP_index1.SCALEOPT=1;
Rx_OCP_index1.workspace=4000;


*PROCESO ITERATIVO
set iter "Iterations" /i1*i100/;
set iner1 "Inner iterations 1" /i1/;
set iner2 "Inner iterations 2" /i1*i100/;
set extvar "External variables" /x1*x18/;
set scat(extvar) "external var: catalyst" /x4*x18/;
set sf1(extvar) "external var: feed 1" /x1/;
set sf2(extvar) "external var: feed 2" /x2/;
set sbup(extvar) "external var: boil up" /x3/;
set extineq "Convex polyhedron inequality equations" /g1*g19/;


set neigh /d1*d36/;
table directions(neigh,extvar)
$include neighbors.txt
;

scalar hstep "step" /1/;
parameter hiner1(iner1) "step: inner iterations 1";
loop(iner1,
if(ord(iner1) eq 1,
hiner1(iner1)=0;
else
hiner1(iner1)=hiner1(iner1-1)+(hstep/((card(iner1))-1));
);
);
parameter hiner2(iner2) "step: inner iterations 2";
loop(iner2,
if(ord(iner2) eq 1,
hiner2(iner2)=hstep;
else
hiner2(iner2)=hiner2(iner2-1)+hstep;
);
);
parameter extvarhstep(extvar,neigh) "Neighborhood considering hstep";
extvarhstep(extvar,neigh)=hstep*directions(neigh,extvar);

***Iterative process parameters and scalars
*stop1: stopping criterion in S4 for local optimality
scalar stop1 "Main stopping criterion";
*stop2: Stopping criterion that decides if the objective function value
*of a neighbor was propperly calculated (e.g, if the nlp solver finds an infeassible
*solution, the additional convergence procedure is executed).
scalar stop2 "Stoppong criterion 2";
*stop3: Stopping criterion for the line search
scalar stop3 "Stoppong criterion 3";
*stop4: Sttoping criterion that decides if the objective function value
*in the line search was propperly calculated.
scalar stop4 "Stoppong criterion 4";
scalar count;
parameter xvalue(iter,extvar) "value of x at each iteration";
parameter xvalueselect(iter,extvar) "selected value of x for the next iteration";
parameter gval(extineq) "value of inequality constraints at x";
parameter dvs(iter,neigh) "f(neigh)-f(x) at the iteration";
parameter dmin(iter) "Minimum value of dvs";
parameter fvalue(iter) "Objective function value";
parameter fplushvalue(iter,neigh) "Objective function to calculate dvs";
parameter fvalueiner(iter,iner2) "Objective function: inner iterations";
parameter selectd(iter,neigh) "ds: Selected direction for line search";
parameter mstatS2(iter)"model status: Feasibility of initialization point";
parameter mstatS4(iter,neigh,iner1) "model status: Feasibility of neighbors";
parameter mstatS5(iter,neigh,iner1) "model status: Feasibility of neighnors with convergence procedure" ;
parameter mstatS7(iter,iner2,iner1) "model status: Feasibility in line search";
parameter mstatS8(iter,iner2,iner1) "model status: Feasibility in line search with convergence procedure";
parameter objval "Objective function value";
parameter xvalinit(extvar) "Initialization of external variables";
parameter CPUtime;
scalar CPUtimeactual;
parameter xvalue_p(extvar) "partial value of the external variable";
parameter fex(extvar)"partial value of the f functions";
scalar cuent_gineq "count for inequality constraint";
parameter yaux(extvar,net) "auxiliar variable for reformualtion";

CPUtimeactual=0;
xvalue(iter,extvar)=0;
xvalueselect(iter,extvar)=0;
gval(extineq)=0;
dvs(iter,neigh)=0;
dmin(iter)=0;
fvalue(iter)=0;
fplushvalue(iter,neigh)=0;
fvalueiner(iter,iner2)=0;
selectd(iter,neigh)=0;
mstatS2(iter)=0;
mstatS4(iter,neigh,iner1)=0;
mstatS5(iter,neigh,iner1)=0;
mstatS7(iter,iner2,iner1)=0;
mstatS8(iter,iner2,iner1)=0;
xvalue_p(extvar)=0;
fex(extvar)=0;
yaux(extvar,net)=0;

***S1. Initialization of external variabkes
xvalinit("x1")=6;
xvalinit("x2")=0;
xvalinit("x3")=4;
xvalinit(scat)=1;
xvalinit("x4")=0;

xvalue('i1',extvar)=xvalinit(extvar);


stop1=0;
loop(iter$((stop1 eq 0) and (ord(iter) ne card(iter))),

***S2. Feasibility of external variables
         if(ord(iter) eq 1,
         xvalue_p(extvar)=xvalue(iter,extvar);

         fex(sf1)=xvalue_p("x1");
         fex(sf2)=xvalue_p("x1")+sum(scat,xvalue_p(scat))+xvalue_p("x2");
         fex(sbup)= xvalue_p("x1")+sum(scat,xvalue_p(scat))+xvalue_p("x2")+xvalue_p("x3");
         loop(scat,
         if(ord(scat) eq 1,
         fex(scat)=xvalue_p("x1")+xvalue_p("x4");
         else
         fex(scat)=xvalue_p(scat)+fex(scat-1);
         );
         );

         gval("g1")=2-xvalue_p("x1");
         gval("g2")=-xvalue_p("x3");
         gval("g3")=-xvalue_p("x4");
         gval("g4")=-xvalue_p("x2");
         gval("g5")=xvalue_p("x1")+sum(scat,xvalue_p(scat))+xvalue_p("x2")+xvalue_p("x3")-44;
         cuent_gineq=2;
         loop(extineq$(ord(extineq) ge 6),
         gval(extineq)=1+(sum(scat$(ord(scat) eq cuent_gineq-1),fex(scat)))-(sum(scat$(ord(scat) eq cuent_gineq),fex(scat)));
         cuent_gineq=cuent_gineq+1;
         );
                 if( sum(extineq$(gval(extineq) gt 0 ),gval(extineq)) gt 0,
                 stop1=1;
                 );

         );
***S1. Initialization if continuous variables
         if(ord(iter) eq 1,
         execute_loadpoint "al_arfaj2002_opt";
         else
         execute_loadpoint "DLR";
         );
***S2. Feasibilyty of continuous variables

yaux(extvar,net)=0;
yaux(extvar,net)$(ord(net) eq floor(fex(extvar)))=1-mod(fex(extvar),1);
yaux(extvar,net)$(ord(net) eq 1+floor(fex(extvar)))=mod(fex(extvar),1);

yf(Net,'1')=sum(sf1,yaux(sf1,net));
yf(Net,'2')=sum(sf2,yaux(sf2,net));
yb(Net)=sum(sbup,yaux(sbup,net));
yc(Net)=sum(scat,yaux(scat,net));
par(Net)$((ord(Net) ne card(Net)) and (ord(Net) ne 1))=(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yr(Net1)))+(yb(Net))-(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yb(Net1)));
ye(Net)$((ord(Net) ne card(Net)) and (ord(Net) ne 1))=par(Net)*(1-yc(Net))*CASE+par(Net)*(1-CASE);

solve Rx_OCP_index1 using nlp minimizing zobj;
mstatS2(iter)=Rx_OCP_index1.modelstat;
         if((mstatS2(iter) eq 3 or mstatS2(iter) eq 4 or mstatS2(iter) eq 5 or mstatS2(iter) eq 6 or mstatS2(iter) eq 11 or mstatS2(iter) eq 12 or mstatS2(iter) eq 13 or mstatS2(iter) eq 14 or mstatS2(iter) eq 18 or mstatS2(iter) eq 19) and ord(iter) eq 1,
         stop1=1;
         elseif (mstatS2(iter) eq 3 or mstatS2(iter) eq 4 or mstatS2(iter) eq 5 or mstatS2(iter) eq 6 or mstatS2(iter) eq 11 or mstatS2(iter) eq 12 or mstatS2(iter) eq 13 or mstatS2(iter) eq 14 or mstatS2(iter) eq 18 or mstatS2(iter) eq 19),
         fvalue(iter)=zobj.l;
         else
         execute_unload "DLR";
         fvalue(iter)=zobj.l;
                 if(ord(iter) eq 1,
                 execute_unload "ALG_INDUSTRIAL2";
                 );
         );


         loop(neigh,
***S3. k-Neighborhood

                  xvalue_p(extvar)=xvalue(iter,extvar)+extvarhstep(extvar,neigh);

                  fex(sf1)=xvalue_p("x1");
                  fex(sf2)=xvalue_p("x1")+sum(scat,xvalue_p(scat))+xvalue_p("x2");
                  fex(sbup)= xvalue_p("x1")+sum(scat,xvalue_p(scat))+xvalue_p("x2")+xvalue_p("x3");
                  loop(scat,
                  if(ord(scat) eq 1,
                  fex(scat)=xvalue_p("x1")+xvalue_p("x4");
                  else
                  fex(scat)=xvalue_p(scat)+fex(scat-1);
                  );
                  );

                  gval("g1")=2-xvalue_p("x1");
                  gval("g2")=-xvalue_p("x3");
                  gval("g3")=-xvalue_p("x4");
                  gval("g4")=-xvalue_p("x2");
                  gval("g5")=xvalue_p("x1")+sum(scat,xvalue_p(scat))+xvalue_p("x2")+xvalue_p("x3")-44;
                  cuent_gineq=2;
                  loop(extineq$(ord(extineq) ge 6),
                  gval(extineq)=1+(sum(scat$(ord(scat) eq cuent_gineq-1),fex(scat)))-(sum(scat$(ord(scat) eq cuent_gineq),fex(scat)));
                  cuent_gineq=cuent_gineq+1;
                  );

                          if( sum(extineq$(gval(extineq) gt 0 ),gval(extineq)) gt 0 ,
                          dvs(iter,neigh)=1;
                          else
                          dvs(iter,neigh)=-1;
                          );

                          if(dvs(iter,neigh) le 0,
                          stop2=0;

                                 loop(iner1$(stop2 eq 0),

                                         if(ord(iner1) eq 1,
***S4. Local optimality: Computation of the objective function value for the neighbors
                                         execute_loadpoint "DLR";
                                         yaux(extvar,net)=0;
                                         yaux(extvar,net)$(ord(net) eq floor(fex(extvar)))=1-mod(fex(extvar),1);
                                         yaux(extvar,net)$(ord(net) eq 1+floor(fex(extvar)))=mod(fex(extvar),1);

                                         yf(Net,'1')=sum(sf1,yaux(sf1,net));
                                         yf(Net,'2')=sum(sf2,yaux(sf2,net));
                                         yb(Net)=sum(sbup,yaux(sbup,net));
                                         yc(Net)=sum(scat,yaux(scat,net));
                                         par(Net)$((ord(Net) ne card(Net)) and (ord(Net) ne 1))=(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yr(Net1)))+(yb(Net))-(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yb(Net1)));
                                         ye(Net)$((ord(Net) ne card(Net)) and (ord(Net) ne 1))=par(Net)*(1-yc(Net))*CASE+par(Net)*(1-CASE);

                                         solve Rx_OCP_index1 using nlp minimizing zobj;
                                         mstatS4(iter,neigh,iner1)=Rx_OCP_index1.modelstat;

                                                 if(mstatS4(iter,neigh,iner1) eq 3 or mstatS4(iter,neigh,iner1) eq 4 or mstatS4(iter,neigh,iner1) eq 5 or mstatS4(iter,neigh,iner1) eq 6 or mstatS4(iter,neigh,iner1) eq 11 or mstatS4(iter,neigh,iner1) eq 12 or mstatS4(iter,neigh,iner1) eq 13 or mstatS4(iter,neigh,iner1) eq 14 or mstatS4(iter,neigh,iner1) eq 18 or mstatS4(iter,neigh,iner1) eq 19,
                                                 stop2=0;
                                                 dvs(iter,neigh)=1;
                                                 else
                                                 stop2=1;
                                                 fplushvalue(iter,neigh)=zobj.l;
                                                 dvs(iter,neigh)=(fplushvalue(iter,neigh)-fvalue(iter))/(hstep);
                                                 );

                                         else
***S4. Local optimality: Computation of the objective function value for the neighbors with the convergence procedure

                                                 if(ord(iner1) eq 2,
                                                 execute_loadpoint "DLR";
                                                 else
                                                 execute_loadpoint "DLR1";
                                                 );
                                         xvalue_p(extvar)=xvalue(iter,extvar)+((extvarhstep(extvar,neigh))/(hstep))*(hiner1(iner1));

                                         fex(sf1)=xvalue_p("x1");
                                         fex(sf2)=xvalue_p("x1")+sum(scat,xvalue_p(scat))+xvalue_p("x2");
                                         fex(sbup)= xvalue_p("x1")+sum(scat,xvalue_p(scat))+xvalue_p("x2")+xvalue_p("x3");
                                         loop(scat,
                                         if(ord(scat) eq 1,
                                         fex(scat)=xvalue_p("x1")+xvalue_p("x4");
                                         else
                                         fex(scat)=xvalue_p(scat)+fex(scat-1);
                                         );
                                         );



                                         yaux(extvar,net)=0;
                                         yaux(extvar,net)$(ord(net) eq floor(fex(extvar)))=1-mod(fex(extvar),1);
                                         yaux(extvar,net)$(ord(net) eq 1+floor(fex(extvar)))=mod(fex(extvar),1);

                                         yf(Net,'1')=sum(sf1,yaux(sf1,net));
                                         yf(Net,'2')=sum(sf2,yaux(sf2,net));
                                         yb(Net)=sum(sbup,yaux(sbup,net));
                                         yc(Net)=sum(scat,yaux(scat,net));
                                         par(Net)$((ord(Net) ne card(Net)) and (ord(Net) ne 1))=(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yr(Net1)))+(yb(Net))-(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yb(Net1)));
                                         ye(Net)$((ord(Net) ne card(Net)) and (ord(Net) ne 1))=par(Net)*(1-yc(Net))*CASE+par(Net)*(1-CASE)



                                         solve Rx_OCP_index1 using nlp minimizing zobj;
                                         mstatS5(iter,neigh,iner1)=Rx_OCP_index1.modelstat;
                                         execute_unload "DLR1";

                                                 if(ord(iner1) eq card(iner1),
                                                         if(mstatS5(iter,neigh,iner1) eq 3 or mstatS5(iter,neigh,iner1) eq 4 or mstatS5(iter,neigh,iner1) eq 5 or mstatS5(iter,neigh,iner1) eq 6 or mstatS5(iter,neigh,iner1) eq 11 or mstatS5(iter,neigh,iner1) eq 12 or mstatS5(iter,neigh,iner1) eq 13 or mstatS5(iter,neigh,iner1) eq 14 or mstatS5(iter,neigh,iner1) eq 18 or mstatS5(iter,neigh,iner1) eq 19 ,
                                                         dvs(iter,neigh)=1;
                                                         else
                                                         fplushvalue(iter,neigh)=zobj.l;
                                                         dvs(iter,neigh)=(fplushvalue(iter,neigh)-fvalue(iter))/(hstep);
                                                         );
                                                 );

                                         );

                                 );

                         );

         );

***S4. Local optimality: Stopping criterion verification
dvs(iter,neigh)=round(dvs(iter,neigh),8);
dmin(iter)=smin((neigh),dvs(iter,neigh));
         if(dmin(iter)  ge 0 ,
         stop1=1;
         );
***S5. Steepest descent: Selection of ds
count=0;
         loop(neigh,
                         if(dvs(iter,neigh) eq dmin(iter) and count eq 0,
                         selectd(iter,neigh)=1;
                         count=count+1;
                         );
         );

stop3=0;
         loop(iner2$(stop1=0 and stop3=0),
***S6. Line search: Feasibility of external variables
         stop4=0;
         xvalue_p(extvar)=xvalue(iter,extvar)+(sum(neigh,selectd(iter,neigh)*directions(neigh,extvar)))*hiner2(iner2);

         fex(sf1)=xvalue_p("x1");
         fex(sf2)=xvalue_p("x1")+sum(scat,xvalue_p(scat))+xvalue_p("x2");
         fex(sbup)= xvalue_p("x1")+sum(scat,xvalue_p(scat))+xvalue_p("x2")+xvalue_p("x3");
         loop(scat,
         if(ord(scat) eq 1,
         fex(scat)=xvalue_p("x1")+xvalue_p("x4");
         else
         fex(scat)=xvalue_p(scat)+fex(scat-1);
         );
         );

         gval("g1")=2-xvalue_p("x1");
         gval("g2")=-xvalue_p("x3");
         gval("g3")=-xvalue_p("x4");
         gval("g4")=-xvalue_p("x2");
         gval("g5")=xvalue_p("x1")+sum(scat,xvalue_p(scat))+xvalue_p("x2")+xvalue_p("x3")-44;
         cuent_gineq=2;
         loop(extineq$(ord(extineq) ge 6),
         gval(extineq)=1+(sum(scat$(ord(scat) eq cuent_gineq-1),fex(scat)))-(sum(scat$(ord(scat) eq cuent_gineq),fex(scat)));
         cuent_gineq=cuent_gineq+1;
         );

                 if( sum(extineq$(gval(extineq) gt 0 ),gval(extineq)) gt 0  ,
                 stop4=1;
                 stop3=1;
                 );

                 loop(iner1$(stop4 eq 0),
                         if(ord(iner1) eq 1 and ord(iner2) eq 1 ,
                         execute_loadpoint "DLR";
                         elseif ord(iner1) eq 2 and ord(iner2) eq 1,
                         execute_loadpoint "DLR";
                         else
                         execute_loadpoint "DLR2";
                         );


                           if(ord(iner1) eq 1 ,
***S6. Line search: Computation of the objective function value
                          yaux(extvar,net)=0;
                          yaux(extvar,net)$(ord(net) eq floor(fex(extvar)))=1-mod(fex(extvar),1);
                          yaux(extvar,net)$(ord(net) eq 1+floor(fex(extvar)))=mod(fex(extvar),1);

                          yf(Net,'1')=sum(sf1,yaux(sf1,net));
                          yf(Net,'2')=sum(sf2,yaux(sf2,net));
                          yb(Net)=sum(sbup,yaux(sbup,net));
                          yc(Net)=sum(scat,yaux(scat,net));
                          par(Net)$((ord(Net) ne card(Net)) and (ord(Net) ne 1))=(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yr(Net1)))+(yb(Net))-(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yb(Net1)));
                          ye(Net)$((ord(Net) ne card(Net)) and (ord(Net) ne 1))=par(Net)*(1-yc(Net))*CASE+par(Net)*(1-CASE);

                           solve Rx_OCP_index1 using nlp minimizing zobj;
                           mstatS7(iter,iner2,iner1)=Rx_OCP_index1.modelstat;
                           execute_unload "DLR2";

                                 if( mstatS7(iter,iner2,iner1) eq 3 or mstatS7(iter,iner2,iner1) eq 4 or mstatS7(iter,iner2,iner1) eq 5 or mstatS7(iter,iner2,iner1) eq 6 or mstatS7(iter,iner2,iner1) eq 11 or mstatS7(iter,iner2,iner1) eq 12 or mstatS7(iter,iner2,iner1) eq 13 or mstatS7(iter,iner2,iner1) eq 14 or mstatS7(iter,iner2,iner1) eq 18 or mstatS7(iter,iner2,iner1) eq 19  ,
                                 fvalueiner(iter,iner2)=1E+10;
                                 stop4=0;
                                 else
                                 stop4=1;
                                 fvalueiner(iter,iner2)=zobj.l;
                                 );
***S6. Line search: Computation of the objective function value  with the convergence procedure
                           else
                           xvalue_p(extvar)=xvalue(iter,extvar)+(sum(neigh,selectd(iter,neigh)*directions(neigh,extvar)))*(hiner2(iner2)+hiner1(iner1)-hstep);

                           fex(sf1)=xvalue_p("x1");
                           fex(sf2)=xvalue_p("x1")+sum(scat,xvalue_p(scat))+xvalue_p("x2");
                           fex(sbup)= xvalue_p("x1")+sum(scat,xvalue_p(scat))+xvalue_p("x2")+xvalue_p("x3");
                           loop(scat,
                           if(ord(scat) eq 1,
                           fex(scat)=xvalue_p("x1")+xvalue_p("x4");
                           else
                           fex(scat)=xvalue_p(scat)+fex(scat-1);
                           );
                           );



                           yaux(extvar,net)=0;
                           yaux(extvar,net)$(ord(net) eq floor(fex(extvar)))=1-mod(fex(extvar),1);
                           yaux(extvar,net)$(ord(net) eq 1+floor(fex(extvar)))=mod(fex(extvar),1);

                           yf(Net,'1')=sum(sf1,yaux(sf1,net));
                           yf(Net,'2')=sum(sf2,yaux(sf2,net));
                           yb(Net)=sum(sbup,yaux(sbup,net));
                           yc(Net)=sum(scat,yaux(scat,net));
                           par(Net)$((ord(Net) ne card(Net)) and (ord(Net) ne 1))=(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yr(Net1)))+(yb(Net))-(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yb(Net1)));
                           ye(Net)$((ord(Net) ne card(Net)) and (ord(Net) ne 1))=par(Net)*(1-yc(Net))*CASE+par(Net)*(1-CASE)

                           solve Rx_OCP_index1 using nlp minimizing zobj;
                           mstatS8(iter,iner2,iner1)=Rx_OCP_index1.modelstat;
                           execute_unload "DLR2";

                                 if(ord(iner1) eq card(iner1),
                                 stop4=1;

                                         if(mstatS8(iter,iner2,iner1) eq 3 or mstatS8(iter,iner2,iner1) eq 4 or mstatS8(iter,iner2,iner1) eq 5 or mstatS8(iter,iner2,iner1) eq 6 or mstatS8(iter,iner2,iner1) eq 11 or mstatS8(iter,iner2,iner1) eq 12 or mstatS8(iter,iner2,iner1) eq 13 or mstatS8(iter,iner2,iner1) eq 14 or mstatS8(iter,iner2,iner1) eq 18 or mstatS8(iter,iner2,iner1) eq 19  ,
                                         fvalueiner(iter,iner2)=1E+10;
                                         else
                                         fvalueiner(iter,iner2)=zobj.l;
                                         );

                                 );

                           );
                 );
***S6. Line search: Stopping criterion
                 if(ord(iner2) ne 1 and  fvalueiner(iter,iner2) ge fvalueiner(iter,iner2-1),
                 stop3=1;
                 );
                 if(stop3=0 and stop1=0,
                 objval=zobj.l;
                 execute_unload "DLR";
                 execute_unload "ALG_INDUSTRIAL2";
                         loop(extvar,
                          xvalueselect(iter,extvar)=xvalue(iter,extvar)+(sum(neigh,selectd(iter,neigh)*directions(neigh,extvar)))*hiner2(iner2);
                         );
                 );


         );



***Initial value of x for the next iteration
         if(stop1 =0,
         xvalue(iter+1,extvar)= xvalueselect(iter,extvar);
         else
         xvalue(iter+1,extvar)=xvalue(iter,extvar);
         );



);
CPUtime=timeElapsed-CPUtimeactual;
CPUtimeactual=timeElapsed;
execute_unload "ALG_INDUSTRIAL_INFO2";




