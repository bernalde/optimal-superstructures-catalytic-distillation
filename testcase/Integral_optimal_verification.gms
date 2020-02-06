$ontext
================================================================================
Verificación de optimalidad del optimo s-local encontrado con el codigo
ALG_SEPARABLE.gms
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
set Net "Todas las etapas de la columna reales y sobrantes incluyendo con y rev" /1*22/;
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
parameter FB "Flujo de alimentación de butenos [mol/min]" /5.774/
parameter zb(N,j,comp) "Porcentaje molar en la alimentación de butenos";
zb(N,j,'iButene')=30;
zb(N,j,'nButene')=100-zb(N,j,'iButene');
zb(N,j,'Ethanol')=0;
zb(N,j,'ETBE')=0;

*Alimentacion etanol
parameter
FE  "Flujo de alimentación de etanol [mol-h]" /1.7118/
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
Pop     "Presión de operación condensador [bar]"        /9.5/
TaliB   "Temperatura de alimentación de butenos [K]"    /323/
TaliE   "Temperatura de alimentación etanol [K]"        /342.38/
xBetbe  "Composición molar de ETBE en fondos deseada"   /83/
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
*                           Cálculo de reacción química
*-------------------------------------------------------------------------------
Parameter
Nu(comp) "Coeficientes estequiométricos en la reacción"
/
iButene -1
Ethanol -1
nButene 0
ETBE    1
/
mcat  "Masa del catalizador"     /0.4/
;
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
scalar NCmax "numero maximo de etapas reactivas" /3/;

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
D.up=0.3;
D.lo=0.1;

hw.up=0.1;
hw.lo=0.0001;

hs.up=0.3;
hs.lo=0.1;

htotal.up=10;

at.lo=1e-6;
at.up=0.1;

ad.lo=0.0001;
ad.up=0.01;

lw.up=1;

A0.lo=1e-12;
A0.up=0.01;

RR.lo(N,j)=1;
RR.up(N,j)=10;

Qr.up(N,j)=400;
Qr.lo(N,j)=100;

Qc.up(N,j)=900;

L.up(N,j,Net)=200;
V.up(N,j,Net)=200;

BR.up(N,j)=10;

P.up(N,j,Net)= 10;
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

Krate.up(N,j,Net)=1000000000;
Krate.lo(N,j,Net)=-10;
Krate.scale(N,j,Net)=10000000;

Ka.up(N,j,Net)=100;

Rx.up(N,j,Net)=100;
Rx.lo(N,j,Net)=-100000;
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

depHvib.lo(N,j,comp,Net)=-10;
depHvib.up(N,j,comp,Net)=10;

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

DP.up(N,j,Net)=0.1;

DPq.up(N,j,Net)=0.1;

uhv.up(N,j,Net)=10;
uhv.lo(N,j,Net)=0.4;

hcl.up(N,j,Net)=0.1;
hcl.lo(N,j,Net)=1e-6;

Lload.up(N,j,Net)=0.008333;
Ffactor.up(N,j,Net)=3.575;
dPcat.up(N,j,net)=10;

x.lo(N,j,comp,Net)= 0;
x.up(N,j,comp,Net)= 100;

y.lo(N,j,comp,Net)= 0;
y.up(N,j,comp,Net)= 100;

z.lo(N,j,Net)= 0.5 ;
z.up(N,j,Net)=  1.3  ;

Temp.lo(N,j,Net)= 200;
Temp.up(N,j,Net)= 417.89;

Tcritm.up(N,j,Net)=600;
Tcritm.lo(N,j,Net)=417.9;

bEOS.lo(N,j,Net)=min(biEOS('iButene'),biEOS('Ethanol'),biEOS('nButene'),biEOS('ETBE'));
aEOS.lo(N,j,Net)=1e-6;

*Variables en las restricciones de inundacion que pueden generar problemas con DICPOT:
unv.lo(N,j,Net)=0.01;
unv.up(N,j,Net)=0.89;

ul.up(N,j,Net)=30;
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
option nlp=msnlp;
Rx_OCP_index1.threads=0;
Rx_OCP_index1.SCALEOPT=1;
Rx_OCP_index1.OPTFILE=1;
Rx_OCP_index1.workspace=4000;
Rx_OCP_index1.reslim=3600;


set extvar "variables externas" /x1*x6/;
set extineq "restricciones desigualdad externas" /g1*g7/;
set direc "direccions" /1*728/;

scalar Nboil "etapa de boil up";
scalar kparam "etapas internas" /3/;

set etapabup "etapa de alimentacion de boil up" /9/;
parameter Nboilloop(etapabup)
/
9  9
/  ;
parameter veci1(direc)
/
1        -1
2        -1
3        -1
4        -1
5        -1
6        -1
7        -1
8        -1
9        -1
10        -1
11        -1
12        -1
13        -1
14        -1
15        -1
16        -1
17        -1
18        -1
19        -1
20        -1
21        -1
22        -1
23        -1
24        -1
25        -1
26        -1
27        -1
28        -1
29        -1
30        -1
31        -1
32        -1
33        -1
34        -1
35        -1
36        -1
37        -1
38        -1
39        -1
40        -1
41        -1
42        -1
43        -1
44        -1
45        -1
46        -1
47        -1
48        -1
49        -1
50        -1
51        -1
52        -1
53        -1
54        -1
55        -1
56        -1
57        -1
58        -1
59        -1
60        -1
61        -1
62        -1
63        -1
64        -1
65        -1
66        -1
67        -1
68        -1
69        -1
70        -1
71        -1
72        -1
73        -1
74        -1
75        -1
76        -1
77        -1
78        -1
79        -1
80        -1
81        -1
82        -1
83        -1
84        -1
85        -1
86        -1
87        -1
88        -1
89        -1
90        -1
91        -1
92        -1
93        -1
94        -1
95        -1
96        -1
97        -1
98        -1
99        -1
100        -1
101        -1
102        -1
103        -1
104        -1
105        -1
106        -1
107        -1
108        -1
109        -1
110        -1
111        -1
112        -1
113        -1
114        -1
115        -1
116        -1
117        -1
118        -1
119        -1
120        -1
121        -1
122        -1
123        -1
124        -1
125        -1
126        -1
127        -1
128        -1
129        -1
130        -1
131        -1
132        -1
133        -1
134        -1
135        -1
136        -1
137        -1
138        -1
139        -1
140        -1
141        -1
142        -1
143        -1
144        -1
145        -1
146        -1
147        -1
148        -1
149        -1
150        -1
151        -1
152        -1
153        -1
154        -1
155        -1
156        -1
157        -1
158        -1
159        -1
160        -1
161        -1
162        -1
163        -1
164        -1
165        -1
166        -1
167        -1
168        -1
169        -1
170        -1
171        -1
172        -1
173        -1
174        -1
175        -1
176        -1
177        -1
178        -1
179        -1
180        -1
181        -1
182        -1
183        -1
184        -1
185        -1
186        -1
187        -1
188        -1
189        -1
190        -1
191        -1
192        -1
193        -1
194        -1
195        -1
196        -1
197        -1
198        -1
199        -1
200        -1
201        -1
202        -1
203        -1
204        -1
205        -1
206        -1
207        -1
208        -1
209        -1
210        -1
211        -1
212        -1
213        -1
214        -1
215        -1
216        -1
217        -1
218        -1
219        -1
220        -1
221        -1
222        -1
223        -1
224        -1
225        -1
226        -1
227        -1
228        -1
229        -1
230        -1
231        -1
232        -1
233        -1
234        -1
235        -1
236        -1
237        -1
238        -1
239        -1
240        -1
241        -1
242        -1
243        -1
244        0
245        0
246        0
247        0
248        0
249        0
250        0
251        0
252        0
253        0
254        0
255        0
256        0
257        0
258        0
259        0
260        0
261        0
262        0
263        0
264        0
265        0
266        0
267        0
268        0
269        0
270        0
271        0
272        0
273        0
274        0
275        0
276        0
277        0
278        0
279        0
280        0
281        0
282        0
283        0
284        0
285        0
286        0
287        0
288        0
289        0
290        0
291        0
292        0
293        0
294        0
295        0
296        0
297        0
298        0
299        0
300        0
301        0
302        0
303        0
304        0
305        0
306        0
307        0
308        0
309        0
310        0
311        0
312        0
313        0
314        0
315        0
316        0
317        0
318        0
319        0
320        0
321        0
322        0
323        0
324        0
325        0
326        0
327        0
328        0
329        0
330        0
331        0
332        0
333        0
334        0
335        0
336        0
337        0
338        0
339        0
340        0
341        0
342        0
343        0
344        0
345        0
346        0
347        0
348        0
349        0
350        0
351        0
352        0
353        0
354        0
355        0
356        0
357        0
358        0
359        0
360        0
361        0
362        0
363        0
364        0
365        0
366        0
367        0
368        0
369        0
370        0
371        0
372        0
373        0
374        0
375        0
376        0
377        0
378        0
379        0
380        0
381        0
382        0
383        0
384        0
385        0
386        0
387        0
388        0
389        0
390        0
391        0
392        0
393        0
394        0
395        0
396        0
397        0
398        0
399        0
400        0
401        0
402        0
403        0
404        0
405        0
406        0
407        0
408        0
409        0
410        0
411        0
412        0
413        0
414        0
415        0
416        0
417        0
418        0
419        0
420        0
421        0
422        0
423        0
424        0
425        0
426        0
427        0
428        0
429        0
430        0
431        0
432        0
433        0
434        0
435        0
436        0
437        0
438        0
439        0
440        0
441        0
442        0
443        0
444        0
445        0
446        0
447        0
448        0
449        0
450        0
451        0
452        0
453        0
454        0
455        0
456        0
457        0
458        0
459        0
460        0
461        0
462        0
463        0
464        0
465        0
466        0
467        0
468        0
469        0
470        0
471        0
472        0
473        0
474        0
475        0
476        0
477        0
478        0
479        0
480        0
481        0
482        0
483        0
484        0
485        0
486        1
487        1
488        1
489        1
490        1
491        1
492        1
493        1
494        1
495        1
496        1
497        1
498        1
499        1
500        1
501        1
502        1
503        1
504        1
505        1
506        1
507        1
508        1
509        1
510        1
511        1
512        1
513        1
514        1
515        1
516        1
517        1
518        1
519        1
520        1
521        1
522        1
523        1
524        1
525        1
526        1
527        1
528        1
529        1
530        1
531        1
532        1
533        1
534        1
535        1
536        1
537        1
538        1
539        1
540        1
541        1
542        1
543        1
544        1
545        1
546        1
547        1
548        1
549        1
550        1
551        1
552        1
553        1
554        1
555        1
556        1
557        1
558        1
559        1
560        1
561        1
562        1
563        1
564        1
565        1
566        1
567        1
568        1
569        1
570        1
571        1
572        1
573        1
574        1
575        1
576        1
577        1
578        1
579        1
580        1
581        1
582        1
583        1
584        1
585        1
586        1
587        1
588        1
589        1
590        1
591        1
592        1
593        1
594        1
595        1
596        1
597        1
598        1
599        1
600        1
601        1
602        1
603        1
604        1
605        1
606        1
607        1
608        1
609        1
610        1
611        1
612        1
613        1
614        1
615        1
616        1
617        1
618        1
619        1
620        1
621        1
622        1
623        1
624        1
625        1
626        1
627        1
628        1
629        1
630        1
631        1
632        1
633        1
634        1
635        1
636        1
637        1
638        1
639        1
640        1
641        1
642        1
643        1
644        1
645        1
646        1
647        1
648        1
649        1
650        1
651        1
652        1
653        1
654        1
655        1
656        1
657        1
658        1
659        1
660        1
661        1
662        1
663        1
664        1
665        1
666        1
667        1
668        1
669        1
670        1
671        1
672        1
673        1
674        1
675        1
676        1
677        1
678        1
679        1
680        1
681        1
682        1
683        1
684        1
685        1
686        1
687        1
688        1
689        1
690        1
691        1
692        1
693        1
694        1
695        1
696        1
697        1
698        1
699        1
700        1
701        1
702        1
703        1
704        1
705        1
706        1
707        1
708        1
709        1
710        1
711        1
712        1
713        1
714        1
715        1
716        1
717        1
718        1
719        1
720        1
721        1
722        1
723        1
724        1
725        1
726        1
727        1
728        1


/;
parameter veci2(direc)
/
1        -1
2        -1
3        -1
4        -1
5        -1
6        -1
7        -1
8        -1
9        -1
10        -1
11        -1
12        -1
13        -1
14        -1
15        -1
16        -1
17        -1
18        -1
19        -1
20        -1
21        -1
22        -1
23        -1
24        -1
25        -1
26        -1
27        -1
28        -1
29        -1
30        -1
31        -1
32        -1
33        -1
34        -1
35        -1
36        -1
37        -1
38        -1
39        -1
40        -1
41        -1
42        -1
43        -1
44        -1
45        -1
46        -1
47        -1
48        -1
49        -1
50        -1
51        -1
52        -1
53        -1
54        -1
55        -1
56        -1
57        -1
58        -1
59        -1
60        -1
61        -1
62        -1
63        -1
64        -1
65        -1
66        -1
67        -1
68        -1
69        -1
70        -1
71        -1
72        -1
73        -1
74        -1
75        -1
76        -1
77        -1
78        -1
79        -1
80        -1
81        -1
82        0
83        0
84        0
85        0
86        0
87        0
88        0
89        0
90        0
91        0
92        0
93        0
94        0
95        0
96        0
97        0
98        0
99        0
100        0
101        0
102        0
103        0
104        0
105        0
106        0
107        0
108        0
109        0
110        0
111        0
112        0
113        0
114        0
115        0
116        0
117        0
118        0
119        0
120        0
121        0
122        0
123        0
124        0
125        0
126        0
127        0
128        0
129        0
130        0
131        0
132        0
133        0
134        0
135        0
136        0
137        0
138        0
139        0
140        0
141        0
142        0
143        0
144        0
145        0
146        0
147        0
148        0
149        0
150        0
151        0
152        0
153        0
154        0
155        0
156        0
157        0
158        0
159        0
160        0
161        0
162        0
163        1
164        1
165        1
166        1
167        1
168        1
169        1
170        1
171        1
172        1
173        1
174        1
175        1
176        1
177        1
178        1
179        1
180        1
181        1
182        1
183        1
184        1
185        1
186        1
187        1
188        1
189        1
190        1
191        1
192        1
193        1
194        1
195        1
196        1
197        1
198        1
199        1
200        1
201        1
202        1
203        1
204        1
205        1
206        1
207        1
208        1
209        1
210        1
211        1
212        1
213        1
214        1
215        1
216        1
217        1
218        1
219        1
220        1
221        1
222        1
223        1
224        1
225        1
226        1
227        1
228        1
229        1
230        1
231        1
232        1
233        1
234        1
235        1
236        1
237        1
238        1
239        1
240        1
241        1
242        1
243        1
244        -1
245        -1
246        -1
247        -1
248        -1
249        -1
250        -1
251        -1
252        -1
253        -1
254        -1
255        -1
256        -1
257        -1
258        -1
259        -1
260        -1
261        -1
262        -1
263        -1
264        -1
265        -1
266        -1
267        -1
268        -1
269        -1
270        -1
271        -1
272        -1
273        -1
274        -1
275        -1
276        -1
277        -1
278        -1
279        -1
280        -1
281        -1
282        -1
283        -1
284        -1
285        -1
286        -1
287        -1
288        -1
289        -1
290        -1
291        -1
292        -1
293        -1
294        -1
295        -1
296        -1
297        -1
298        -1
299        -1
300        -1
301        -1
302        -1
303        -1
304        -1
305        -1
306        -1
307        -1
308        -1
309        -1
310        -1
311        -1
312        -1
313        -1
314        -1
315        -1
316        -1
317        -1
318        -1
319        -1
320        -1
321        -1
322        -1
323        -1
324        -1
325        0
326        0
327        0
328        0
329        0
330        0
331        0
332        0
333        0
334        0
335        0
336        0
337        0
338        0
339        0
340        0
341        0
342        0
343        0
344        0
345        0
346        0
347        0
348        0
349        0
350        0
351        0
352        0
353        0
354        0
355        0
356        0
357        0
358        0
359        0
360        0
361        0
362        0
363        0
364        0
365        0
366        0
367        0
368        0
369        0
370        0
371        0
372        0
373        0
374        0
375        0
376        0
377        0
378        0
379        0
380        0
381        0
382        0
383        0
384        0
385        0
386        0
387        0
388        0
389        0
390        0
391        0
392        0
393        0
394        0
395        0
396        0
397        0
398        0
399        0
400        0
401        0
402        0
403        0
404        0
405        1
406        1
407        1
408        1
409        1
410        1
411        1
412        1
413        1
414        1
415        1
416        1
417        1
418        1
419        1
420        1
421        1
422        1
423        1
424        1
425        1
426        1
427        1
428        1
429        1
430        1
431        1
432        1
433        1
434        1
435        1
436        1
437        1
438        1
439        1
440        1
441        1
442        1
443        1
444        1
445        1
446        1
447        1
448        1
449        1
450        1
451        1
452        1
453        1
454        1
455        1
456        1
457        1
458        1
459        1
460        1
461        1
462        1
463        1
464        1
465        1
466        1
467        1
468        1
469        1
470        1
471        1
472        1
473        1
474        1
475        1
476        1
477        1
478        1
479        1
480        1
481        1
482        1
483        1
484        1
485        1
486        -1
487        -1
488        -1
489        -1
490        -1
491        -1
492        -1
493        -1
494        -1
495        -1
496        -1
497        -1
498        -1
499        -1
500        -1
501        -1
502        -1
503        -1
504        -1
505        -1
506        -1
507        -1
508        -1
509        -1
510        -1
511        -1
512        -1
513        -1
514        -1
515        -1
516        -1
517        -1
518        -1
519        -1
520        -1
521        -1
522        -1
523        -1
524        -1
525        -1
526        -1
527        -1
528        -1
529        -1
530        -1
531        -1
532        -1
533        -1
534        -1
535        -1
536        -1
537        -1
538        -1
539        -1
540        -1
541        -1
542        -1
543        -1
544        -1
545        -1
546        -1
547        -1
548        -1
549        -1
550        -1
551        -1
552        -1
553        -1
554        -1
555        -1
556        -1
557        -1
558        -1
559        -1
560        -1
561        -1
562        -1
563        -1
564        -1
565        -1
566        -1
567        0
568        0
569        0
570        0
571        0
572        0
573        0
574        0
575        0
576        0
577        0
578        0
579        0
580        0
581        0
582        0
583        0
584        0
585        0
586        0
587        0
588        0
589        0
590        0
591        0
592        0
593        0
594        0
595        0
596        0
597        0
598        0
599        0
600        0
601        0
602        0
603        0
604        0
605        0
606        0
607        0
608        0
609        0
610        0
611        0
612        0
613        0
614        0
615        0
616        0
617        0
618        0
619        0
620        0
621        0
622        0
623        0
624        0
625        0
626        0
627        0
628        0
629        0
630        0
631        0
632        0
633        0
634        0
635        0
636        0
637        0
638        0
639        0
640        0
641        0
642        0
643        0
644        0
645        0
646        0
647        0
648        1
649        1
650        1
651        1
652        1
653        1
654        1
655        1
656        1
657        1
658        1
659        1
660        1
661        1
662        1
663        1
664        1
665        1
666        1
667        1
668        1
669        1
670        1
671        1
672        1
673        1
674        1
675        1
676        1
677        1
678        1
679        1
680        1
681        1
682        1
683        1
684        1
685        1
686        1
687        1
688        1
689        1
690        1
691        1
692        1
693        1
694        1
695        1
696        1
697        1
698        1
699        1
700        1
701        1
702        1
703        1
704        1
705        1
706        1
707        1
708        1
709        1
710        1
711        1
712        1
713        1
714        1
715        1
716        1
717        1
718        1
719        1
720        1
721        1
722        1
723        1
724        1
725        1
726        1
727        1
728        1


/

;
parameter veci3(direc)
/
1        -1
2        -1
3        -1
4        -1
5        -1
6        -1
7        -1
8        -1
9        -1
10        -1
11        -1
12        -1
13        -1
14        -1
15        -1
16        -1
17        -1
18        -1
19        -1
20        -1
21        -1
22        -1
23        -1
24        -1
25        -1
26        -1
27        -1
28        0
29        0
30        0
31        0
32        0
33        0
34        0
35        0
36        0
37        0
38        0
39        0
40        0
41        0
42        0
43        0
44        0
45        0
46        0
47        0
48        0
49        0
50        0
51        0
52        0
53        0
54        0
55        1
56        1
57        1
58        1
59        1
60        1
61        1
62        1
63        1
64        1
65        1
66        1
67        1
68        1
69        1
70        1
71        1
72        1
73        1
74        1
75        1
76        1
77        1
78        1
79        1
80        1
81        1
82        -1
83        -1
84        -1
85        -1
86        -1
87        -1
88        -1
89        -1
90        -1
91        -1
92        -1
93        -1
94        -1
95        -1
96        -1
97        -1
98        -1
99        -1
100        -1
101        -1
102        -1
103        -1
104        -1
105        -1
106        -1
107        -1
108        -1
109        0
110        0
111        0
112        0
113        0
114        0
115        0
116        0
117        0
118        0
119        0
120        0
121        0
122        0
123        0
124        0
125        0
126        0
127        0
128        0
129        0
130        0
131        0
132        0
133        0
134        0
135        0
136        1
137        1
138        1
139        1
140        1
141        1
142        1
143        1
144        1
145        1
146        1
147        1
148        1
149        1
150        1
151        1
152        1
153        1
154        1
155        1
156        1
157        1
158        1
159        1
160        1
161        1
162        1
163        -1
164        -1
165        -1
166        -1
167        -1
168        -1
169        -1
170        -1
171        -1
172        -1
173        -1
174        -1
175        -1
176        -1
177        -1
178        -1
179        -1
180        -1
181        -1
182        -1
183        -1
184        -1
185        -1
186        -1
187        -1
188        -1
189        -1
190        0
191        0
192        0
193        0
194        0
195        0
196        0
197        0
198        0
199        0
200        0
201        0
202        0
203        0
204        0
205        0
206        0
207        0
208        0
209        0
210        0
211        0
212        0
213        0
214        0
215        0
216        0
217        1
218        1
219        1
220        1
221        1
222        1
223        1
224        1
225        1
226        1
227        1
228        1
229        1
230        1
231        1
232        1
233        1
234        1
235        1
236        1
237        1
238        1
239        1
240        1
241        1
242        1
243        1
244        -1
245        -1
246        -1
247        -1
248        -1
249        -1
250        -1
251        -1
252        -1
253        -1
254        -1
255        -1
256        -1
257        -1
258        -1
259        -1
260        -1
261        -1
262        -1
263        -1
264        -1
265        -1
266        -1
267        -1
268        -1
269        -1
270        -1
271        0
272        0
273        0
274        0
275        0
276        0
277        0
278        0
279        0
280        0
281        0
282        0
283        0
284        0
285        0
286        0
287        0
288        0
289        0
290        0
291        0
292        0
293        0
294        0
295        0
296        0
297        0
298        1
299        1
300        1
301        1
302        1
303        1
304        1
305        1
306        1
307        1
308        1
309        1
310        1
311        1
312        1
313        1
314        1
315        1
316        1
317        1
318        1
319        1
320        1
321        1
322        1
323        1
324        1
325        -1
326        -1
327        -1
328        -1
329        -1
330        -1
331        -1
332        -1
333        -1
334        -1
335        -1
336        -1
337        -1
338        -1
339        -1
340        -1
341        -1
342        -1
343        -1
344        -1
345        -1
346        -1
347        -1
348        -1
349        -1
350        -1
351        -1
352        0
353        0
354        0
355        0
356        0
357        0
358        0
359        0
360        0
361        0
362        0
363        0
364        0
365        0
366        0
367        0
368        0
369        0
370        0
371        0
372        0
373        0
374        0
375        0
376        0
377        0
378        1
379        1
380        1
381        1
382        1
383        1
384        1
385        1
386        1
387        1
388        1
389        1
390        1
391        1
392        1
393        1
394        1
395        1
396        1
397        1
398        1
399        1
400        1
401        1
402        1
403        1
404        1
405        -1
406        -1
407        -1
408        -1
409        -1
410        -1
411        -1
412        -1
413        -1
414        -1
415        -1
416        -1
417        -1
418        -1
419        -1
420        -1
421        -1
422        -1
423        -1
424        -1
425        -1
426        -1
427        -1
428        -1
429        -1
430        -1
431        -1
432        0
433        0
434        0
435        0
436        0
437        0
438        0
439        0
440        0
441        0
442        0
443        0
444        0
445        0
446        0
447        0
448        0
449        0
450        0
451        0
452        0
453        0
454        0
455        0
456        0
457        0
458        0
459        1
460        1
461        1
462        1
463        1
464        1
465        1
466        1
467        1
468        1
469        1
470        1
471        1
472        1
473        1
474        1
475        1
476        1
477        1
478        1
479        1
480        1
481        1
482        1
483        1
484        1
485        1
486        -1
487        -1
488        -1
489        -1
490        -1
491        -1
492        -1
493        -1
494        -1
495        -1
496        -1
497        -1
498        -1
499        -1
500        -1
501        -1
502        -1
503        -1
504        -1
505        -1
506        -1
507        -1
508        -1
509        -1
510        -1
511        -1
512        -1
513        0
514        0
515        0
516        0
517        0
518        0
519        0
520        0
521        0
522        0
523        0
524        0
525        0
526        0
527        0
528        0
529        0
530        0
531        0
532        0
533        0
534        0
535        0
536        0
537        0
538        0
539        0
540        1
541        1
542        1
543        1
544        1
545        1
546        1
547        1
548        1
549        1
550        1
551        1
552        1
553        1
554        1
555        1
556        1
557        1
558        1
559        1
560        1
561        1
562        1
563        1
564        1
565        1
566        1
567        -1
568        -1
569        -1
570        -1
571        -1
572        -1
573        -1
574        -1
575        -1
576        -1
577        -1
578        -1
579        -1
580        -1
581        -1
582        -1
583        -1
584        -1
585        -1
586        -1
587        -1
588        -1
589        -1
590        -1
591        -1
592        -1
593        -1
594        0
595        0
596        0
597        0
598        0
599        0
600        0
601        0
602        0
603        0
604        0
605        0
606        0
607        0
608        0
609        0
610        0
611        0
612        0
613        0
614        0
615        0
616        0
617        0
618        0
619        0
620        0
621        1
622        1
623        1
624        1
625        1
626        1
627        1
628        1
629        1
630        1
631        1
632        1
633        1
634        1
635        1
636        1
637        1
638        1
639        1
640        1
641        1
642        1
643        1
644        1
645        1
646        1
647        1
648        -1
649        -1
650        -1
651        -1
652        -1
653        -1
654        -1
655        -1
656        -1
657        -1
658        -1
659        -1
660        -1
661        -1
662        -1
663        -1
664        -1
665        -1
666        -1
667        -1
668        -1
669        -1
670        -1
671        -1
672        -1
673        -1
674        -1
675        0
676        0
677        0
678        0
679        0
680        0
681        0
682        0
683        0
684        0
685        0
686        0
687        0
688        0
689        0
690        0
691        0
692        0
693        0
694        0
695        0
696        0
697        0
698        0
699        0
700        0
701        0
702        1
703        1
704        1
705        1
706        1
707        1
708        1
709        1
710        1
711        1
712        1
713        1
714        1
715        1
716        1
717        1
718        1
719        1
720        1
721        1
722        1
723        1
724        1
725        1
726        1
727        1
728        1


/
;

parameter veci4(direc)
/
1        -1
2        -1
3        -1
4        -1
5        -1
6        -1
7        -1
8        -1
9        -1
10        0
11        0
12        0
13        0
14        0
15        0
16        0
17        0
18        0
19        1
20        1
21        1
22        1
23        1
24        1
25        1
26        1
27        1
28        -1
29        -1
30        -1
31        -1
32        -1
33        -1
34        -1
35        -1
36        -1
37        0
38        0
39        0
40        0
41        0
42        0
43        0
44        0
45        0
46        1
47        1
48        1
49        1
50        1
51        1
52        1
53        1
54        1
55        -1
56        -1
57        -1
58        -1
59        -1
60        -1
61        -1
62        -1
63        -1
64        0
65        0
66        0
67        0
68        0
69        0
70        0
71        0
72        0
73        1
74        1
75        1
76        1
77        1
78        1
79        1
80        1
81        1
82        -1
83        -1
84        -1
85        -1
86        -1
87        -1
88        -1
89        -1
90        -1
91        0
92        0
93        0
94        0
95        0
96        0
97        0
98        0
99        0
100        1
101        1
102        1
103        1
104        1
105        1
106        1
107        1
108        1
109        -1
110        -1
111        -1
112        -1
113        -1
114        -1
115        -1
116        -1
117        -1
118        0
119        0
120        0
121        0
122        0
123        0
124        0
125        0
126        0
127        1
128        1
129        1
130        1
131        1
132        1
133        1
134        1
135        1
136        -1
137        -1
138        -1
139        -1
140        -1
141        -1
142        -1
143        -1
144        -1
145        0
146        0
147        0
148        0
149        0
150        0
151        0
152        0
153        0
154        1
155        1
156        1
157        1
158        1
159        1
160        1
161        1
162        1
163        -1
164        -1
165        -1
166        -1
167        -1
168        -1
169        -1
170        -1
171        -1
172        0
173        0
174        0
175        0
176        0
177        0
178        0
179        0
180        0
181        1
182        1
183        1
184        1
185        1
186        1
187        1
188        1
189        1
190        -1
191        -1
192        -1
193        -1
194        -1
195        -1
196        -1
197        -1
198        -1
199        0
200        0
201        0
202        0
203        0
204        0
205        0
206        0
207        0
208        1
209        1
210        1
211        1
212        1
213        1
214        1
215        1
216        1
217        -1
218        -1
219        -1
220        -1
221        -1
222        -1
223        -1
224        -1
225        -1
226        0
227        0
228        0
229        0
230        0
231        0
232        0
233        0
234        0
235        1
236        1
237        1
238        1
239        1
240        1
241        1
242        1
243        1
244        -1
245        -1
246        -1
247        -1
248        -1
249        -1
250        -1
251        -1
252        -1
253        0
254        0
255        0
256        0
257        0
258        0
259        0
260        0
261        0
262        1
263        1
264        1
265        1
266        1
267        1
268        1
269        1
270        1
271        -1
272        -1
273        -1
274        -1
275        -1
276        -1
277        -1
278        -1
279        -1
280        0
281        0
282        0
283        0
284        0
285        0
286        0
287        0
288        0
289        1
290        1
291        1
292        1
293        1
294        1
295        1
296        1
297        1
298        -1
299        -1
300        -1
301        -1
302        -1
303        -1
304        -1
305        -1
306        -1
307        0
308        0
309        0
310        0
311        0
312        0
313        0
314        0
315        0
316        1
317        1
318        1
319        1
320        1
321        1
322        1
323        1
324        1
325        -1
326        -1
327        -1
328        -1
329        -1
330        -1
331        -1
332        -1
333        -1
334        0
335        0
336        0
337        0
338        0
339        0
340        0
341        0
342        0
343        1
344        1
345        1
346        1
347        1
348        1
349        1
350        1
351        1
352        -1
353        -1
354        -1
355        -1
356        -1
357        -1
358        -1
359        -1
360        -1
361        0
362        0
363        0
364        0
365        0
366        0
367        0
368        0
369        1
370        1
371        1
372        1
373        1
374        1
375        1
376        1
377        1
378        -1
379        -1
380        -1
381        -1
382        -1
383        -1
384        -1
385        -1
386        -1
387        0
388        0
389        0
390        0
391        0
392        0
393        0
394        0
395        0
396        1
397        1
398        1
399        1
400        1
401        1
402        1
403        1
404        1
405        -1
406        -1
407        -1
408        -1
409        -1
410        -1
411        -1
412        -1
413        -1
414        0
415        0
416        0
417        0
418        0
419        0
420        0
421        0
422        0
423        1
424        1
425        1
426        1
427        1
428        1
429        1
430        1
431        1
432        -1
433        -1
434        -1
435        -1
436        -1
437        -1
438        -1
439        -1
440        -1
441        0
442        0
443        0
444        0
445        0
446        0
447        0
448        0
449        0
450        1
451        1
452        1
453        1
454        1
455        1
456        1
457        1
458        1
459        -1
460        -1
461        -1
462        -1
463        -1
464        -1
465        -1
466        -1
467        -1
468        0
469        0
470        0
471        0
472        0
473        0
474        0
475        0
476        0
477        1
478        1
479        1
480        1
481        1
482        1
483        1
484        1
485        1
486        -1
487        -1
488        -1
489        -1
490        -1
491        -1
492        -1
493        -1
494        -1
495        0
496        0
497        0
498        0
499        0
500        0
501        0
502        0
503        0
504        1
505        1
506        1
507        1
508        1
509        1
510        1
511        1
512        1
513        -1
514        -1
515        -1
516        -1
517        -1
518        -1
519        -1
520        -1
521        -1
522        0
523        0
524        0
525        0
526        0
527        0
528        0
529        0
530        0
531        1
532        1
533        1
534        1
535        1
536        1
537        1
538        1
539        1
540        -1
541        -1
542        -1
543        -1
544        -1
545        -1
546        -1
547        -1
548        -1
549        0
550        0
551        0
552        0
553        0
554        0
555        0
556        0
557        0
558        1
559        1
560        1
561        1
562        1
563        1
564        1
565        1
566        1
567        -1
568        -1
569        -1
570        -1
571        -1
572        -1
573        -1
574        -1
575        -1
576        0
577        0
578        0
579        0
580        0
581        0
582        0
583        0
584        0
585        1
586        1
587        1
588        1
589        1
590        1
591        1
592        1
593        1
594        -1
595        -1
596        -1
597        -1
598        -1
599        -1
600        -1
601        -1
602        -1
603        0
604        0
605        0
606        0
607        0
608        0
609        0
610        0
611        0
612        1
613        1
614        1
615        1
616        1
617        1
618        1
619        1
620        1
621        -1
622        -1
623        -1
624        -1
625        -1
626        -1
627        -1
628        -1
629        -1
630        0
631        0
632        0
633        0
634        0
635        0
636        0
637        0
638        0
639        1
640        1
641        1
642        1
643        1
644        1
645        1
646        1
647        1
648        -1
649        -1
650        -1
651        -1
652        -1
653        -1
654        -1
655        -1
656        -1
657        0
658        0
659        0
660        0
661        0
662        0
663        0
664        0
665        0
666        1
667        1
668        1
669        1
670        1
671        1
672        1
673        1
674        1
675        -1
676        -1
677        -1
678        -1
679        -1
680        -1
681        -1
682        -1
683        -1
684        0
685        0
686        0
687        0
688        0
689        0
690        0
691        0
692        0
693        1
694        1
695        1
696        1
697        1
698        1
699        1
700        1
701        1
702        -1
703        -1
704        -1
705        -1
706        -1
707        -1
708        -1
709        -1
710        -1
711        0
712        0
713        0
714        0
715        0
716        0
717        0
718        0
719        0
720        1
721        1
722        1
723        1
724        1
725        1
726        1
727        1
728        1


/
;

parameter veci5(direc)
/
1        -1
2        -1
3        -1
4        0
5        0
6        0
7        1
8        1
9        1
10        -1
11        -1
12        -1
13        0
14        0
15        0
16        1
17        1
18        1
19        -1
20        -1
21        -1
22        0
23        0
24        0
25        1
26        1
27        1
28        -1
29        -1
30        -1
31        0
32        0
33        0
34        1
35        1
36        1
37        -1
38        -1
39        -1
40        0
41        0
42        0
43        1
44        1
45        1
46        -1
47        -1
48        -1
49        0
50        0
51        0
52        1
53        1
54        1
55        -1
56        -1
57        -1
58        0
59        0
60        0
61        1
62        1
63        1
64        -1
65        -1
66        -1
67        0
68        0
69        0
70        1
71        1
72        1
73        -1
74        -1
75        -1
76        0
77        0
78        0
79        1
80        1
81        1
82        -1
83        -1
84        -1
85        0
86        0
87        0
88        1
89        1
90        1
91        -1
92        -1
93        -1
94        0
95        0
96        0
97        1
98        1
99        1
100        -1
101        -1
102        -1
103        0
104        0
105        0
106        1
107        1
108        1
109        -1
110        -1
111        -1
112        0
113        0
114        0
115        1
116        1
117        1
118        -1
119        -1
120        -1
121        0
122        0
123        0
124        1
125        1
126        1
127        -1
128        -1
129        -1
130        0
131        0
132        0
133        1
134        1
135        1
136        -1
137        -1
138        -1
139        0
140        0
141        0
142        1
143        1
144        1
145        -1
146        -1
147        -1
148        0
149        0
150        0
151        1
152        1
153        1
154        -1
155        -1
156        -1
157        0
158        0
159        0
160        1
161        1
162        1
163        -1
164        -1
165        -1
166        0
167        0
168        0
169        1
170        1
171        1
172        -1
173        -1
174        -1
175        0
176        0
177        0
178        1
179        1
180        1
181        -1
182        -1
183        -1
184        0
185        0
186        0
187        1
188        1
189        1
190        -1
191        -1
192        -1
193        0
194        0
195        0
196        1
197        1
198        1
199        -1
200        -1
201        -1
202        0
203        0
204        0
205        1
206        1
207        1
208        -1
209        -1
210        -1
211        0
212        0
213        0
214        1
215        1
216        1
217        -1
218        -1
219        -1
220        0
221        0
222        0
223        1
224        1
225        1
226        -1
227        -1
228        -1
229        0
230        0
231        0
232        1
233        1
234        1
235        -1
236        -1
237        -1
238        0
239        0
240        0
241        1
242        1
243        1
244        -1
245        -1
246        -1
247        0
248        0
249        0
250        1
251        1
252        1
253        -1
254        -1
255        -1
256        0
257        0
258        0
259        1
260        1
261        1
262        -1
263        -1
264        -1
265        0
266        0
267        0
268        1
269        1
270        1
271        -1
272        -1
273        -1
274        0
275        0
276        0
277        1
278        1
279        1
280        -1
281        -1
282        -1
283        0
284        0
285        0
286        1
287        1
288        1
289        -1
290        -1
291        -1
292        0
293        0
294        0
295        1
296        1
297        1
298        -1
299        -1
300        -1
301        0
302        0
303        0
304        1
305        1
306        1
307        -1
308        -1
309        -1
310        0
311        0
312        0
313        1
314        1
315        1
316        -1
317        -1
318        -1
319        0
320        0
321        0
322        1
323        1
324        1
325        -1
326        -1
327        -1
328        0
329        0
330        0
331        1
332        1
333        1
334        -1
335        -1
336        -1
337        0
338        0
339        0
340        1
341        1
342        1
343        -1
344        -1
345        -1
346        0
347        0
348        0
349        1
350        1
351        1
352        -1
353        -1
354        -1
355        0
356        0
357        0
358        1
359        1
360        1
361        -1
362        -1
363        -1
364        0
365        0
366        1
367        1
368        1
369        -1
370        -1
371        -1
372        0
373        0
374        0
375        1
376        1
377        1
378        -1
379        -1
380        -1
381        0
382        0
383        0
384        1
385        1
386        1
387        -1
388        -1
389        -1
390        0
391        0
392        0
393        1
394        1
395        1
396        -1
397        -1
398        -1
399        0
400        0
401        0
402        1
403        1
404        1
405        -1
406        -1
407        -1
408        0
409        0
410        0
411        1
412        1
413        1
414        -1
415        -1
416        -1
417        0
418        0
419        0
420        1
421        1
422        1
423        -1
424        -1
425        -1
426        0
427        0
428        0
429        1
430        1
431        1
432        -1
433        -1
434        -1
435        0
436        0
437        0
438        1
439        1
440        1
441        -1
442        -1
443        -1
444        0
445        0
446        0
447        1
448        1
449        1
450        -1
451        -1
452        -1
453        0
454        0
455        0
456        1
457        1
458        1
459        -1
460        -1
461        -1
462        0
463        0
464        0
465        1
466        1
467        1
468        -1
469        -1
470        -1
471        0
472        0
473        0
474        1
475        1
476        1
477        -1
478        -1
479        -1
480        0
481        0
482        0
483        1
484        1
485        1
486        -1
487        -1
488        -1
489        0
490        0
491        0
492        1
493        1
494        1
495        -1
496        -1
497        -1
498        0
499        0
500        0
501        1
502        1
503        1
504        -1
505        -1
506        -1
507        0
508        0
509        0
510        1
511        1
512        1
513        -1
514        -1
515        -1
516        0
517        0
518        0
519        1
520        1
521        1
522        -1
523        -1
524        -1
525        0
526        0
527        0
528        1
529        1
530        1
531        -1
532        -1
533        -1
534        0
535        0
536        0
537        1
538        1
539        1
540        -1
541        -1
542        -1
543        0
544        0
545        0
546        1
547        1
548        1
549        -1
550        -1
551        -1
552        0
553        0
554        0
555        1
556        1
557        1
558        -1
559        -1
560        -1
561        0
562        0
563        0
564        1
565        1
566        1
567        -1
568        -1
569        -1
570        0
571        0
572        0
573        1
574        1
575        1
576        -1
577        -1
578        -1
579        0
580        0
581        0
582        1
583        1
584        1
585        -1
586        -1
587        -1
588        0
589        0
590        0
591        1
592        1
593        1
594        -1
595        -1
596        -1
597        0
598        0
599        0
600        1
601        1
602        1
603        -1
604        -1
605        -1
606        0
607        0
608        0
609        1
610        1
611        1
612        -1
613        -1
614        -1
615        0
616        0
617        0
618        1
619        1
620        1
621        -1
622        -1
623        -1
624        0
625        0
626        0
627        1
628        1
629        1
630        -1
631        -1
632        -1
633        0
634        0
635        0
636        1
637        1
638        1
639        -1
640        -1
641        -1
642        0
643        0
644        0
645        1
646        1
647        1
648        -1
649        -1
650        -1
651        0
652        0
653        0
654        1
655        1
656        1
657        -1
658        -1
659        -1
660        0
661        0
662        0
663        1
664        1
665        1
666        -1
667        -1
668        -1
669        0
670        0
671        0
672        1
673        1
674        1
675        -1
676        -1
677        -1
678        0
679        0
680        0
681        1
682        1
683        1
684        -1
685        -1
686        -1
687        0
688        0
689        0
690        1
691        1
692        1
693        -1
694        -1
695        -1
696        0
697        0
698        0
699        1
700        1
701        1
702        -1
703        -1
704        -1
705        0
706        0
707        0
708        1
709        1
710        1
711        -1
712        -1
713        -1
714        0
715        0
716        0
717        1
718        1
719        1
720        -1
721        -1
722        -1
723        0
724        0
725        0
726        1
727        1
728        1

/
;
parameter veci6(direc)
/
1        -1
2        0
3        1
4        -1
5        0
6        1
7        -1
8        0
9        1
10        -1
11        0
12        1
13        -1
14        0
15        1
16        -1
17        0
18        1
19        -1
20        0
21        1
22        -1
23        0
24        1
25        -1
26        0
27        1
28        -1
29        0
30        1
31        -1
32        0
33        1
34        -1
35        0
36        1
37        -1
38        0
39        1
40        -1
41        0
42        1
43        -1
44        0
45        1
46        -1
47        0
48        1
49        -1
50        0
51        1
52        -1
53        0
54        1
55        -1
56        0
57        1
58        -1
59        0
60        1
61        -1
62        0
63        1
64        -1
65        0
66        1
67        -1
68        0
69        1
70        -1
71        0
72        1
73        -1
74        0
75        1
76        -1
77        0
78        1
79        -1
80        0
81        1
82        -1
83        0
84        1
85        -1
86        0
87        1
88        -1
89        0
90        1
91        -1
92        0
93        1
94        -1
95        0
96        1
97        -1
98        0
99        1
100        -1
101        0
102        1
103        -1
104        0
105        1
106        -1
107        0
108        1
109        -1
110        0
111        1
112        -1
113        0
114        1
115        -1
116        0
117        1
118        -1
119        0
120        1
121        -1
122        0
123        1
124        -1
125        0
126        1
127        -1
128        0
129        1
130        -1
131        0
132        1
133        -1
134        0
135        1
136        -1
137        0
138        1
139        -1
140        0
141        1
142        -1
143        0
144        1
145        -1
146        0
147        1
148        -1
149        0
150        1
151        -1
152        0
153        1
154        -1
155        0
156        1
157        -1
158        0
159        1
160        -1
161        0
162        1
163        -1
164        0
165        1
166        -1
167        0
168        1
169        -1
170        0
171        1
172        -1
173        0
174        1
175        -1
176        0
177        1
178        -1
179        0
180        1
181        -1
182        0
183        1
184        -1
185        0
186        1
187        -1
188        0
189        1
190        -1
191        0
192        1
193        -1
194        0
195        1
196        -1
197        0
198        1
199        -1
200        0
201        1
202        -1
203        0
204        1
205        -1
206        0
207        1
208        -1
209        0
210        1
211        -1
212        0
213        1
214        -1
215        0
216        1
217        -1
218        0
219        1
220        -1
221        0
222        1
223        -1
224        0
225        1
226        -1
227        0
228        1
229        -1
230        0
231        1
232        -1
233        0
234        1
235        -1
236        0
237        1
238        -1
239        0
240        1
241        -1
242        0
243        1
244        -1
245        0
246        1
247        -1
248        0
249        1
250        -1
251        0
252        1
253        -1
254        0
255        1
256        -1
257        0
258        1
259        -1
260        0
261        1
262        -1
263        0
264        1
265        -1
266        0
267        1
268        -1
269        0
270        1
271        -1
272        0
273        1
274        -1
275        0
276        1
277        -1
278        0
279        1
280        -1
281        0
282        1
283        -1
284        0
285        1
286        -1
287        0
288        1
289        -1
290        0
291        1
292        -1
293        0
294        1
295        -1
296        0
297        1
298        -1
299        0
300        1
301        -1
302        0
303        1
304        -1
305        0
306        1
307        -1
308        0
309        1
310        -1
311        0
312        1
313        -1
314        0
315        1
316        -1
317        0
318        1
319        -1
320        0
321        1
322        -1
323        0
324        1
325        -1
326        0
327        1
328        -1
329        0
330        1
331        -1
332        0
333        1
334        -1
335        0
336        1
337        -1
338        0
339        1
340        -1
341        0
342        1
343        -1
344        0
345        1
346        -1
347        0
348        1
349        -1
350        0
351        1
352        -1
353        0
354        1
355        -1
356        0
357        1
358        -1
359        0
360        1
361        -1
362        0
363        1
364        -1
365        1
366        -1
367        0
368        1
369        -1
370        0
371        1
372        -1
373        0
374        1
375        -1
376        0
377        1
378        -1
379        0
380        1
381        -1
382        0
383        1
384        -1
385        0
386        1
387        -1
388        0
389        1
390        -1
391        0
392        1
393        -1
394        0
395        1
396        -1
397        0
398        1
399        -1
400        0
401        1
402        -1
403        0
404        1
405        -1
406        0
407        1
408        -1
409        0
410        1
411        -1
412        0
413        1
414        -1
415        0
416        1
417        -1
418        0
419        1
420        -1
421        0
422        1
423        -1
424        0
425        1
426        -1
427        0
428        1
429        -1
430        0
431        1
432        -1
433        0
434        1
435        -1
436        0
437        1
438        -1
439        0
440        1
441        -1
442        0
443        1
444        -1
445        0
446        1
447        -1
448        0
449        1
450        -1
451        0
452        1
453        -1
454        0
455        1
456        -1
457        0
458        1
459        -1
460        0
461        1
462        -1
463        0
464        1
465        -1
466        0
467        1
468        -1
469        0
470        1
471        -1
472        0
473        1
474        -1
475        0
476        1
477        -1
478        0
479        1
480        -1
481        0
482        1
483        -1
484        0
485        1
486        -1
487        0
488        1
489        -1
490        0
491        1
492        -1
493        0
494        1
495        -1
496        0
497        1
498        -1
499        0
500        1
501        -1
502        0
503        1
504        -1
505        0
506        1
507        -1
508        0
509        1
510        -1
511        0
512        1
513        -1
514        0
515        1
516        -1
517        0
518        1
519        -1
520        0
521        1
522        -1
523        0
524        1
525        -1
526        0
527        1
528        -1
529        0
530        1
531        -1
532        0
533        1
534        -1
535        0
536        1
537        -1
538        0
539        1
540        -1
541        0
542        1
543        -1
544        0
545        1
546        -1
547        0
548        1
549        -1
550        0
551        1
552        -1
553        0
554        1
555        -1
556        0
557        1
558        -1
559        0
560        1
561        -1
562        0
563        1
564        -1
565        0
566        1
567        -1
568        0
569        1
570        -1
571        0
572        1
573        -1
574        0
575        1
576        -1
577        0
578        1
579        -1
580        0
581        1
582        -1
583        0
584        1
585        -1
586        0
587        1
588        -1
589        0
590        1
591        -1
592        0
593        1
594        -1
595        0
596        1
597        -1
598        0
599        1
600        -1
601        0
602        1
603        -1
604        0
605        1
606        -1
607        0
608        1
609        -1
610        0
611        1
612        -1
613        0
614        1
615        -1
616        0
617        1
618        -1
619        0
620        1
621        -1
622        0
623        1
624        -1
625        0
626        1
627        -1
628        0
629        1
630        -1
631        0
632        1
633        -1
634        0
635        1
636        -1
637        0
638        1
639        -1
640        0
641        1
642        -1
643        0
644        1
645        -1
646        0
647        1
648        -1
649        0
650        1
651        -1
652        0
653        1
654        -1
655        0
656        1
657        -1
658        0
659        1
660        -1
661        0
662        1
663        -1
664        0
665        1
666        -1
667        0
668        1
669        -1
670        0
671        1
672        -1
673        0
674        1
675        -1
676        0
677        1
678        -1
679        0
680        1
681        -1
682        0
683        1
684        -1
685        0
686        1
687        -1
688        0
689        1
690        -1
691        0
692        1
693        -1
694        0
695        1
696        -1
697        0
698        1
699        -1
700        0
701        1
702        -1
703        0
704        1
705        -1
706        0
707        1
708        -1
709        0
710        1
711        -1
712        0
713        1
714        -1
715        0
716        1
717        -1
718        0
719        1
720        -1
721        0
722        1
723        -1
724        0
725        1
726        -1
727        0
728        1
/
;


parameter objval(etapabup);
parameter xvalinit(etapabup,extvar) "punto a evaluar optimalidad integral para cada etapa de boil up";
parameter fvalueinit(etapabup) "valor de funcion objetivo en el punto a evaluar";
parameter fvaluedirec(etapabup,direc) "valor de funcion objetivo en cada direccion del vecindario";
parameter dif_f(etapabup,direc) "diferencia del valor de la funcion objetivo entre los vecinos y el optimo encontrado [$/year]";
parameter gval(extineq) "valor de g en x";
parameter mstat1(etapabup,direc) "status del solver";
parameter cumple(etapabup) "1 si cumple el criterio de optimalidad. 0 de lo contrario";
scalar CPUtime "Tiempo [s]";




loop(etapabup,

Nboil=Nboilloop(etapabup);
xvalinit(etapabup,"x1")=6;
xvalinit(etapabup,"x2")=2;
xvalinit(etapabup,"x3")=0;
xvalinit(etapabup,"x4")=0;
xvalinit(etapabup,"x5")=5;
xvalinit(etapabup,"x6")=17;
fvalueinit(etapabup)=22410.2062755944;

loop(direc,




         if(ord(direc) eq 1,
         gval("g1")=2-(xvalinit(etapabup,"x1"));
         gval("g2")=((xvalinit(etapabup,"x6"))-(xvalinit(etapabup,"x2")))-(xvalinit(etapabup,"x6"));
         gval("g3")=-(xvalinit(etapabup,"x3"));
         gval("g4")=-(xvalinit(etapabup,"x4"));
         gval("g5")=1-(xvalinit(etapabup,"x5"));
         gval("g6")=(xvalinit(etapabup,"x1"))+(xvalinit(etapabup,"x3"))+(xvalinit(etapabup,"x5"))+(xvalinit(etapabup,"x4"))-((xvalinit(etapabup,"x6"))-(xvalinit(etapabup,"x2")))+1;
         gval("g7")=(xvalinit(etapabup,"x6"))-21;
         else
         gval("g1")=2-(xvalinit(etapabup,"x1")+veci1(direc));
         gval("g2")=((xvalinit(etapabup,"x6")+veci6(direc))-(xvalinit(etapabup,"x2")+veci2(direc)))-(xvalinit(etapabup,"x6")+veci6(direc));
         gval("g3")=-(xvalinit(etapabup,"x3")+veci3(direc));
         gval("g4")=-(xvalinit(etapabup,"x4")+veci4(direc));
         gval("g5")=1-(xvalinit(etapabup,"x5")+veci5(direc));
         gval("g6")=(xvalinit(etapabup,"x1")+veci1(direc))+(xvalinit(etapabup,"x3")+veci3(direc))+(xvalinit(etapabup,"x5")+veci5(direc))+(xvalinit(etapabup,"x4")+veci4(direc))-((xvalinit(etapabup,"x6")+veci6(direc))-(xvalinit(etapabup,"x2")+veci2(direc)))+1;
         gval("g7")=(xvalinit(etapabup,"x6")+veci6(direc))-21;
         );



         if(sum(extineq$(gval(extineq) gt 0 ),gval(extineq)) gt 0,
         fvaluedirec(etapabup,direc)=1e+8;
         else
         yf(Net,F)=0;
         yf(Net,'1')$(ord(net) eq (xvalinit(etapabup,"x1")+veci1(direc)))=1;
         yf(Net,'2')$(ord(net) eq (xvalinit(etapabup,"x6")+veci6(direc))-(xvalinit(etapabup,"x2")+veci2(direc)))=1;
         yc(Net)=0;
         yc(net)$(ord(Net) eq (xvalinit(etapabup,"x1")+veci1(direc))+(xvalinit(etapabup,"x3")+veci3(direc)))=1;
         yc(net)$(ord(Net) eq ((xvalinit(etapabup,"x6")+veci6(direc))-(xvalinit(etapabup,"x2")+veci2(direc)))-(xvalinit(etapabup,"x4")+veci4(direc)))=1;
         yc(net)$(ord(Net) eq (xvalinit(etapabup,"x1")+veci1(direc))+(xvalinit(etapabup,"x3")+veci3(direc))+(xvalinit(etapabup,"x5")+veci5(direc)))=1;
         yr(Net)$((ord(Net) ne card(Net)) and (ord(Net) ne 1))=0;
         yr('2')=1;
         yb(Net)$((ord(Net) ne card(Net)) and (ord(Net) ne 1))=0;
         yb(Net)$(ord(Net) eq xvalinit(etapabup,"x6")+veci6(direc))=1;
         par(Net)$((ord(Net) ne card(Net)) and (ord(Net) ne 1))=(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yr(Net1)))+(yb(Net))-(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yb(Net1)));
         ye(Net)$((ord(Net) ne card(Net)) and (ord(Net) ne 1))=par(Net)*(1-yc(Net))*CASE+par(Net)*(1-CASE);

                  if(Nboilloop(etapabup) eq 9,
                  execute_loadpoint "ALG_8";
                  solve Rx_OCP_index1 using nlp minimizing zobj;
                  mstat1(etapabup,direc)=Rx_OCP_index1.modelstat;

                           if(mstat1(etapabup,direc) eq 3 or mstat1(etapabup,direc) eq 4 or mstat1(etapabup,direc) eq 5 or mstat1(etapabup,direc) eq 6 or mstat1(etapabup,direc) eq 11 or mstat1(etapabup,direc) eq 12 or mstat1(etapabup,direc) eq 13 or mstat1(etapabup,direc) eq 14 or mstat1(etapabup,direc) eq 18 or mstat1(etapabup,direc) eq 19,
                           fvaluedirec(etapabup,direc)=1e+8;
                           else
                           fvaluedirec(etapabup,direc)=zobj.l;
                           );


                  );


         );


dif_f(etapabup,direc)= fvaluedirec(etapabup,direc)-fvalueinit(etapabup);


);


if(sum(direc$(dif_f(etapabup,direc) lt 0),dif_f(etapabup,direc)) lt 0,
cumple(etapabup)=0;
else
cumple(etapabup)=1;
);

);
CPUtime=timeElapsed;

execute_unload "integral_optimalverification";
execute_unload "integral_optimalverification_INFO",cumple,CPUtime,dif_f;

