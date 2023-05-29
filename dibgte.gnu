# Definicions per dibuixar ground tracks
# Adaptat per orbites Molniya (mnpr02)
#
# CAL CARREGAR: ctants.gnu
#
# CAL DEFINIR:
# - ipert : 1 per tenir en compte la pert. del J2, 0 per no tenir-la en compte
# - lonesp : longitud del meridia del lloc que espia
# - dlapo : correcció per baixar la latitud de l'apogeu
# - xc : excentricitat de l'orbita
# - T : periode de l'orbita
# - i : inclinacio de l'orbita
# - tf : durada de la simulacio (segons)
# - nt : nombre de punts de la simulacio
# - ind : punt inicial de l'animacio (es comen,ca amb ind=0)

# Equacio de l'orbita
r(nu)=a*(1-xc**2)/(1+xc*cos(nu))

# Terme J2 del potencial gravitatori terrestre
J2=ipert*0.0010826267

# Moviment mig
n=2*pi/T

# Semieix a partir de la 3a llei de Kepler
a=(mutr/n**2)**(1./3)

# Parametre a partir de semieix i excentricitat
p=a*(1-xc**2)

# Retrogradacio de l'argument del node asc per segon (Vallado 3rd ed p. 645)
DOmg=-3*n*(Rtr**2)*J2/(2*(p**2))*cos(i)

# Retrogradacio de l'argument del peri per segon (Vallado 3rd ed p. 647)
Domg=3*n*(Rtr**2)*J2/(4*(p**2))*(4-5*sin(i)**2)

# Argument del "node asc" com a funcio del temps
Omg(t)=DOmg*t
# En fer això, estic imposant que l'eix x apunti cap al node ascendent,
# en comptes del vernal equinox.

# Argument del peri com a funcio del temps
omg(t)=-pi/2-dlapo+Domg*t
# Això és perquè, al pla orbital, l'apogeu de l'el·lipse estigui a la
# part positiva de l'eix y. Després de totes les rotacions, l'eix y
# apuntarà al meridià del lloc que espia.

# Posicio a sobre de l'orbita
xnu(nu)=r(nu)*cos(nu)
ynu(nu)=r(nu)*sin(nu)
znu(nu)=0

# Rotacio al pla orbital donada per argument del peri
xomg(nu,t)=cos(omg(t))*xnu(nu)-sin(omg(t))*ynu(nu)
yomg(nu,t)=sin(omg(t))*xnu(nu)+cos(omg(t))*ynu(nu)
zomg(nu,t)=znu(nu)

# Rotacio del pla orbital donada per la inclinacio
xi(nu,t)=xomg(nu,t)
yi(nu,t)=cos(i)*yomg(nu,t)-sin(i)*zomg(nu,t)
zi(nu,t)=sin(i)*yomg(nu,t)+cos(i)*zomg(nu,t)

# Rotacio del pla orbital donada per l'argument del node
xOmg(nu,t)=cos(Omg(t))*xi(nu,t)-sin(Omg(t))*yi(nu,t)
yOmg(nu,t)=sin(Omg(t))*xi(nu,t)+cos(Omg(t))*yi(nu,t)
zOmg(nu,t)=zi(nu,t)

# "Ascensio recta" del satel.lit (xyz amb y al meridia del lloc 
# que espia a temps zero [nomes canvia si ipert==1])
alf(nu,t)=atan2(yOmg(nu,t),xOmg(nu,t))

# Declinacio del satel.lit
dlt(nu,t)=asin(zOmg(nu,t)/r(nu))

# "Temps sideri" (xyz com abans) del lloc que espia
the(t)=pi/2+omgtr*t

# Determinacio principal d'un angle (a [-pi:pi])
detpr(th)=th-floor((th+pi)/(2*pi))*(2*pi)

# Latitud i longitud en graus
latd(t,nu)=dlt(nu,t)/pi*180
lond(t,nu)=detpr( lonesp+alf(nu,t)-the(t) )/pi*180

# Parametres de dibuix
set xrange [-180:180]
set yrange [-90:90]
set size ratio -1
#unset key

# Comanda per generar t M E
cmdpl=sprintf( \
   "kplt2nu %.16G=e %.16G=T %.16G=M0 %.16G=tf %d=nt > gtrk.txt",\
   xc,T,pi,tf,nt)
set zeroaxis

system cmdpl

plot 'wm.txt' w l,'gtrk.txt' u (lond($1,$3)):(latd($1,$3)) w p lt 3 pt 2
