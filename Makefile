OPT=-g -Wall

# Biblioteca

bisnwt.o : bisnwt.c
	gcc $(OPT) -c bisnwt.c

# Utilitat

kplt2nu : kplt2nu.c bisnwt.o
	gcc $(OPT) -o kplt2nu kplt2nu.c bisnwt.o -lm

# Test amb exp(x)-2
prova_expm2 : prova_expm2.c bisnwt.o
	gcc $(OPT) -o prova_expm2 prova_expm2.c bisnwt.o -lm

# Neteja

clean :
	rm -f bisnwt.o

realclean : clean
	rm -f kplt2nu
