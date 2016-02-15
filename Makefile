sim:		main.o
		g++ -o sim main.o

main.o:		main.cpp main.h GA.h AS.h S.h
		g++ -Wall -c main.cpp


