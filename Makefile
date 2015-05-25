sim:		src/main.o
		g++ -o sim src/main.o

main.o:		src/main.cpp src/main.h src/GA.h src/AS.h src/S.h
		g++ -Wall -c src/main.cpp


