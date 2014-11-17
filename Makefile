sim:		source_code/main.o
		g++ -o sim souurce_code/main.o

main.o:		source_code/main.cpp source_code/main.h source_code/GA.h source_code/AS.h source_code/S.h
		 g++ -Wall -c source_code/main.cpp


