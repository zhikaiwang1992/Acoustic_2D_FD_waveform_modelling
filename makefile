#! /bin/sh

# Make wave modeling programs
INC=/home/wzk/software/su/include
LIK=/home/wzk/software/su/lib
LIB=-lsu -lpar -lcwp -lm

main : main.c FDTD_2D_for.c
	gcc -O2 main.c -o main -I$(INC) -L$(LIK)  $(LIB)

