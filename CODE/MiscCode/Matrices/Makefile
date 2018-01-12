CC = gcc
CFLAGS = -I/usr/include/python2.7 -Wall -fPIC -O3


matrixNR_wrap.c:
	swig -python matrixNR.i

matrixNR: matrixNR_wrap.o
	ld -shared matrixNR.o matrixNR_wrap.o -o _matrixNR.so 

%.o : %.c
	$(CC) -c matrixNR.c $^ $(CFLAGS)


clean:
	rm -f ./*.o *~ ./matrixNR.py ./matrixNR_wrap.c ./*.pyc ./*.so

