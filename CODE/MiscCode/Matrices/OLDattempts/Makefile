CC = gcc
CFLAGS = -I/usr/include/python2.7 -Wall -fPIC -O3


matrix1_wrap.c:
	swig -python matrix1.i

matrix1: matrix1_wrap.o
	ld -shared matrix1.o matrix1_wrap.o -o _matrix1.so 

%.o : %.c
	$(CC) -c matrix1.c $^ $(CFLAGS)


clean:
	rm -f ./*.o *~ ./matrix1.py ./matrix1_wrap.c ./*.pyc ./*.so

