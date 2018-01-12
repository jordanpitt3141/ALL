CC = gcc
CFLAGS = -I/usr/include/python2.7 -Wall -fPIC -O3


SWWE_wrap.c:
	swig -python SWWE.i

SWWE: SWWE_wrap.o
	ld -shared SWWE.o SWWE_wrap.o -o _SWWE.so 

%.o : %.c
	$(CC) -c SWWE.c $^ $(CFLAGS)


clean:
	rm -f ./*.o *~ ./SWWE.py ./SWWE_wrap.c ./*.pyc ./*.so

