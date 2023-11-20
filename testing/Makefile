CC = gcc

default: libexamples.a

libexamples.a: examples.o
	ar rcs $@ $^
	
examples.o: examples.c examples.h
	$(CC) -c $<

clean:
	rm *.o *.a
	

default: pyexamples

pyexamples: setup.py pyexamples.pyx libexamples.a
	python3 setup.py build_ext --inplace && rm -f pyexamples.c && rm -Rf build

libexamples.a:
	make -C libexamples.a

clean:
	rm *.so