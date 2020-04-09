LIB_DIR=/usr/local/lib
NTL_INC=/usr/local/include/NTL
GMPL_INC=/usr/local/include

all: clean
	g++ -Wall -Wpedantic -g -O5 -I${NTL_INC} -I${GMPL_INC} ky_v2.cpp -o ky_v2 -std=c++11 -L${LIB_DIR} -lntl -lgmp -lm -lpthread

ky_v2: clean
	g++ -Wall -Wpedantic -g -O5 -I${NTL_INC} -I${GMPL_INC} ky_v2.cpp -o ky_v2 -std=c++11 -L${LIB_DIR} -lntl -lgmp -lm -lpthread

test: clean
	g++ -Wall -Wpedantic -g -O5 test.cpp -o test -std=c++11 -L${LIB_DIR} 	

debug: clean
	g++ -Wall -Wpedantic -g -I${NTL_INC} -I${GMPL_INC} ky_v2.cpp -o ky_v2 -std=c++11 -L${LIB_DIR} -lntl -lgmp -lm -lpthread
	gdb hanmat_mt

clean:
	rm -rf ky_v2 test
