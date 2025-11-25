CXX = g++
HOME= /usr/local/include
LIB_HOME = ../../flowstar/flowstar-toolbox
LIBS = -lflowstar -lmpfr -lgmp -lgsl -lgslcblas -lm -lglpk
CFLAGS = -I . -I $(HOME) -g -O3 -std=c++11
LINK_FLAGS = -g -L$(LIB_HOME) -L/usr/local/lib
OUTPUT_DIR = outputs

all: $(OUTPUT_DIR) d_inv_pendulum

$(OUTPUT_DIR):
	mkdir -p $(OUTPUT_DIR)

d_inv_pendulum: d_inv_pendulum.o
	g++ -O3 -w $(LINK_FLAGS) -o $@ $^ $(LIBS)

%.0: %.cc
	$(CXX) -O3 -c $(CFLAGS) -o $@ $<

%.o: %.cpp
	$(CXX) -O3 -c $(CFLAGS) -o $@ $<

%.o: %.c
	$(CXX) -O3 -c $(CFLAGS) -o $@ $<

clean:
	rm -f *.o d_inv_pendulum
	rm -f $(OUTPUT_DIR)/*.eps $(OUTPUT_DIR)/*.plt