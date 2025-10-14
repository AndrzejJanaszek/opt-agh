# Kompilator i flagi
CXX = g++
CXXFLAGS = #-Wall -Wextra -O2

# Pliki
TARGET = main.a
OBJS = main.o matrix.o ode_solver.o opt_alg.o solution.o user_funs.o

# Domyślny cel
all: $(TARGET)

# Budowanie programu
$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $(TARGET)

# Reguły dla plików .o
main.o: main.cpp opt_alg.h matrix.h ode_solver.h solution.h user_funs.h
	$(CXX) $(CXXFLAGS) -c main.cpp -o main.o

matrix.o: matrix.cpp matrix.h
	$(CXX) $(CXXFLAGS) -c matrix.cpp -o matrix.o

ode_solver.o: ode_solver.cpp ode_solver.h
	$(CXX) $(CXXFLAGS) -c ode_solver.cpp -o ode_solver.o

opt_alg.o: opt_alg.cpp opt_alg.h
	$(CXX) $(CXXFLAGS) -c opt_alg.cpp -o opt_alg.o

solution.o: solution.cpp solution.h
	$(CXX) $(CXXFLAGS) -c solution.cpp -o solution.o

user_funs.o: user_funs.cpp user_funs.h
	$(CXX) $(CXXFLAGS) -c user_funs.cpp -o user_funs.o

# Uruchamianie programu
run: $(TARGET)
	./$(TARGET)

# Czyszczenie
clean:
	rm -f *.o $(TARGET)