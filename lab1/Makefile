main.a: main.o matrix.o ode_solver.o opt_alg.o solution.o user_funs.o
	g++ main.o matrix.o ode_solver.o opt_alg.o solution.o user_funs.o -o main.a

main.o: main.cpp opt_alg.h matrix.h ode_solver.h solution.h user_funs.h
	g++ -c main.cpp -o main.o

matrix.o: matrix.cpp matrix.h
	g++ -c matrix.cpp -o matrix.o

ode_solver.o: ode_solver.cpp ode_solver.h
	g++ -c ode_solver.cpp -o ode_solver.o

opt_alg.o: opt_alg.cpp opt_alg.h
	g++ -c opt_alg.cpp -o opt_alg.o

solution.o: solution.cpp solution.h
	g++ -c solution.cpp -o solution.o

user_funs.o: user_funs.cpp user_funs.h
	g++ -c user_funs.cpp -o user_funs.o

run:
	./main.a

clean:
	rm -f *.o main.a