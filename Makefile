CC=g++
CFLAGS= -Wall -Ofast -std=c++17  -flto -pipe -funit-at-a-time -fopenmp -lz -Isparsepp -flto
LDFLAGS=-flto -lpthread -fopenmp -lz  -Isparsepp  -flto
LIBS=utils.h Bloom.h  ExponentialBloom.h bcardi.h best.h bestpart.h
EXEC=bcardi best


best: main.o best.o  Bloom.o ExponentialBloom.o utils.o bestpart.o
	$(CC) -o $@ $^ $(LDFLAGS)

# bcardi: main.o bcardi.o  Bloom.o ExponentialBloom.o utils.o
# 	$(CC) -o $@ $^ $(LDFLAGS)

main.o: main.cpp $(LIBS)
	$(CC) -o $@ -c $< $(CFLAGS)

best.o: best.cpp $(LIBS)
	$(CC) -o $@ -c $< $(CFLAGS)

bestpart.o: BestPart.cpp $(LIBS)
	$(CC) -o $@ -c $< $(CFLAGS)

bcardi.o: bcardi.cpp $(LIBS)
	$(CC) -o $@ -c $< $(CFLAGS)

Bloom.o: Bloom.cpp $(LIBS)
	$(CC) -o $@ -c $< $(CFLAGS) 

ExponentialBloom.o: ExponentialBloom.cpp $(LIBS)
	$(CC) -o $@ -c $< $(CFLAGS) 

utils.o: utils.cpp $(LIBS)
	$(CC) -o $@ -c $< $(CFLAGS) 

clean:
	rm -rf *.o
	rm -rf $(EXEC)


rebuild: clean $(EXEC)
