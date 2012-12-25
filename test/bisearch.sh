g++ -O0 c/bisearch.cpp -o bisearch
time ./bisearch 1> /dev/null
g++ -O1 c/bisearch.cpp -o bisearch
time ./bisearch 1> /dev/null
g++ -O2 c/bisearch.cpp -o bisearch
time ./bisearch 1> /dev/null
g++ -O3 c/bisearch.cpp -o bisearch
time ./bisearch 1> out_bisearch_time.txt
