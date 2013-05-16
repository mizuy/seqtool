#ifndef _PRIMER_H_
#define _PRIMER_H_

#include <string>

void print_primer_pair(const std::string& fw, const std::string& rv, bool end_annealing);
void print_primer_bar(const std::string& fw, const std::string& rv, int index);

std::tuple<int,int> annealing_score(const std::string& fw, const std::string& rv);
std::tuple<int,int> end_annealing_score(const std::string& fw, const std::string& rv);

#endif
