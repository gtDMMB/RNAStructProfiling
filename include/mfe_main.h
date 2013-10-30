#ifndef _MFE_MAIN_H_
#define _MFE_MAIN_H_

int mfe_main(int argc, char** argv);
//double calculate_mfe(int argc, char** argv);
double calculate_mfe(std::string seq);
void init_fold(const char* seq);
void parse_mfe_options(int argc, char** argv);
void free_fold(int len);

#endif
