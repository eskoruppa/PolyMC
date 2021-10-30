#ifndef __DUMP_ENDBEAD_H_INCLUDED__
#define __DUMP_ENDBEAD_H_INCLUDED__
#include "Dump.h"

////////////////////////////////////////////////////////////////////////////////

#define DUMP_ENDBEAD_DEFAULT_STORE_NUM 1000


class Dump_Endbead;
class Dump_Endbead: public Dump
{
protected:
    int  store_counter = 0;
    int  store_num     = 0;
    arma::mat store_topol;

public:
    Dump_Endbead(Chain * ch, int N_dump, const std::string& filename, bool append=true);
    Dump_Endbead(Chain * ch, int N_dump, int N_to_file, const std::string& filename, bool append=true);
    ~Dump_Endbead();
    void init_store(int num);
    void write2file();

    void prod_dump();
    void final_dump();
};


class Dump_Endbead_Stats;
class Dump_Endbead_Stats: public Dump
{
protected:
    int    sum_dump_every;
    double sum;
    double sqsum;
    int    count_sum;

public:
    Dump_Endbead_Stats(Chain * ch, int N_dump, const std::string& filename, bool append=true);
    ~Dump_Endbead_Stats();
    void prod_dump();
    void final_dump();
};


class Dump_Endbead_Ceff;
class Dump_Endbead_Ceff: public Dump
{
protected:
    int    sum_dump_every;
    double sum;
    double sqsum;
    int    count_sum;

public:
    Dump_Endbead_Ceff(Chain * ch, int N_dump, const std::string& filename, bool append=true);
    ~Dump_Endbead_Ceff();
    void prod_dump();
    void final_dump();
};


#endif
