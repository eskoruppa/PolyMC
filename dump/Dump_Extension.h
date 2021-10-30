#ifndef __DUMP_EXTENSION_H__
#define __DUMP_EXTENSION_H__
#include "Dump.h"

#define DEXT_STORE_EXTS 1000

////////////////////////////////////////////////////////////////////////////////

class Dump_Extension;

class Dump_Extension: public Dump
{
protected:
    int iID,fID;
    arma::colvec store_exts;
    unsigned store_counter;

public:
    Dump_Extension(Chain * ch, int N_dump, const std::string& filename, bool append=true);
    Dump_Extension(Chain * ch, int N_dump, const std::string& filename, int start_id, int end_id, bool append=true);
    ~Dump_Extension();

    void init_store();
    void write2file();

    void prod_dump();
    void final_dump();
};

#endif
