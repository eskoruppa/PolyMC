#ifndef __DUMP_HatCurve_H__
#define __DUMP_HatCurve_H__
#include "Dump.h"

////////////////////////////////////////////////////////////////////////////////
/*
    Dump the extension of Writhe of a Tweezer Simulations
    Thus far only works in the constant linking number ensemble
*/

class Dump_HatCurve;
class Dump_HatCurve: public Dump
{
protected:
    int N_dump_to_file;
    double z,zsq;
    double x,xsq;
    double y,ysq;
    double Wr,Wrsq;
    double Tw,Twsq;
    int counter;

public:
    Dump_HatCurve(Chain * ch, int N_dump, const std::string& filename, int N_dump_to_file, bool append=true);
    ~Dump_HatCurve();

    void prod_dump();
    void final_dump();
};


class Dump_HatCurveStatistics;
class Dump_HatCurveStatistics: public Dump
{
public:
    Dump_HatCurveStatistics(Chain * ch, int N_dump, std::string filename, bool append=true);
    ~Dump_HatCurveStatistics();

    void prod_dump();
    void final_dump();
};



#endif
