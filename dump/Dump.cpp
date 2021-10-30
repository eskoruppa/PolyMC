#include "Dump.h"

Dump::Dump(Chain * ch, int N_dump, const std::string& filename,bool append) :
// Set State Variables
chain(ch),
BPS(*ch->get_BPS()),
triads(ch->get_triads()),
pos(ch->get_bp_pos()),
N_step(N_dump),
step(0),
fn(filename),
num_bp(ch->get_num_bp()),
num_bps(ch->get_num_bps()),
app(append),
disc_len(ch->get_disc_len()) {}

Dump::~Dump() {}

void Dump::dump() {
    step++;
    if (step%N_step==0) {
        prod_dump();
    }
}

std::string Dump::get_filename() {
    return fn;
}

//void Dump::init() {}
void Dump::prod_dump() {}
void Dump::final_dump() {}

bool Dump::fexists(const std::string& filename)
{
  std::ifstream ifile(filename);
  return (bool)ifile;
}
