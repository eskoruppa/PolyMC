#ifndef __INIT_DUMP__
#define __INIT_DUMP__

#include <vector>

// include input functions
#include "../Input/Argparse.h"
#include "../Input/InputRead.h"
#include "../Input/InputChoice.h"
#include "../Input/GenInputFile.h"

// include Dumps
#include "Dump.h"
#include "Dump_Conf.h"
#include "Dump_DistanceMap.h"
#include "Dump_Endbead.h"
#include "Dump_EndToEndDistance.h"
#include "Dump_Energy.h"
#include "Dump_Extension.h"
#include "Dump_ForceExtension.h"
#include "Dump_Full_Coup_Matrix.h"
#include "Dump_HatCurve.h"
#include "Dump_Helicity.h"
#include "Dump_Kinkstats.h"
#include "Dump_Linkingnumber.h"
#include "Dump_MSD.h"
#include "Dump_PersistenceLength.h"
#include "Dump_PlasmidEndpoints.h"
#include "Dump_PlecStats.h"
#include "Dump_ProximityPlec.h"
#include "Dump_Restart.h"
#include "Dump_SolenoidCorrelationFunction.h"
#include "Dump_State.h"
#include "Dump_Stiff.h"
#include "Dump_SupercoilingPhase.h"
#include "Dump_TangentCorrelation.h"
#include "Dump_Thetas.h"
#include "Dump_Torque.h"
#include "Dump_TorsionalStiffness.h"
#include "Dump_Twist.h"
#include "Dump_TorqueTrap.h"
#include "Dump_WritheMap.h"
#include "Dump_xyz.h"
#include "Dump_xyz_extreme_ext.h"


std::vector<Dump*> init_dump_cmdargs(const std::vector<std::string> & argv, GenInputFile * geninfile, Chain * _chain, std::string mode, double EV_rad, std::string inputfile="");

#endif
