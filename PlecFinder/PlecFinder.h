#ifndef __INCLUDED_PLECFINDER__
#define __INCLUDED_PLECFINDER__

#include "../Chain.h"
#include "../SO3Methods.h"
#include "../ExtraFuncs.h"
#include "../ExtraArma.h"
#include "BranchBoxes.h"

#define _USE_MATH_DEFINES

#include <armadillo> // armadillo
#include <iostream>  // input output
#include <string>    // string
#include <algorithm>
#include <stdexcept> // exceptions
#include <fstream>   // read file
#include <cmath>     // math lib
#include <vector>    // vector
#include <cassert>   // assert

//#define PLECFINDER_DEBUG


class PlecFinder;

class PlecFinder {
protected:

    double   WRTHRESHOLD                 = 0.25;
    double   DWPEAKTHRESHOLD             = 4;
    double   BRANCH_THRES_MEANFAC        = 3;
    double   BRANCH_THRES_MAXFAC         = 0.20;
    double   BRANCH_CONNECT_THRESHOLD    = 10;
    double   BRANCH_MINFRAC_CONTAINED_WR = 0.75;
    int      BRANCH_MAX_ITERATIONS       = 1;
    double   BRANCH_ITERATION_FAC        = 0.6667;
    double   BRANCH_ENDS_AVG_DIST        = 20;
    double   ENTRY_REFINEMENT_MAX_REMOVE = 0.02;


    Chain * chain;
    bool   closed;
    double disc_len;
    int    num_bp;
    int    num_bps;

    int sgnWM;
    arma::mat WM;
    arma::mat pWM;
    arma::mat diag;

    std::vector<BranchBox*> bboxes;
    std::vector<int> entry_bboxes_ids;


public:

    PlecFinder(Chain * ch);
    PlecFinder(Chain * ch,const arma::mat & WM);
    ~PlecFinder();

    std::vector<double> getPlecStats();


private:
    void calBBoxes();

    arma::mat calPosWM();
    arma::mat diagTransform(const arma::mat & WM);
    std::vector<int> phi(int i, int j);
    std::vector<int> invphi(int k, int l);

    std::vector<int> phi(std::vector<int> P);
    std::vector<int> invphi(std::vector<int> Q);

    std::vector<std::vector<int>> findPeaks(const arma::mat & diag);
    bool local_minimum(const arma::colvec & DW,int id);

    std::vector<BranchBox*> findBranches(const arma::mat & diag,const std::vector<std::vector<int>> & intervals);
    std::vector<std::vector<std::vector<int>>> BranchesInInterval(const arma::mat & diag, const std::vector<int> & interval);

    void sortBboxes(std::vector<BranchBox*> * bboxes);
    void purgeInactiveBboxes(std::vector<BranchBox*> * bboxes);
    void purgeBboxesAtID(std::vector<BranchBox*> * bboxes,int id);

    // Remove BranchBox Overlap
    void removeBboxOverlap(std::vector<BranchBox*> * bboxes);
    std::vector<std::vector<int>> __balanceXOverlap(const std::vector<int> & xlimA,const std::vector<int> & ylimA,const std::vector<int> & xlimB,const std::vector<int> & ylimB);
    std::vector<std::vector<int>> __balanceXOverlapL2R(const std::vector<int> & xlim1,const std::vector<int> & ylim1,const std::vector<int> & xlim2,const std::vector<int> & ylim2, double * totalwr);

    // Resolve Inconsistencies between bboxes
    void simpleConflictRemoval(std::vector<BranchBox*> * bboxes);
    bool isDownstreamBranch(BranchBox * refbbox, BranchBox * checkbbox);

    // Find Entry BBoxes and refine the entry limits
    std::vector<int> findEntryBboxes(const std::vector<BranchBox*> & bboxes);
    void refineAllEntryBBoxes(const std::vector<BranchBox*> & bboxes,const std::vector<int> & entry_bboxes_ids, double max_remove);
    void refineEntryBBox(BranchBox * bbox,double max_remove);


};

#endif
