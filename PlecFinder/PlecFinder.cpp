#include "PlecFinder.h"

PlecFinder::PlecFinder(Chain * ch) :
chain(ch),
closed(ch->topology_closed()),
disc_len(ch->get_disc_len()),
num_bp(ch->get_num_bp()),
num_bps(ch->get_num_bps())
{
    chain->langowski_writhe_elements(&WM,closed);
    calBBoxes();
}

PlecFinder::PlecFinder(Chain * ch,const arma::mat & WM) :
chain(ch),
closed(ch->topology_closed()),
disc_len(chain->get_disc_len()),
WM(WM)
{
    calBBoxes();
}

PlecFinder::~PlecFinder() {
}


void PlecFinder::calBBoxes() {
    pWM  = calPosWM();
//    WM   = pWM;
    diag = diagTransform(pWM);

    std::vector<std::vector<int>> intervals = findPeaks(diag);
    bboxes = findBranches(diag,intervals);
    removeBboxOverlap(&bboxes);
    simpleConflictRemoval(&bboxes);
    entry_bboxes_ids = findEntryBboxes(bboxes);
    refineAllEntryBBoxes(bboxes,entry_bboxes_ids,ENTRY_REFINEMENT_MAX_REMOVE);


    #ifdef PLECFINDER_DEBUG
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    // DEBUGING
    std::ofstream out;
    out.open("dump/plecout", std::ofstream::out | std::ofstream::trunc);
    out << " Final " << std::endl;
    for (int i=0;i<bboxes.size();i++) {
        out << bboxes[i]->xlim[0] << " " << bboxes[i]->xlim[1] << " " << bboxes[i]->ylim[0] << " " << bboxes[i]->ylim[1] << " " << bboxes[i]->get_writhe() << std::endl;
    }
    out.close();

    std::cout << "##############" << std::endl;
    std::cout << "Remove Overlap" << std::endl;
    for (int i=0;i<bboxes.size();i++) {
        std::cout << bboxes[i]->xlim[0] << " " << bboxes[i]->xlim[1] << " " << bboxes[i]->ylim[0] << " " << bboxes[i]->ylim[1] << " " << bboxes[i]->get_writhe() << std::endl;
    }
    std::cout << "##############" << std::endl;
    std::cout << "##############" << std::endl;
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    #endif
}


arma::mat PlecFinder::calPosWM() {
    if (arma::accu(WM) < 0) {
        WM    = -WM;
        sgnWM = -1;
    }
    else {
        sgnWM = 1;
    }
    pWM = WM;
    for (unsigned i=0;i<pWM.n_rows;i++) {
        for (unsigned j=0;j<pWM.n_cols;j++) {
            if (pWM(i,j) < 0) {
                pWM(i,j) = 0;
            }
        }
    }
    return pWM;
}

arma::mat PlecFinder::diagTransform(const arma::mat & WM) {
    unsigned N = WM.n_cols;
    arma::mat diag = arma::zeros(N,int(std::ceil(N*0.5)));
    for (unsigned i=0;i<N;i++) {
        for (unsigned j=i+1;j<N;j++) {
            diag( int(std::floor((i+j)/2)) , int(std::floor((j-i)/2)) ) += 2*WM(i,j);
        }
    }
    return diag;
}

std::vector<int> PlecFinder::phi(int i, int j) {
    return {int(std::floor((i+j)/2)), int(std::floor((j-i)/2))};
}

std::vector<int> PlecFinder::phi(std::vector<int> P) {
    return phi(P[0],P[1]);
}

std::vector<int> PlecFinder::invphi(int k, int l) {
    return {int(k-l),int(k+l)};
}

std::vector<int> PlecFinder::invphi(std::vector<int> Q) {
    return invphi(Q[0], Q[1]);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////


std::vector<std::vector<int>> PlecFinder::findPeaks(const arma::mat & diag) {
/*
    This function identifies the positions of the branches along the diagonal of the writhe map.
    It identifies peaks in the y-integrated phi-map.
*/

    // Integrate the y axis of the diag (phi) map
    arma::colvec DW = arma::sum(diag,1);
    std::vector<std::vector<int>> intervals;

    /*
        Gauge the Width of the peaks for box averaging
    */
    std::vector<std::vector<int>> gauge_peaks;
    double mDW     = arma::mean(DW);
    double thres   = DWPEAKTHRESHOLD*mDW;
    bool   inpeak  = false;
    double psum;
    int    pfrom = 0;
    for (int i=0;i<DW.n_elem;i++) {

        if (DW(i) > thres) {
            if (!inpeak) {
                inpeak=true;
                pfrom = i;
            }
        }
        else {
            if (inpeak) {
                inpeak=false;
                psum = arma::sum(DW(arma::span(pfrom, i-1)));
                if (psum>WRTHRESHOLD) {
                    gauge_peaks.push_back({pfrom,i});
                }
            }
        }
    }

    // Assume that the chain contains no plectonemes if no gauge peaks are found
    if (gauge_peaks.size() == 0) return intervals;

    /*
        Smooth Diagonal Writhe
    */
    arma::colvec smoothedDW = arma::zeros(DW.n_elem);
    double mpeakwidth = 0;
    for (int i=0;i<gauge_peaks.size();i++) {
        mpeakwidth += gauge_peaks[i][1]-gauge_peaks[i][0];
    }
    mpeakwidth  /= gauge_peaks.size();
    int hwbox = int(std::ceil(mpeakwidth*0.5));

    int idfrom,idto;
    for (int i=0;i<DW.n_elem;i++) {
        idfrom = larger(i-hwbox,0);
        idto   = smaller(i+hwbox+1,DW.size());
        smoothedDW(i) = arma::sum(DW(arma::span(idfrom,idto-1)));
    }

    /*
        detect branches
    */

    inpeak = false;
    pfrom  = 0;
    for (int i=0;i<smoothedDW.n_elem;i++) {
        if (inpeak) {
            if (smoothedDW(i) < WRTHRESHOLD) {
                inpeak = false;
                if (arma::sum(DW(arma::span(pfrom,i-1)))>WRTHRESHOLD) {
                    intervals.push_back({pfrom,i});
                }
            }
            else if (local_minimum(smoothedDW,i)) {
                if (arma::sum(DW(arma::span(pfrom,i-1)))>WRTHRESHOLD) {
                    intervals.push_back({pfrom,i});
                }
                pfrom = i;
            }
        }
        else if (smoothedDW(i) > WRTHRESHOLD) {
            inpeak = true;
            pfrom  = i;
        }
    }
    if (inpeak) {
        if (arma::sum(DW(arma::span(pfrom,smoothedDW.n_elem)))>WRTHRESHOLD) {
            intervals.push_back({pfrom,smoothedDW.n_elem+1});
        }
    }
    return intervals;
}

bool PlecFinder::local_minimum(const arma::colvec & DW,int id) {
    if (id > 0 && id < DW.n_elem-1) {
        if (DW(id)<DW(id-1) && DW(id)<DW(id+1)) {
            return true;
        }
    }
    return false;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////


std::vector<BranchBox*> PlecFinder::findBranches(const arma::mat & diag,const std::vector<std::vector<int>> & intervals) {

    std::vector<BranchBox*> bboxes;
    std::vector<float>      midx;
    std::vector<int>        xlim,ylim;
    for (unsigned i=0;i<intervals.size();i++) {
        std::vector<int> interval = intervals[i];
        std::vector<std::vector<std::vector<int>>> termini = BranchesInInterval(diag,interval);

        for (int j=0;j<termini.size();j++) {
            xlim = termini[j][0];
            ylim = termini[j][1];

            BranchBox* bbox = new BranchBox(xlim,ylim,&WM);
            bboxes.push_back(bbox);
        }
    }
//    sortBboxes(&bboxes);
    return bboxes;
}

std::vector<std::vector<std::vector<int>>> PlecFinder::BranchesInInterval(const arma::mat & diag, const std::vector<int> & interval) {

    arma::mat diagint   = diag(arma::span(interval[0],interval[1]-1),arma::span(0,diag.n_cols-1));
    arma::colvec horsum = arma::sum(diagint,0);

    double meanmap   = arma::mean(horsum) * BRANCH_THRES_MEANFAC;
    double maxmap    = arma::max (horsum) * BRANCH_THRES_MAXFAC;
    double basethres = smaller(meanmap,maxmap);

    double total_wr         = arma::sum(horsum);
    double accounted_writhe = 0;
    int    iteration        = 0;

    double thres,segwr;
    bool   inpeak;
    int    pfrom;

    std::vector<std::vector<int>> segs;
    std::vector<int> seg;

    while (iteration < BRANCH_MAX_ITERATIONS) {

        /*
            Find directly connected segments
        */

        thres = basethres*std::pow(BRANCH_ITERATION_FAC,iteration);
        iteration++;

        std::vector<std::vector<int>> pieces;
        inpeak = false;
        pfrom  = 0;
        for (int i=0;i<horsum.n_elem;i++) {
            if (horsum(i) > thres) {
                if (!inpeak) {
                    inpeak=true;
                    pfrom = i;
                }
            }
            else if (inpeak) {
                inpeak = false;
                pieces.push_back({pfrom,i});
            }
        }

        /*
            Connect Segments
        */
        std::vector<std::vector<int>> consegs;
        consegs.push_back(pieces[0]);
        for (int i=1;i<pieces.size();i++) {
            if (pieces[i][0]-consegs[consegs.size()-1][1] < BRANCH_CONNECT_THRESHOLD) {
                consegs[consegs.size()-1][1] = pieces[i][1];
            }
            else {
                consegs.push_back(pieces[i]);
            }
        }

        segs.clear();
        segs.shrink_to_fit();
        accounted_writhe = 0;
        for (int i=0;i<consegs.size();i++) {
            seg = consegs[i];
            segwr = arma::sum(horsum(arma::span(seg[0],seg[1]-1)));
            if (segwr > WRTHRESHOLD) {
                segs.push_back(seg);
                accounted_writhe += segwr;
            }
        }
        if (accounted_writhe/total_wr >= BRANCH_MINFRAC_CONTAINED_WR) {
            break;
        }
    }

    /*
        Define start and end point
    */
    std::vector<std::vector<std::vector<int>>> termini;
    std::vector<int> P1,P2,Pn;
    arma::colvec tempvec;
    int x_from,x_to,y_from,y_to;
    int add = 2;

    for (unsigned i=0;i<segs.size();i++) {
        seg = segs[i];

        tempvec = diagint(arma::span(0,diagint.n_rows-1),seg[0]);
        P1 = invphi(interval[0] + arma_argmax(tempvec) , seg[0]);

        tempvec = diagint(arma::span(0,diagint.n_rows-1),seg[1]);
        P2 = invphi(interval[0] + arma_argmax(tempvec) , seg[1]);

        for (int coord=seg[0];coord<=seg[1];coord++) {
            tempvec = diagint(arma::span(0,diagint.n_rows-1),coord);
            Pn = invphi(interval[0] + arma_argmax(tempvec) , coord);
            if (Pn[0] > P1[0]) P1[0] = Pn[0];
            if (Pn[0] < P2[0]) P2[0] = Pn[0];

            if (Pn[1] < P1[1]) P1[1] = Pn[1];
            if (Pn[1] > P2[1]) P2[1] = Pn[1];
        }
//        P1 = invphi(phi(P1));
//        P2 = invphi(phi(P2));

        x_to   = P1[0]+add;
        y_from = P1[1]-add;

        x_from = P2[0]-add;
        y_to   = P2[1]+add;

        if (x_from < 0) x_from = 0;
        if (y_from < 0) y_from = 0;
        if (x_to > WM.n_cols) x_to = WM.n_cols;
        if (y_to > WM.n_cols) y_to = WM.n_cols;

        termini.push_back({{x_from,x_to},{y_from,y_to}});
    }
    return termini;
}


/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

void PlecFinder::sortBboxes(std::vector<BranchBox*> * bboxes) {
    /*
        Sort bboxes by order of xlim[0]
    */
    if (bboxes->size()==0) return;
    BranchBox* tmpbbox;
    for (int i=0;i<bboxes->size()-1;i++) {
        for (int j=i+1;j<bboxes->size();j++) {
            if ((*bboxes)[j]->xlim[0] < (*bboxes)[i]->xlim[0]) {
                tmpbbox      = (*bboxes)[j];
                (*bboxes)[j] = (*bboxes)[i];
                (*bboxes)[i] = tmpbbox;
            }
        }
    }
}

void PlecFinder::purgeInactiveBboxes(std::vector<BranchBox*> * bboxes) {
    if (bboxes->size()==0) return;
    for (int i=bboxes->size()-1;i>=0;i--) {
        if (!(*bboxes)[i]->is_active()) {
            bboxes->erase(bboxes->begin()+i);
            // delete (*bboxes)[i]; // leads to Segmentation Fault
        }
    }
}

void PlecFinder::purgeBboxesAtID(std::vector<BranchBox*> * bboxes,int id) {
    bboxes->erase(bboxes->begin()+id);
    // delete (*bboxes)[i]; // leads to Segmentation Fault
}


/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////


void PlecFinder::removeBboxOverlap(std::vector<BranchBox*> * bboxes) {

    if (bboxes->size()==0) return;
    for (int i=0;i<bboxes->size()-1;i++) {
        if (!(*bboxes)[i]->is_active()) continue;
        std::vector<int> xlimA, ylimA;
        xlimA = (*bboxes)[i]->xlim;
        ylimA = (*bboxes)[i]->ylim;

        for (int j=i+1;j<bboxes->size();j++) {
            if (!(*bboxes)[j]->is_active()) continue;
            if (!(*bboxes)[i]->is_active()) break;

            std::vector<int> xlimB,ylimB;
//            xlimA = (*bboxes)[i]->xlim;
//            ylimA = (*bboxes)[i]->ylim;

            xlimB = (*bboxes)[j]->xlim;
            ylimB = (*bboxes)[j]->ylim;

            std::vector<std::vector<int>> lims;
            lims = __balanceXOverlap(xlimA,ylimA,xlimB,ylimB);
            (*bboxes)[i]->set_xlim(lims[0]);
            (*bboxes)[j]->set_xlim(lims[1]);
            xlimA = lims[0];
            xlimB = lims[1];
            if ((!(*bboxes)[i]->is_active()) || (!(*bboxes)[j]->is_active())) continue;

            lims = __balanceXOverlap(ylimA,xlimA,ylimB,xlimB);
            (*bboxes)[i]->set_ylim(lims[0]);
            (*bboxes)[j]->set_ylim(lims[1]);
            ylimA = lims[0];
            ylimB = lims[1];
            if ((!(*bboxes)[i]->is_active()) || (!(*bboxes)[j]->is_active())) continue;

            lims = __balanceXOverlap(ylimA,xlimA,xlimB,ylimB);
            (*bboxes)[i]->set_ylim(lims[0]);
            (*bboxes)[j]->set_xlim(lims[1]);
            ylimA = lims[0];
            xlimB = lims[1];
            if ((!(*bboxes)[i]->is_active()) || (!(*bboxes)[j]->is_active())) continue;

            lims = __balanceXOverlap(xlimA,ylimA,ylimB,xlimB);
            (*bboxes)[i]->set_xlim(lims[0]);
            (*bboxes)[j]->set_ylim(lims[1]);
            xlimA = lims[0];
            ylimB = lims[1];
            // if ((!(*bboxes)[i]->is_active()) || (!(*bboxes)[j]->is_active())) continue;
        }
    }

    for (int i=0;i<bboxes->size();i++) {
        if ((*bboxes)[i]->get_writhe() < WRTHRESHOLD) {
            (*bboxes)[i]->active = false;
        }
    }
    purgeInactiveBboxes(bboxes);
    sortBboxes(bboxes);
}

std::vector<std::vector<int>> PlecFinder::__balanceXOverlap(const std::vector<int> & xlimA,const std::vector<int> & ylimA,const std::vector<int> & xlimB,const std::vector<int> & ylimB) {
    double totalwrLR,totalwrRL;
    std::vector<std::vector<int>> limsLR,limsRL;

    limsLR = __balanceXOverlapL2R(xlimA,ylimA,xlimB,ylimB,&totalwrLR);
    limsRL = __balanceXOverlapL2R(xlimB,ylimB,xlimA,ylimA,&totalwrRL);
    if (totalwrLR >= totalwrRL) {
        return limsLR;
    }
    else {
        return {limsRL[1],limsRL[0]};
    }
}

std::vector<std::vector<int>> PlecFinder::__balanceXOverlapL2R(const std::vector<int> & xlim1,const std::vector<int> & ylim1,const std::vector<int> & xlim2,const std::vector<int> & ylim2, double * totalwr) {

    // if there is no overlap
    if (xlim1[1] <= xlim2[0]) {
        *totalwr = 0;
        return {xlim1,xlim2};
    }
    if (xlim2[1] <= xlim1[0]) {
        *totalwr = 1e10;
        return {xlim1,xlim2};
    }

    int startid,endid;
    double box1_left_writhe,box1_right_writhe,box2_right_writhe;

    // find the starting point of the overlap region
    if (xlim1[0] < xlim2[0]) {
        startid = xlim2[0];
        box1_left_writhe = arma::accu(WM(arma::span(xlim1[0],startid-1),arma::span(ylim1[0],ylim1[1]-1)));
    }
    else {
        startid = xlim1[0];
        box1_left_writhe = 0;
    }

    // find the end point of the overlap region

    if ( xlim1[1] < xlim2[1] ) {
        endid = xlim1[1];
        box1_right_writhe = 0;
        box2_right_writhe = arma::accu(WM(arma::span(endid,xlim2[1]-1),arma::span(ylim2[0],ylim2[1]-1)));
    }
    else {
        endid = xlim2[1];
        box1_right_writhe = arma::accu(WM(arma::span(endid,xlim1[1]-1),arma::span(ylim1[0],ylim1[1]-1)));
        box2_right_writhe = 0;
    }

    // calculate the vertical sums along the overlap region
    arma::colvec mid1,mid2;
    mid1 = arma::sum(WM(arma::span(startid,endid-1),arma::span(ylim1[0],ylim1[1]-1)),1);
    mid2 = arma::sum(WM(arma::span(startid,endid-1),arma::span(ylim2[0],ylim2[1]-1)),1);
    mid1[mid1.n_elem-1] += box1_right_writhe;

    // search for the overlap cut that maximizes the writhe
    double midwr,maxwr;
    int    maxid;
    midwr = arma::sum(mid2);
    maxwr = midwr;
    maxid = 0;
    for (int i=0;i<mid1.n_elem;i++) {
        midwr = midwr + mid1(i) - mid2(i);
        if (midwr > maxwr ) {
            maxwr = midwr;
            maxid = i+1;
        }
    }

    std::vector<int> ret_xlim1,ret_xlim2;
    ret_xlim1 = xlim1;
    ret_xlim2 = xlim2;

    if ((maxid == mid1.n_elem) && (box1_right_writhe != 0)) {
        ret_xlim2[1] = ret_xlim2[0];
    }
    else {
        ret_xlim1[1] = startid + maxid;
        ret_xlim2[0] = ret_xlim1[1];
    }

    *totalwr = maxwr + box1_left_writhe + box2_right_writhe;
    return {ret_xlim1,ret_xlim2};
}



/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
// Remove comflicting bboxes
/*
    Comflicting bboxes are due to inter-plectoneme writhe contributions that are
    inconsistent with plectoneme downstream order
*/

void PlecFinder::simpleConflictRemoval(std::vector<BranchBox*> * bboxes) {
    if (bboxes->size()==0) return;
    arma::colvec conflicts = arma::zeros(bboxes->size());
    for (int i=0;i<bboxes->size()-1;i++) {
        for (int j=i+1;j<bboxes->size();j++) {
            if (!isDownstreamBranch((*bboxes)[i],(*bboxes)[j]) && (*bboxes)[i]->ylim[1] > (*bboxes)[j]->xlim[0] ) {
                double score = (*bboxes)[j]->get_writhe()/(*bboxes)[i]->get_writhe();
                if (score > 0) {
                    conflicts(i) += score;
                    conflicts(j) += 1./score;
                }
            }
        }
    }
    if (arma::accu(conflicts) > 1e-10) {
        int rem = arma_argmax(conflicts);
        purgeBboxesAtID(bboxes,rem);
        simpleConflictRemoval(bboxes);
    }
}

bool PlecFinder::isDownstreamBranch(BranchBox * refbbox, BranchBox * checkbbox) {
    if ((refbbox->xlim[1] < checkbbox->xlim[1]) && (refbbox->ylim[0] > checkbbox->ylim[0])) {
        return true;
    }
    return false;
}

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////


std::vector<int> PlecFinder::findEntryBboxes(const std::vector<BranchBox*> & bboxes) {
    /*
        returns the is of the bboxes at the root of plectonemes
    */
    std::vector<int> entry_bboxes_ids;
    int id   = 0;
    int last;
    while (id < bboxes.size()) {
        last = id;
        entry_bboxes_ids.push_back(last);
        id++;
        while ((id < bboxes.size()) && isDownstreamBranch(bboxes[last],bboxes[id]) ) {
            id++;
        }
    }
    return entry_bboxes_ids;
}



void PlecFinder::refineAllEntryBBoxes(const std::vector<BranchBox*> & bboxes,const std::vector<int> & entry_bboxes_ids, double max_remove) {
    for (unsigned i=0;i<entry_bboxes_ids.size();i++) {
        refineEntryBBox(bboxes[entry_bboxes_ids[i]],max_remove);
    }
}

void PlecFinder::refineEntryBBox(BranchBox * bbox,double max_remove) {
    double Wr,accu;
    int new_x,new_y;
    Wr    = bbox->get_writhe();
    /*
        refine xlims
    */
    new_x = bbox->xlim[0];
    accu  = 0;
    for (int i=bbox->xlim[0];i<bbox->xlim[1];i++) {
        accu += arma::accu(WM(i,arma::span(bbox->ylim[0],bbox->ylim[1]-1)));
        if (accu/Wr < max_remove) {
            new_x = i;
        }
        else break;
    }
    bbox->set_xlim({new_x,bbox->xlim[1]});

    /*
        refine ylims
    */
    new_y = bbox->ylim[1];
    accu  = 0;
    for (int i=bbox->ylim[1]-1;i>bbox->xlim[0];i--) {
        accu += arma::accu(WM(arma::span(bbox->xlim[0],bbox->xlim[1]-1),i));
        if (accu/Wr < max_remove) {
            new_y = i;
        }
        else break;
    }
    bbox->set_ylim({bbox->ylim[0],new_y});
    bbox->get_writhe();
}


/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////


std::vector<double> PlecFinder::getPlecStats() {

    double Wr,Tw,Lk,z;
    double pWr,pTw,pLk;
    int    num_plecs;
    int    plecsegs;
    double partWr;

//    for (unsigned i=0;i<WM.n_rows;i++) {
//        for (unsigned j=0;j<WM.n_cols;j++) {
//            if (WM(i,j) < 0) {
//                std::cout << "WM smaller than zero!" << std::endl;
//            }
//        }
//    }

    Wr = arma::accu(WM)*sgnWM;
    Tw = chain->cal_twist();
    Lk = Wr + Tw;
    z  = arma::dot(chain->get_force_dir(),chain->get_bp_pos()->col(num_bp-1)-chain->get_bp_pos()->col(0));

    int id1,id2;

    std::cout << "~~~~~~~~~~~~~~~~~~~~" << std::endl;
    pWr = 0; pTw = 0,plecsegs=0;
    for (unsigned i=0;i<entry_bboxes_ids.size();i++) {
        id1 = bboxes[entry_bboxes_ids[i]]->xlim[0];
        id2 = bboxes[entry_bboxes_ids[i]]->ylim[1];
        partWr    =   arma::accu(WM(arma::span(id1,id2-1),arma::span(id1,id2   -1)));
        partWr   += 2*arma::accu(WM(arma::span(0  ,id1-1),arma::span(id1,id2   -1))) ;
        partWr   += 2*arma::accu(WM(arma::span(id1,id2-1),arma::span(id2,WM.n_rows-1)));
        pWr      += partWr *sgnWM;
        pTw      += chain->cal_twist(id1,id2);
        plecsegs += id2-id1;

//        std::cout << pWr << " " << id1 << " " << id2 << std::endl;
        std::cout << id1 << " " << id2 << std::endl;
    }
    std::cout << "~~~~~~~~~~~~~~~~~~~~" << std::endl;
    pLk = pWr + pTw;
    num_plecs = entry_bboxes_ids.size();

    return {z, Lk, Wr, Tw, pLk, pTw, pWr, plecsegs, num_plecs};
}

