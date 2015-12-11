#ifndef ForestElectrons_h
#define ForestElectrons_h

#define maxForestElectrons 10

#include <iostream>
#include <vector>
#include "TBranch.h"

class ForestElectrons {
public :
   ForestElectrons(){};
   ~ForestElectrons(){};

   // Declaration of leaf types
   // Event info
   Int_t           Run;
   Int_t           Event;
   Int_t           Lumi;
   Int_t           CentBin;
   Float_t         vx;
   Float_t         vy;
   Float_t         vz;
   // GenParticle info
   Int_t           Gen_nptl;
   Int_t           Gen_pid[maxForestElectrons];   //[Gen_nptl]
   Int_t           Gen_mom[maxForestElectrons];   //[Gen_nptl] pid mother
   Int_t           Gen_status[maxForestElectrons];   //[Gen_nptl]
   Float_t         Gen_p[maxForestElectrons];   //[Gen_nptl]
   Float_t         Gen_pt[maxForestElectrons];   //[Gen_nptl]
   Float_t         Gen_eta[maxForestElectrons];   //[Gen_nptl]
   Float_t         Gen_phi[maxForestElectrons];   //[Gen_nptl]
   // RecoElectron info
   Int_t           Glb_nptl;
   std::vector<Int_t>           *Glb_charge=0;   //[Glb_nptl]
   std::vector<Float_t>         *Glb_p=0;   //[Glb_nptl]
   std::vector<Float_t>         *Glb_pt=0;   //[Glb_nptl]
   std::vector<Float_t>         *Glb_eta=0;   //[Glb_nptl]
   std::vector<Float_t>         *Glb_phi=0;   //[Glb_nptl]
   std::vector<Float_t>         *Glb_dxy=0;   //[Glb_nptl]
   std::vector<Float_t>         *Glb_dz=0;   //[Glb_nptl]
   std::vector<Int_t>           *Glb_nValMuHits=0;   //[Glb_nptl]
   std::vector<Int_t>           *Glb_nValTrkHits=0;   //[Glb_nptl]
   std::vector<Int_t>           *Glb_nValPixHits=0;   //[Glb_nptl]
   std::vector<Int_t>           *Glb_trkLayerWMeas=0;   //[Glb_nptl]
   std::vector<Int_t>           *Glb_nMatchedStations=0;   //[Glb_nptl]
   std::vector<Int_t>           *Glb_nTrkFound=0;   //[Glb_nptl]
   std::vector<Float_t>         *Glb_glbChi2_ndof=0;   //[Glb_nptl]
   std::vector<Float_t>         *Glb_trkChi2_ndof=0;   //[Glb_nptl]
   std::vector<Int_t>           *Glb_pixLayerWMeas=0;   //[Glb_nptl]
   std::vector<Float_t>         *Glb_trkDxy=0;   //[Glb_nptl]
   std::vector<Float_t>         *Glb_trkDz=0;   //[Glb_nptl]
   std::vector<Int_t>           *Glb_trkQuality=0; //[Glb_nptl]
   std::vector<Int_t>           *Glb_isGood=0; //[Glb_nptl]
   std::vector<Float_t>         *Glb_pfChIso=0; //[Glb_nptl]
   std::vector<Float_t>         *Glb_pfPhoIso=0; //[Glb_nptl]
   std::vector<Float_t>         *Glb_pfNeuIso=0; //[Glb_nptl]
   std::vector<Float_t>         *Glb_pfPUIso=0; //[Glb_nptl]
   std::vector<Int_t>           *Glb_isArbitrated=0;   //[Glb_nptl]
   //Needed ? 
   Int_t                        Sta_nptl;
   std::vector<Int_t>           *Sta_charge=0;   //[Sta_nptl]
   std::vector<Float_t>         *Sta_p=0;   //[Sta_nptl]
   std::vector<Float_t>         *Sta_pt=0;   //[Sta_nptl]
   std::vector<Float_t>         *Sta_eta=0;   //[Sta_nptl]
   std::vector<Float_t>         *Sta_phi=0;   //[Sta_nptl]
   std::vector<Float_t>         *Sta_dxy=0;   //[Sta_nptl]
   std::vector<Float_t>         *Sta_dz=0;   //[Sta_nptl]
   // Dimuon
   Int_t                        Di_npair;
   std::vector<Float_t>         *Di_vProb=0;   //[Di_npair]
   std::vector<Float_t>         *Di_mass=0;   //[Di_npair]
   std::vector<Float_t>         *Di_e=0;   //[Di_npair]
   std::vector<Float_t>         *Di_pt=0;   //[Di_npair]
   std::vector<Float_t>         *Di_pt1=0;   //[Di_npair]
   std::vector<Float_t>         *Di_pt2=0;   //[Di_npair]
   std::vector<Float_t>         *Di_eta=0;   //[Di_npair]
   std::vector<Float_t>         *Di_eta1=0;   //[Di_npair]
   std::vector<Float_t>         *Di_eta2=0;   //[Di_npair]
   std::vector<Float_t>         *Di_phi=0;   //[Di_npair]
   std::vector<Float_t>         *Di_phi1=0;   //[Di_npair]
   std::vector<Float_t>         *Di_phi2=0;   //[Di_npair]
   std::vector<Int_t>           *Di_charge1=0;   //[Di_npair]
   std::vector<Int_t>           *Di_charge2=0;   //[Di_npair]
   std::vector<Int_t>           *Di_isArb1=0;   //[Di_npair]
   std::vector<Int_t>           *Di_isArb2=0;   //[Di_npair]
   std::vector<Float_t>         *Di_nTrkHit1=0;   //[Di_npair]
   std::vector<Float_t>         *Di_nTrkHit2=0;   //[Di_npair]
   std::vector<Float_t>         *Di_trkChi2_1=0;   //[Di_npair]
   std::vector<Float_t>         *Di_trkChi2_2=0;   //[Di_npair]
   std::vector<Float_t>         *Di_glbChi2_1=0;   //[Di_npair]
   std::vector<Float_t>         *Di_glbChi2_2=0;   //[Di_npair]
   std::vector<Float_t>         *Di_dxy1=0;   //[Di_npair]
   std::vector<Float_t>         *Di_dxy2=0;   //[Di_npair]
   std::vector<Float_t>         *Di_dz1=0;   //[Di_npair]
   std::vector<Float_t>         *Di_dz2=0;   //[Di_npair]

   // List of branches
   // Event info
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_cbin;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!
   // GenParticle info
   TBranch        *b_Gen_nptl;   //!
   TBranch        *b_Gen_pid;   //!
   TBranch        *b_Gen_mom;   //!
   TBranch        *b_Gen_status;   //!
   TBranch        *b_Gen_p;   //!
   TBranch        *b_Gen_pt;   //!
   TBranch        *b_Gen_eta;   //!
   TBranch        *b_Gen_phi;   //!
   // RecoElectron info
   TBranch        *b_Glb_nptl;   //!
   TBranch        *b_Glb_charge;   //!
   TBranch        *b_Glb_p;   //!
   TBranch        *b_Glb_pt;   //!
   TBranch        *b_Glb_eta;   //!
   TBranch        *b_Glb_phi;   //!
   TBranch        *b_Glb_dxy;   //!
   TBranch        *b_Glb_dz;   //!
   TBranch        *b_Glb_nValMuHits;   //!
   TBranch        *b_Glb_nValTrkHits;   //!
   TBranch        *b_Glb_nValPixHits;   //!
   TBranch        *b_Glb_trkLayerWMeas;   //!
   TBranch        *b_Glb_nMatchedStations;   //!
   TBranch        *b_Glb_nTrkFound;   //!
   TBranch        *b_Glb_glbChi2_ndof;   //!
   TBranch        *b_Glb_trkChi2_ndof;   //!
   TBranch        *b_Glb_pixLayerWMeas;   //!
   TBranch        *b_Glb_trkDxy;   //!
   TBranch        *b_Glb_trkDz;   //!
   TBranch        *b_Glb_trkQuality;
   TBranch        *b_Glb_isGood;
   TBranch        *b_Glb_pfChIso;
   TBranch        *b_Glb_pfPhoIso;
   TBranch        *b_Glb_pfNeuIso;
   TBranch        *b_Glb_pfPUIso;
   TBranch        *b_Glb_isArbitrated;   //!
   // Needed?
   TBranch        *b_Sta_nptl;   //!
   TBranch        *b_Sta_charge;   //!
   TBranch        *b_Sta_p;   //!
   TBranch        *b_Sta_pt;   //!
   TBranch        *b_Sta_eta;   //!
   TBranch        *b_Sta_phi;   //!
   TBranch        *b_Sta_dxy;   //!
   TBranch        *b_Sta_dz;   //!
   //Dimuons
   TBranch        *b_Di_npair;   //!
   TBranch        *b_Di_vProb;   //!
   TBranch        *b_Di_mass;   //!
   TBranch        *b_Di_e;   //!
   TBranch        *b_Di_pt;   //!
   TBranch        *b_Di_pt1;   //!
   TBranch        *b_Di_pt2;   //!
   TBranch        *b_Di_eta;   //!
   TBranch        *b_Di_eta1;   //!
   TBranch        *b_Di_eta2;   //!
   TBranch        *b_Di_phi;   //!
   TBranch        *b_Di_phi1;   //!
   TBranch        *b_Di_phi2;   //!
   TBranch        *b_Di_charge1;   //!
   TBranch        *b_Di_charge2;   //!
   TBranch        *b_Di_isArb1;   //!
   TBranch        *b_Di_isArb2;   //!
   TBranch        *b_Di_nTrkHit1;   //!
   TBranch        *b_Di_nTrkHit2;   //!
   TBranch        *b_Di_trkChi2_1;   //!
   TBranch        *b_Di_trkChi2_2;   //!
   TBranch        *b_Di_glbChi2_1;   //!
   TBranch        *b_Di_glbChi2_2;   //!
   TBranch        *b_Di_dxy1;   //!
   TBranch        *b_Di_dxy2;   //!
   TBranch        *b_Di_dz1;   //!
   TBranch        *b_Di_dz2;   //!
};
#endif
